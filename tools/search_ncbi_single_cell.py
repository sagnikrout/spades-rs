import urllib.request
import xml.etree.ElementTree as ET
import json

# Search SRA for single cell E. coli MDA Illumina paired
query = '("Escherichia coli"[Organism] AND ("single cell"[All Fields] OR "single-cell"[All Fields]) AND "MDA"[All Fields] AND "Illumina"[Platform] AND "paired"[Layout])'
url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={urllib.parse.quote(query)}&retmode=json&retmax=20"

print(f"Querying NCBI: {url}")
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
with urllib.request.urlopen(req) as response:
    data = json.loads(response.read().decode())

id_list = data['esearchresult']['idlist']
print(f"Found {len(id_list)} SRA IDs: {id_list}")

if id_list:
    summary_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id={','.join(id_list)}&retmode=json"
    req2 = urllib.request.Request(summary_url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req2) as resp2:
        sum_data = json.loads(resp2.read().decode())
    
    for uid in id_list:
        item = sum_data['result'][uid]
        exp = item.get('expxml', '')
        runs = item.get('runs', '')
        print(f"\nUID: {uid}")
        print(f"Runs: {runs}")
        # Parse expxml for title and platform
        root = ET.fromstring(f"<root>{exp}</root>")
        title = root.findtext('.//Title')
        bioproject = root.findtext('.//Bioproject')
        print(f"Title: {title}")
        print(f"BioProject: {bioproject}")
