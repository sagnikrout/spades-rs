use intelligent_pascal::fastq::parse_reads_from_file;
use intelligent_pascal::gfa::write_graph_gfa;
use intelligent_pascal::graph::Unitig;
use std::fs::File;
use std::io::{Read, Write};

#[test]
fn test_gfa_export_format() {
    let temp_dir = std::env::temp_dir();
    let gfa_path = temp_dir.join("test_output.gfa");

    let k = 5;
    let unitigs = vec![
        Unitig {
            id: 0,
            sequence: b"GGCTACGT".to_vec(),
            mean_coverage: 15.5,
            kmers_count: 4,
        },
        Unitig {
            id: 1,
            sequence: b"ACGTCCTA".to_vec(),
            mean_coverage: 25.0,
            kmers_count: 4,
        },
    ];

    write_graph_gfa(k, &unitigs, &gfa_path).expect("GFA export failed");

    let mut content = String::new();
    let mut file = File::open(&gfa_path).expect("Failed to read GFA");
    file.read_to_string(&mut content).unwrap();

    let lines: Vec<&str> = content.lines().collect();

    // Line 1: Header
    assert_eq!(lines[0], "H\tVN:Z:1.0");

    // S lines: Segments
    assert!(lines.contains(&"S\t0\tGGCTACGT\tLN:i:8\tRC:f:15.5"));
    assert!(lines.contains(&"S\t1\tACGTCCTA\tLN:i:8\tRC:f:25.0"));

    // L lines: Overlap of k-1 = 4 bases (ACGT)
    assert!(lines.contains(&"L\t0\t+\t1\t+\t4M"));

    let _ = std::fs::remove_file(gfa_path);
}

#[test]
fn test_fastq_parsing_plain() {
    let temp_dir = std::env::temp_dir();
    let fq_path = temp_dir.join("test_reads.fq");

    let fastq_content = b"@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTTGGGG\n+\nHHHHHHHH\n@read3\nCCCCAAAA\n+\nJJJJJJJJ\n";
    {
        let mut f = File::create(&fq_path).unwrap();
        f.write_all(fastq_content).unwrap();
    }

    let parsed = parse_reads_from_file(&fq_path).expect("FASTQ parse failed");
    assert_eq!(parsed.len(), 3);
    assert_eq!(parsed[0], b"ACGTACGT");
    assert_eq!(parsed[1], b"TTTTGGGG");
    assert_eq!(parsed[2], b"CCCCAAAA");

    let _ = std::fs::remove_file(fq_path);
}

#[test]
fn test_fasta_multiline_parsing() {
    let temp_dir = std::env::temp_dir();
    let fa_path = temp_dir.join("test_reads.fa");

    let fasta_content = b">seq1\nACGT\nACGT\n>seq2\nGGGG\nTTTT\nAAAA\n";
    {
        let mut f = File::create(&fa_path).unwrap();
        f.write_all(fasta_content).unwrap();
    }

    let parsed = parse_reads_from_file(&fa_path).expect("FASTA parse failed");
    assert_eq!(parsed.len(), 2);
    assert_eq!(parsed[0], b"ACGTACGT");
    assert_eq!(parsed[1], b"GGGGTTTTAAAA");

    let _ = std::fs::remove_file(fa_path);
}

#[test]
fn test_empty_file_parsing() {
    let temp_dir = std::env::temp_dir();
    let empty_path = temp_dir.join("test_empty.fa");
    File::create(&empty_path).unwrap();

    let parsed = parse_reads_from_file(&empty_path).expect("Empty file parse failed");
    assert!(parsed.is_empty());

    let _ = std::fs::remove_file(empty_path);
}
