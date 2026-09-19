#!/usr/bin/env bash
# spades-rs: Fast, zero-dependency installer for Linux, macOS, and Google Colab
set -euo pipefail

REPO="sagnikrout/spades-rs"
LATEST_URL="https://github.com/${REPO}/releases/latest/download"

# 1. Detect OS & Architecture
OS="$(uname -s | tr '[:upper:]' '[:lower:]')"
ARCH="$(uname -m)"

case "${OS}" in
  linux)
    case "${ARCH}" in
      x86_64|amd64)   ASSET="spades-rs-linux-x86_64-musl" ;;
      aarch64|arm64)  ASSET="spades-rs-linux-arm64" ;;
      *)
        echo "Error: Unsupported architecture: ${ARCH}. spades-rs requires a 64-bit system (x86_64 or aarch64)." >&2
        exit 1
        ;;
    esac
    ;;
  darwin)
    case "${ARCH}" in
      arm64|aarch64)  ASSET="spades-rs-macos-arm64" ;;
      *)
        echo "Error: Intel macOS (${ARCH}) is not supported with prebuilt binaries. Please build with 'cargo build --release'." >&2
        exit 1
        ;;
    esac
    ;;
  *)
    echo "Error: Unsupported OS: ${OS}. Supported: Linux, macOS (Apple Silicon)." >&2
    exit 1
    ;;
esac

# 2. Determine target installation directory
if [ -w "/usr/local/bin" ]; then
  DEST_DIR="/usr/local/bin"
elif [ -n "${PREFIX:-}" ] && [ -w "${PREFIX}" ]; then
  DEST_DIR="${PREFIX}"
else
  DEST_DIR="${HOME}/.local/bin"
  mkdir -p "${DEST_DIR}"
  if [[ ":${PATH}:" != *":${DEST_DIR}:"* ]]; then
    export PATH="${DEST_DIR}:${PATH}"
  fi
fi

BIN_PATH="${DEST_DIR}/spades-rs"
ALIAS_PATH="${DEST_DIR}/spades.py"

echo "Installing spades-rs (${ASSET}) into ${DEST_DIR}..."

# 3. Download release binary
DOWNLOAD_URL="${LATEST_URL}/${ASSET}"
if command -v curl >/dev/null 2>&1; then
  curl -fsSL "${DOWNLOAD_URL}" -o "${BIN_PATH}"
elif command -v wget >/dev/null 2>&1; then
  wget -qO "${BIN_PATH}" "${DOWNLOAD_URL}"
else
  echo "Error: Neither curl nor wget found." >&2
  exit 1
fi

chmod +x "${BIN_PATH}"
ln -sf "${BIN_PATH}" "${ALIAS_PATH}" 2>/dev/null || cp -f "${BIN_PATH}" "${ALIAS_PATH}"

echo "Successfully installed spades-rs to ${BIN_PATH}"
echo "Compatibility alias created: ${ALIAS_PATH}"
echo ""
"${BIN_PATH}" --help | head -n 8
echo ""
echo "Done! You can now run 'spades-rs assemble -h' or 'spades.py assemble -h'."
