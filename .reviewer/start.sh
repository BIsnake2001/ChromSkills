#!/usr/bin/env bash
set -Eeuo pipefail

echo ""
echo "============================================================"
echo " ChromSkills Reviewer Setup"
echo "============================================================"
echo ""

# Go to repository root
if [ -d "/workspaces/ChromSkills" ]; then
  cd /workspaces/ChromSkills
elif [ -d "/workspace/ChromSkills" ]; then
  cd /workspace/ChromSkills
elif [ -d "/repo" ]; then
  cd /repo
else
  cd "$(pwd)"
fi

echo "Current directory:"
pwd
echo ""

# ------------------------------------------------------------
# DeepSeek Anthropic-compatible configuration
# ------------------------------------------------------------
export ANTHROPIC_BASE_URL="https://api.deepseek.com/anthropic"
export ANTHROPIC_MODEL="deepseek-v4-pro[1m]"
export ANTHROPIC_SMALL_FAST_MODEL="deepseek-v4-flash"

# Optional but useful for Claude Code model routing
export ANTHROPIC_DEFAULT_OPUS_MODEL="deepseek-v4-pro[1m]"
export ANTHROPIC_DEFAULT_SONNET_MODEL="deepseek-v4-pro[1m]"
export ANTHROPIC_DEFAULT_HAIKU_MODEL="deepseek-v4-flash"
export CLAUDE_CODE_SUBAGENT_MODEL="deepseek-v4-flash"

echo "DeepSeek environment variables configured."
echo ""

# ------------------------------------------------------------
# Ask reviewer for API key
# ------------------------------------------------------------
if [ -z "${ANTHROPIC_AUTH_TOKEN:-}" ]; then
  echo "Please paste your DeepSeek API key."
  echo "Input will be hidden."
  echo ""

  read -r -s -p "DEEPSEEK_API_KEY: " REVIEWER_DEEPSEEK_API_KEY
  echo ""

  if [ -z "${REVIEWER_DEEPSEEK_API_KEY:-}" ]; then
    echo "ERROR: No API key was provided."
    exit 1
  fi

  export ANTHROPIC_AUTH_TOKEN="${REVIEWER_DEEPSEEK_API_KEY}"
fi

echo "API key configured for this session."
echo ""

# ------------------------------------------------------------
# Download review data
# ------------------------------------------------------------
DATA_DIR="${CHROMSKILLS_DATA_DIR:-$PWD/data/zenodo_1324070}"
mkdir -p "$DATA_DIR"

echo "============================================================"
echo " Downloading review data"
echo "============================================================"
echo ""
echo "Data directory:"
echo "  $DATA_DIR"
echo ""

download_file() {
  local url="$1"
  local filename
  filename="$(basename "$url")"
  local output_path="$DATA_DIR/$filename"

  if [ -s "$output_path" ]; then
    echo "Already exists, skipping: $filename"
    return 0
  fi

  echo "Downloading: $filename"

  if command -v wget >/dev/null 2>&1; then
    wget -c -O "$output_path" "$url"
  elif command -v curl >/dev/null 2>&1; then
    curl -L --fail --continue-at - -o "$output_path" "$url"
  else
    echo "ERROR: neither wget nor curl is available."
    exit 1
  fi
}

download_file "https://zenodo.org/record/1324070/files/wt_H3K4me3_rep1.bam"
download_file "https://zenodo.org/record/1324070/files/wt_H3K4me3_rep2.bam"
download_file "https://zenodo.org/record/1324070/files/wt_H3K27me3_rep1.bam"
download_file "https://zenodo.org/record/1324070/files/wt_H3K27me3_rep2.bam"
download_file "https://zenodo.org/record/1324070/files/wt_input_rep1.bam"
download_file "https://zenodo.org/record/1324070/files/wt_input_rep2.bam"

export CHROMSKILLS_DATA_DIR="$DATA_DIR"
export DATA_DIR="$DATA_DIR"

echo ""
echo "Data download finished."
echo "CHROMSKILLS_DATA_DIR=$CHROMSKILLS_DATA_DIR"
echo ""

# ------------------------------------------------------------
# Run install-claude.sh
# ------------------------------------------------------------
echo "============================================================"
echo " Installing ChromSkills Claude tools"
echo "============================================================"
echo ""

INSTALLER=""

if [ -f "$HOME/scripts/install-claude.sh" ]; then
  INSTALLER="$HOME/scripts/install-claude.sh"
elif [ -f "/root/scripts/install-claude.sh" ]; then
  INSTALLER="/root/scripts/install-claude.sh"
fi

if [ -z "$INSTALLER" ]; then
  echo "ERROR: Cannot find install-claude.sh."
  echo ""
  echo "Expected:"
  echo "  ~/scripts/install-claude.sh"
  echo "or:"
  echo "  /root/scripts/install-claude.sh"
  echo ""
  echo "Please check whether the Codespace is using the ChromSkills Docker image."
  exit 1
fi

chmod +x "$INSTALLER"

(
  cd "$(dirname "$INSTALLER")"
  ./install-claude.sh
)

echo ""
echo "install-claude.sh finished."
echo ""

# ------------------------------------------------------------
# Start Claude / MCP / chat
# ------------------------------------------------------------
echo "============================================================"
echo " Starting Claude"
echo "============================================================"
echo ""

if ! command -v claude >/dev/null 2>&1; then
  echo "ERROR: claude command was not found."
  exit 1
fi

echo "Initializing MCP..."
claude /mcp || true

echo ""
echo "Opening ChromSkills chat..."
echo ""

if command -v chat >/dev/null 2>&1; then
  exec chat
else
  exec claude
fi