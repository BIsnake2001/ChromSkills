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
# Add DeepSeek Anthropic-compatible configuration to ~/.bashrc if not already present
CONFIG_LINE1='export ANTHROPIC_BASE_URL="https://api.deepseek.com/anthropic"'
CONFIG_LINE2='export ANTHROPIC_MODEL="deepseek-v4-pro[1m]"'
CONFIG_LINE3='export ANTHROPIC_SMALL_FAST_MODEL="deepseek-v4-flash"'

if ! grep -Fxq "$CONFIG_LINE1" ~/.bashrc; then
  echo "$CONFIG_LINE1" >> ~/.bashrc
fi
if ! grep -Fxq "$CONFIG_LINE2" ~/.bashrc; then
  echo "$CONFIG_LINE2" >> ~/.bashrc
fi
if ! grep -Fxq "$CONFIG_LINE3" ~/.bashrc; then
  echo "$CONFIG_LINE3" >> ~/.bashrc
fi

# Reload ~/.bashrc for current shell session if running interactively
DATA_DIR="$HOME/ChromOmics/ChIPseq"
mkdir -p "$DATA_DIR"


echo "DeepSeek environment variables configured."
echo ""


# ------------------------------------------------------------
# Download review data
# ------------------------------------------------------------


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

# download_file "https://zenodo.org/record/1324070/files/wt_H3K4me3_rep1.bam"
# download_file "https://zenodo.org/record/1324070/files/wt_H3K4me3_rep2.bam"
# download_file "https://zenodo.org/record/1324070/files/wt_H3K27me3_rep1.bam"
# download_file "https://zenodo.org/record/1324070/files/wt_H3K27me3_rep2.bam"
# download_file "https://zenodo.org/record/1324070/files/wt_input_rep1.bam"
# download_file "https://zenodo.org/record/1324070/files/wt_input_rep2.bam"

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
cd "$DATA_DIR"

if ! command -v claude >/dev/null 2>&1; then
  echo "ERROR: claude command was not found."
  exit 1
fi

echo ""
echo "then run the following steps in a new terminal window ..."
echo ""
echo "cd ~/ChromOmics/ChIPseq"
echo "source ~/.bashrc"
echo "export ANTHROPIC_API_KEY=\"your_api_key\""
echo "claude /mcp"


