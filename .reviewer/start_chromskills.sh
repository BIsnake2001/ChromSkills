#!/usr/bin/env bash
set -euo pipefail

echo ""
echo "=================================================="
echo " ChromSkills reviewer environment"
echo "=================================================="
echo ""
echo "You are now inside the ChromSkills Docker container."
echo ""
echo "Repository is mounted at:"
echo "  /repo"
echo ""
echo "Working directory is:"
echo "  /work"
echo ""

cd /work

echo "Installing / updating ChromSkills MCP tools and skills..."
if [ -x "$HOME/scripts/install-claude.sh" ]; then
  cd "$HOME/scripts"
  ./install-claude.sh
else
  echo "Warning: $HOME/scripts/install-claude.sh was not found or not executable."
fi

cd /work

echo ""
echo "=================================================="
echo " Coding agent configuration"
echo "=================================================="
echo ""
echo "For Anthropic Claude backend, paste your API key below."
echo "The input will be hidden."
echo ""
echo "If you prefer DeepSeek / MiniMax / another Anthropic-compatible"
echo "backend, press Enter here and manually export:"
echo ""
echo "  export ANTHROPIC_BASE_URL=..."
echo "  export ANTHROPIC_AUTH_TOKEN=..."
echo ""
echo "Then run:"
echo ""
echo "  claude /mcp"
echo "  chat"
echo ""

if [ -z "${ANTHROPIC_API_KEY:-}" ] && [ -z "${ANTHROPIC_AUTH_TOKEN:-}" ]; then
  read -r -s -p "ANTHROPIC_API_KEY, or press Enter to skip: " USER_ANTHROPIC_API_KEY
  echo ""
  if [ -n "$USER_ANTHROPIC_API_KEY" ]; then
    export ANTHROPIC_API_KEY="$USER_ANTHROPIC_API_KEY"
  fi
fi

echo ""
echo "Current directory:"
pwd

echo ""
echo "Available files:"
ls -lah

echo ""
echo "=================================================="
echo " Next step"
echo "=================================================="
echo ""

if [ -n "${ANTHROPIC_API_KEY:-}" ] || [ -n "${ANTHROPIC_AUTH_TOKEN:-}" ]; then
  echo "Initializing Claude MCP..."
  claude /mcp || true

  echo ""
  echo "Now start ChromSkills by running:"
  echo ""
  echo "  chat"
else
  echo "No API key was configured."
  echo ""
  echo "Configure your backend first, for example:"
  echo ""
  echo "  export ANTHROPIC_API_KEY=your_key_here"
  echo "  claude /mcp"
  echo "  chat"
fi

echo ""
exec /bin/bash