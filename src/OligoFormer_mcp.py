"""OligoFormer MCP server: siRNA efficacy, off-target/toxicity screening and mismatch tools."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from fastmcp import FastMCP

from tools.cli_wrapper import cli_wrapper_mcp

mcp = FastMCP(name="OligoFormer")
mcp.mount(cli_wrapper_mcp)

if __name__ == "__main__":
    mcp.run()
