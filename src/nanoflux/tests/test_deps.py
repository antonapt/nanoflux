import shutil


def find_tool(tool_name: str) -> bool:
    return shutil.which(tool_name) is not None


def test_tools() -> None:
    tools = ["minimap2", "samtools", "modkit", "bedtools"]

    for tool in tools:
        assert find_tool(tool), f"Tool is not installed or not in PATH: {tool}"
