"""Registry of external design tools.

Adapters are discovered by name. Nothing here is imported as a library: each tool
is a subprocess with its own environment, so their dependency pins stay isolated
from this package and from each other.
"""
from __future__ import annotations

from .base import Tool, ToolEnv, ToolError
from .siren import SirenTool
from .sirnadesigner import SirnaDesignerTool
from .sirnaforge import SirnaforgeTool

REGISTRY: dict[str, type[Tool]] = {
    SirnaDesignerTool.name: SirnaDesignerTool,
    SirenTool.name: SirenTool,
    SirnaforgeTool.name: SirnaforgeTool,
}

DEFAULT_TOOLS = ("sirnadesigner", "siren", "sirnaforge")


def available_tools() -> list[str]:
    return sorted(REGISTRY)


def get_tool(name: str, env: ToolEnv | None = None) -> Tool:
    try:
        return REGISTRY[name](env)
    except KeyError:
        raise ToolError(
            f"unknown tool {name!r}; known: {', '.join(available_tools())}") from None


def survey(env: ToolEnv | None = None) -> list[tuple[str, bool, str]]:
    """Availability of every registered tool, without running any of them."""
    out = []
    for name in available_tools():
        usable, detail = get_tool(name, env).availability()
        out.append((name, usable, detail))
    return out


__all__ = [
    "DEFAULT_TOOLS",
    "REGISTRY",
    "Tool",
    "ToolEnv",
    "ToolError",
    "available_tools",
    "get_tool",
    "survey",
]
