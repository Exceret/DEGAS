"""
R pkg cli-like message print with timestamp
"""

import sys
from datetime import datetime
from typing import Literal

#: symbol -> (display character, ANSI color code) table, built once at import time
_SYMBOLS: dict[str, tuple[str, str]] = {
    "info": ("ℹ", "\033[36m"),  # Cyan
    "success": ("✔", "\033[32m"),  # Green
    "warning": ("⚠", "\033[33m"),  # Yellow
    "error": ("✖", "\033[31m"),  # Red
    "debug": ("◼", "\033[35m"),  # Magenta
}


def ts_print(
    message: str,
    symbol: Literal["info", "success", "warning", "error", "debug"] | None = None,
    color: bool | None = True,
) -> None:
    """
    Print with Timestamp

    Args:
        message: A message to print
        symbol: CLI symbol: 'info', 'success', 'warning', 'error', 'debug'. Default: None
        color: Whether to print with color. Default: True
    """
    # 未显式指定时自动检测：重定向到文件/管道时关闭颜色，避免日志中出现转义序列。
    if color is None:
        color = sys.stdout.isatty()

    # 无效 symbol 直接报错，而不是静默退化为无符号输出（避免拼写错误难以察觉）。
    if symbol is not None and symbol not in _SYMBOLS:
        raise ValueError(f"Unsupported symbol: {symbol!r}")

    timestamp = datetime.now().astimezone().strftime("[%Y/%m/%d %H:%M:%S]")

    if symbol is not None:
        sym, col = _SYMBOLS[symbol]
        symbol_prefix = f"{col}{sym} \033[0m" if color else f"{sym} "
    else:
        symbol_prefix = ""

    print(symbol_prefix + f"{timestamp} {message}")
