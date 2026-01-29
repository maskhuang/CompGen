"""
CompGene Dual Format Logging.

This module provides a dual-format logger that outputs both
human-readable text logs and machine-readable JSON Lines logs.

Source: ADR-004 Log Format (Dual Output)
"""

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from workflow.lib.io import atomic_append


# =============================================================================
# ANSI Color Codes
# =============================================================================

class ANSIColors:
    """ANSI color codes for terminal output."""
    RED = "\033[91m"
    YELLOW = "\033[93m"
    GREEN = "\033[92m"
    CYAN = "\033[96m"
    GRAY = "\033[90m"
    RESET = "\033[0m"
    BOLD = "\033[1m"


# Level to color mapping
LEVEL_COLORS = {
    logging.DEBUG: ANSIColors.GRAY,
    logging.INFO: ANSIColors.GREEN,
    logging.WARNING: ANSIColors.YELLOW,
    logging.ERROR: ANSIColors.RED,
}


# =============================================================================
# Log Path Generation
# =============================================================================

def get_log_paths(
    rule: str,
    wildcards: dict[str, str],
    log_dir: Path = Path("logs")
) -> tuple[Path, Path]:
    """
    Generate log file paths for a rule execution.

    Args:
        rule: The rule name.
        wildcards: Dictionary of wildcard values.
        log_dir: Base directory for logs (default: "logs").

    Returns:
        Tuple of (text_log_path, jsonl_log_path).

    Example:
        >>> log_path, jsonl_path = get_log_paths("orthofinder", {"species_set": "lemur_macaque"})
        >>> print(log_path)
        logs/orthofinder/species_set=lemur_macaque.log
    """
    # Build wildcard string for filename
    if wildcards:
        wildcard_str = "_".join(f"{k}={v}" for k, v in sorted(wildcards.items()))
    else:
        wildcard_str = "default"

    base_path = log_dir / rule / wildcard_str

    return (
        base_path.with_suffix(".log"),
        base_path.with_suffix(".jsonl")
    )


def format_wildcards_for_path(wildcards: dict[str, str]) -> str:
    """
    Format wildcards dictionary into a path-safe string.

    Args:
        wildcards: Dictionary of wildcard values.

    Returns:
        Path-safe string representation.
    """
    if not wildcards:
        return "default"
    return "_".join(f"{k}={v}" for k, v in sorted(wildcards.items()))


# =============================================================================
# DualLogger Class
# =============================================================================

class DualLogger:
    """
    Logger that outputs both human-readable text and JSON Lines format.

    Provides simultaneous output to:
    - .log file: Human-readable text with timestamps and optional colors
    - .jsonl file: Machine-readable JSON Lines for jq/pandas querying

    Attributes:
        rule: The rule name being logged.
        wildcards: Dictionary of wildcard values.
        log_path: Path to the text log file.
        jsonl_path: Path to the JSON Lines log file.
        level: Current logging level.
        use_color: Whether to use ANSI colors in text output.

    Example:
        >>> logger = DualLogger("orthofinder", {"species_set": "lemur"})
        >>> logger.info("Starting analysis", species_count=3)
        >>> logger.error("Analysis failed", error_code="E_TIMEOUT")
    """

    def __init__(
        self,
        rule: str,
        wildcards: Optional[dict[str, str]] = None,
        log_dir: Path = Path("logs"),
        level: int = logging.INFO,
        use_color: bool = True
    ):
        """
        Initialize a DualLogger.

        Args:
            rule: The rule name being logged.
            wildcards: Dictionary of wildcard values (default: empty).
            log_dir: Base directory for logs (default: "logs").
            level: Logging level (default: INFO).
            use_color: Whether to use ANSI colors (default: True).
        """
        self.rule = rule
        self.wildcards = wildcards or {}
        self.level = level
        self.use_color = use_color

        # Generate log paths
        self.log_path, self.jsonl_path = get_log_paths(rule, self.wildcards, log_dir)

        # Ensure log directories exist
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

    def _should_log(self, level: int) -> bool:
        """Check if a message at the given level should be logged."""
        return level >= self.level

    def _format_text_line(self, level: int, msg: str) -> str:
        """
        Format a log message for text output.

        Format: YYYY-MM-DD HH:MM:SS [LEVEL] rule: message
        """
        now = datetime.now()
        timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
        level_name = logging.getLevelName(level)

        if self.use_color:
            color = LEVEL_COLORS.get(level, "")
            reset = ANSIColors.RESET
            return f"{timestamp} {color}[{level_name}]{reset} {self.rule}: {msg}\n"
        else:
            return f"{timestamp} [{level_name}] {self.rule}: {msg}\n"

    def _format_json_line(self, level: int, msg: str, **extra: Any) -> str:
        """
        Format a log message for JSON Lines output.

        Format: {"ts": "ISO8601", "level": "LEVEL", "rule": "rule", "msg": "message", ...}
        """
        now = datetime.now(timezone.utc)
        record = {
            "ts": now.isoformat(),
            "level": logging.getLevelName(level),
            "rule": self.rule,
            "msg": msg,
        }

        # Add wildcards if present
        if self.wildcards:
            record["wildcards"] = self.wildcards

        # Add extra fields
        record.update(extra)

        return json.dumps(record, ensure_ascii=False) + "\n"

    def _log(self, level: int, msg: str, **extra: Any) -> None:
        """
        Internal log method that writes to both outputs.

        Args:
            level: Logging level.
            msg: Log message.
            **extra: Additional fields for JSON output.
        """
        if not self._should_log(level):
            return

        # Write to text log
        text_line = self._format_text_line(level, msg)
        atomic_append(self.log_path, text_line)

        # Write to JSON Lines log
        json_line = self._format_json_line(level, msg, **extra)
        atomic_append(self.jsonl_path, json_line)

    def debug(self, msg: str, **extra: Any) -> None:
        """Log a DEBUG level message."""
        self._log(logging.DEBUG, msg, **extra)

    def info(self, msg: str, **extra: Any) -> None:
        """Log an INFO level message."""
        self._log(logging.INFO, msg, **extra)

    def warning(self, msg: str, **extra: Any) -> None:
        """Log a WARNING level message."""
        self._log(logging.WARNING, msg, **extra)

    def error(self, msg: str, **extra: Any) -> None:
        """Log an ERROR level message."""
        self._log(logging.ERROR, msg, **extra)

    def set_level(self, level: int) -> None:
        """
        Set the logging level.

        Args:
            level: New logging level (e.g., logging.DEBUG).
        """
        self.level = level

    def get_level_from_string(self, level_str: str) -> int:
        """
        Convert a level string to logging level constant.

        Args:
            level_str: Level name (DEBUG, INFO, WARNING, ERROR).

        Returns:
            Logging level constant.
        """
        level_map = {
            "DEBUG": logging.DEBUG,
            "INFO": logging.INFO,
            "WARNING": logging.WARNING,
            "ERROR": logging.ERROR,
        }
        return level_map.get(level_str.upper(), logging.INFO)


# =============================================================================
# Convenience Functions
# =============================================================================

def create_logger(
    rule: str,
    wildcards: Optional[dict[str, str]] = None,
    log_dir: Optional[Path] = None,
    level: str = "INFO",
    use_color: bool = True
) -> DualLogger:
    """
    Create a DualLogger with common defaults.

    Args:
        rule: The rule name being logged.
        wildcards: Dictionary of wildcard values.
        log_dir: Base directory for logs (default: "logs").
        level: Logging level as string (default: "INFO").
        use_color: Whether to use ANSI colors (default: True).

    Returns:
        Configured DualLogger instance.

    Example:
        >>> logger = create_logger("orthofinder", {"species_set": "lemur"})
        >>> logger.info("Starting...")
    """
    if log_dir is None:
        log_dir = Path("logs")

    logger = DualLogger(
        rule=rule,
        wildcards=wildcards,
        log_dir=log_dir,
        level=logging.INFO,
        use_color=use_color
    )

    # Set level from string
    logger.set_level(logger.get_level_from_string(level))

    return logger
