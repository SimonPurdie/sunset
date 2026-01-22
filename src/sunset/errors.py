"""Error handling utilities for sunset rendering."""

import sys
from typing import Optional


class SunsetError(Exception):
    """Base exception for sunset-specific errors."""

    def __init__(self, message: str, suggestions: Optional[list[str]] = None):
        """Initialize with message and optional suggestions.

        Args:
            message: Error description
            suggestions: Optional list of actionable suggestions
        """
        self.message = message
        self.suggestions = suggestions or []
        super().__init__(self._format_message())

    def _format_message(self) -> str:
        """Format error message with suggestions."""
        formatted = f"{self.message}"
        if self.suggestions:
            formatted += "\n\nSuggestions:"
            for suggestion in self.suggestions:
                formatted += f"\n  - {suggestion}"
        return formatted


class BodyNotFoundError(SunsetError):
    """Raised when a body identifier is not recognized."""

    def __init__(self, body_id: str, available_bodies: list[str]):
        message = f"Unknown body: '{body_id}'"
        suggestions = [
            f"Available bodies: {', '.join(sorted(available_bodies))}",
            "Check spelling (body names are case-insensitive)",
        ]
        super().__init__(message, suggestions)


class SunsetNotFoundError(SunsetError):
    """Raised when a sunset location cannot be found for a body."""

    def __init__(self, body_id: str, utc_time: str):
        message = f"Could not find sunset location for {body_id} at {utc_time}"
        suggestions = [
            "Try a different UTC time (sunset occurs at different times for different locations)",
            "For airless bodies, ensure the sun is below the horizon by 1.5° or more",
            "Run without --body flag to let the system auto-select a body with a visible sunset",
        ]
        super().__init__(message, suggestions)


class VisibilityError(SunsetError):
    """Raised when sun is not visible from the body."""

    def __init__(self, body_name: str, justification: str):
        message = f"Sun not visible on {body_name} at this time"
        suggestions = [
            f"Reason: {justification}",
            "Try a different UTC time when the sun might be visible",
            "Run without --body flag to try other bodies instead",
        ]
        super().__init__(message, suggestions)


class TimeParseError(SunsetError):
    """Raised when UTC time cannot be parsed."""

    def __init__(self, utc_time: str):
        message = f"Invalid UTC time format: '{utc_time}'"
        suggestions = [
            "Use ISO-8601 format with 'Z' suffix for UTC (e.g., '2024-01-15T18:00:00Z')",
            "Omit the command to use the current time",
            "Example: --utc-time 2024-01-15T18:00:00Z",
        ]
        super().__init__(message, suggestions)


class AtmosphereProfileError(SunsetError):
    """Raised when atmosphere profile is not available."""

    def __init__(self, body_name: str):
        message = f"Atmosphere profile not yet implemented for {body_name}"
        suggestions = [
            "This body may not have an atmosphere, or its atmospheric modeling is not yet complete",
            "Check BREADCRUMBS.md for implementation notes on this body",
        ]
        super().__init__(message, suggestions)


def print_error(error: Exception) -> None:
    """Print error to stderr with formatted output.

    Args:
        error: Exception to print
    """
    print(f"Error: {error}", file=sys.stderr)

    if isinstance(error, SunsetError):
        if error.suggestions:
            print(file=sys.stderr)


def handle_error(error: Exception, context: Optional[str] = None) -> int:
    """Handle an error with optional context and return exit code.

    Args:
        error: Exception that occurred
        context: Optional description of what was being attempted

    Returns:
        Exit code (1 for error)
    """
    if context:
        print(f"Error while {context}:", file=sys.stderr)

    print_error(error)

    import traceback

    if not isinstance(error, SunsetError):
        traceback.print_exc()

    return 1
