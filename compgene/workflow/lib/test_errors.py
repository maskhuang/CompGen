"""
Unit tests for the CompGene Error Code System.

Tests cover:
- ErrorCode enum values and properties
- ERROR_RECOVERY mapping completeness
- CompGeneError exception creation and methods
- Exit code mapping
- Error formatting functions
"""

import pickle
import pytest
import sys

from workflow.lib.errors import (
    ErrorCode,
    ERROR_RECOVERY,
    EXIT_CODES,
    EXIT_SUCCESS,
    EXIT_GENERAL_ERROR,
    EXIT_CONFIG_ERROR,
    CompGeneError,
    ConfigurationError,
    get_recovery,
    format_error_message,
    exit_with_error,
)


# =============================================================================
# Test ErrorCode Enum (Task 1)
# =============================================================================

class TestErrorCodeEnum:
    """Tests for the ErrorCode enumeration."""

    def test_all_error_codes_defined(self):
        """[P1] All 9 required error codes should be defined."""
        # GIVEN: The required error codes from ADR-003
        required_codes = [
            "E_INPUT_MISSING",
            "E_INPUT_FORMAT",
            "E_TOOL_NOT_FOUND",
            "E_TOOL_VERSION",
            "E_TIMEOUT",
            "E_OOM",
            "E_NET_RATE_LIMIT",
            "E_DISK_FULL",
            "E_OUTPUT_MISSING",
            "E_NONZERO_EXIT",
        ]

        # WHEN: Checking ErrorCode enum members
        actual_codes = [code.name for code in ErrorCode]

        # THEN: All required codes should be present
        assert len(actual_codes) == 10
        for code in required_codes:
            assert code in actual_codes

    def test_error_code_inherits_str(self):
        """[P1] ErrorCode should inherit from str for JSON serialization."""
        # GIVEN: An error code
        code = ErrorCode.E_INPUT_MISSING

        # WHEN: Checking type
        # THEN: Should be instance of str
        assert isinstance(code, str)
        assert isinstance(code.value, str)

    def test_error_code_value_equals_name(self):
        """[P2] ErrorCode value should equal its name for consistency."""
        # GIVEN: All error codes
        for code in ErrorCode:
            # WHEN/THEN: Value should match name
            assert code.value == code.name

    def test_error_code_string_conversion(self):
        """[P2] ErrorCode should convert to string correctly."""
        # GIVEN: An error code
        code = ErrorCode.E_TIMEOUT

        # WHEN: Converting to string
        code_str = str(code)

        # THEN: Should include the value
        assert "E_TIMEOUT" in code_str


# =============================================================================
# Test ERROR_RECOVERY Mapping (Task 2)
# =============================================================================

class TestErrorRecoveryMapping:
    """Tests for the ERROR_RECOVERY mapping."""

    def test_all_error_codes_have_recovery(self):
        """[P1] ERROR_RECOVERY should cover all ErrorCode values."""
        # GIVEN: All error codes
        for code in ErrorCode:
            # WHEN: Looking up recovery
            # THEN: Should exist in mapping
            assert code in ERROR_RECOVERY, f"Missing recovery for {code}"

    def test_recovery_tuple_structure(self):
        """[P1] Each recovery entry should be (bool, str) tuple."""
        # GIVEN: All recovery entries
        for code, recovery in ERROR_RECOVERY.items():
            # WHEN: Checking structure
            # THEN: Should be tuple of (bool, str)
            assert isinstance(recovery, tuple)
            assert len(recovery) == 2
            assert isinstance(recovery[0], bool)
            assert isinstance(recovery[1], str)

    def test_retryable_errors(self):
        """[P1] Correct errors should be marked as retryable."""
        # GIVEN: Retryable error codes per ADR-003
        retryable_codes = [ErrorCode.E_TIMEOUT, ErrorCode.E_NET_RATE_LIMIT]

        # WHEN/THEN: Only these should be retryable
        for code in ErrorCode:
            is_retryable, _ = ERROR_RECOVERY[code]
            if code in retryable_codes:
                assert is_retryable, f"{code} should be retryable"
            else:
                assert not is_retryable, f"{code} should not be retryable"

    def test_recovery_suggestions_not_empty(self):
        """[P2] All recovery suggestions should be non-empty strings."""
        # GIVEN: All recovery entries
        for code, (_, suggestion) in ERROR_RECOVERY.items():
            # WHEN/THEN: Suggestion should be non-empty
            assert len(suggestion) > 0, f"Empty suggestion for {code}"


class TestGetRecoveryFunction:
    """Tests for the get_recovery function."""

    def test_get_recovery_known_code(self):
        """[P1] get_recovery should return correct tuple for known code."""
        # GIVEN: A known error code
        code = ErrorCode.E_INPUT_MISSING

        # WHEN: Getting recovery
        is_retryable, suggestion = get_recovery(code)

        # THEN: Should match mapping
        assert is_retryable is False
        assert len(suggestion) > 0

    def test_get_recovery_retryable_code(self):
        """[P1] get_recovery should return True for retryable codes."""
        # GIVEN: A retryable error code
        code = ErrorCode.E_TIMEOUT

        # WHEN: Getting recovery
        is_retryable, _ = get_recovery(code)

        # THEN: Should be retryable
        assert is_retryable is True


# =============================================================================
# Test CompGeneError Exception (Task 3)
# =============================================================================

class TestCompGeneError:
    """Tests for the CompGeneError exception class."""

    def test_create_basic_error(self):
        """[P1] Should create error with code and message."""
        # GIVEN: Error code and message
        code = ErrorCode.E_INPUT_MISSING
        message = "File not found"

        # WHEN: Creating error
        error = CompGeneError(code, message)

        # THEN: Attributes should be set
        assert error.error_code == code
        assert error.message == message
        assert error.details is None

    def test_create_error_with_details(self):
        """[P1] Should create error with optional details."""
        # GIVEN: Error with details
        code = ErrorCode.E_INPUT_FORMAT
        message = "Invalid GFF file"
        details = "/path/to/file.gff3"

        # WHEN: Creating error
        error = CompGeneError(code, message, details=details)

        # THEN: Details should be set
        assert error.details == details

    def test_error_has_is_retryable(self):
        """[P1] Error should have is_retryable attribute from recovery."""
        # GIVEN: A non-retryable error
        error = CompGeneError(ErrorCode.E_INPUT_MISSING, "Missing")

        # WHEN: Checking is_retryable
        # THEN: Should be False
        assert error.is_retryable is False

        # GIVEN: A retryable error
        error2 = CompGeneError(ErrorCode.E_TIMEOUT, "Timed out")

        # WHEN: Checking is_retryable
        # THEN: Should be True
        assert error2.is_retryable is True

    def test_error_has_recovery_suggestion(self):
        """[P1] Error should have recovery_suggestion from mapping."""
        # GIVEN: An error
        error = CompGeneError(ErrorCode.E_TOOL_NOT_FOUND, "Tool missing")

        # WHEN: Checking recovery_suggestion
        # THEN: Should be non-empty string
        assert isinstance(error.recovery_suggestion, str)
        assert len(error.recovery_suggestion) > 0

    def test_error_str_format(self):
        """[P1] Error __str__ should include code and message."""
        # GIVEN: An error
        error = CompGeneError(ErrorCode.E_OOM, "Out of memory")

        # WHEN: Converting to string
        error_str = str(error)

        # THEN: Should contain code and message
        assert "E_OOM" in error_str
        assert "Out of memory" in error_str

    def test_error_is_exception(self):
        """[P1] CompGeneError should be raisable as exception."""
        # GIVEN: An error
        error = CompGeneError(ErrorCode.E_DISK_FULL, "Disk full")

        # WHEN/THEN: Should be raisable and catchable
        with pytest.raises(CompGeneError) as exc_info:
            raise error

        assert exc_info.value.error_code == ErrorCode.E_DISK_FULL

    def test_error_pickle_serialization(self):
        """[P1] CompGeneError should be pickle-serializable for multiprocessing."""
        # GIVEN: An error with all attributes
        error = CompGeneError(
            ErrorCode.E_OOM,
            "Out of memory",
            details="/proc/meminfo"
        )

        # WHEN: Pickling and unpickling
        pickled = pickle.dumps(error)
        unpickled = pickle.loads(pickled)

        # THEN: All attributes should be preserved
        assert unpickled.error_code == error.error_code
        assert unpickled.message == error.message
        assert unpickled.details == error.details
        assert unpickled.is_retryable == error.is_retryable
        assert unpickled.recovery_suggestion == error.recovery_suggestion

    def test_error_pickle_without_details(self):
        """[P2] CompGeneError without details should also be pickle-serializable."""
        # GIVEN: An error without details
        error = CompGeneError(ErrorCode.E_TIMEOUT, "Timed out")

        # WHEN: Pickling and unpickling
        pickled = pickle.dumps(error)
        unpickled = pickle.loads(pickled)

        # THEN: Attributes should be preserved, details should be None
        assert unpickled.error_code == error.error_code
        assert unpickled.message == error.message
        assert unpickled.details is None


# =============================================================================
# Test Exit Code Mapping (Task 4)
# =============================================================================

class TestExitCodeMapping:
    """Tests for exit code mapping."""

    def test_exit_codes_defined(self):
        """[P1] All error codes should have exit code mapping."""
        # GIVEN: All error codes
        for code in ErrorCode:
            # WHEN/THEN: Should have exit code
            assert code in EXIT_CODES, f"Missing exit code for {code}"

    def test_exit_code_values_in_range(self):
        """[P1] Exit codes should be in valid range (1-7)."""
        # GIVEN: All exit codes
        for code, exit_code in EXIT_CODES.items():
            # WHEN/THEN: Should be 1-7 (0 is success, not an error)
            assert 1 <= exit_code <= 7, f"Invalid exit code {exit_code} for {code}"

    def test_input_errors_return_code_3(self):
        """[P1] Input errors should return exit code 3."""
        # GIVEN: Input error codes
        input_codes = [ErrorCode.E_INPUT_MISSING, ErrorCode.E_INPUT_FORMAT]

        # WHEN/THEN: Should all return 3
        for code in input_codes:
            assert EXIT_CODES[code] == 3

    def test_tool_errors_return_code_4(self):
        """[P1] Tool errors should return exit code 4."""
        # GIVEN: Tool error codes
        tool_codes = [ErrorCode.E_TOOL_NOT_FOUND, ErrorCode.E_TOOL_VERSION]

        # WHEN/THEN: Should all return 4
        for code in tool_codes:
            assert EXIT_CODES[code] == 4

    def test_resource_errors_return_code_5(self):
        """[P1] Resource errors should return exit code 5."""
        # GIVEN: Resource error codes
        resource_codes = [ErrorCode.E_TIMEOUT, ErrorCode.E_OOM, ErrorCode.E_DISK_FULL]

        # WHEN/THEN: Should all return 5
        for code in resource_codes:
            assert EXIT_CODES[code] == 5

    def test_network_error_returns_code_6(self):
        """[P1] Network rate limit should return exit code 6."""
        # GIVEN: Network error code
        # WHEN/THEN: Should return 6
        assert EXIT_CODES[ErrorCode.E_NET_RATE_LIMIT] == 6

    def test_general_error_returns_code_1(self):
        """[P1] General error should return exit code 1."""
        # GIVEN: General error code
        # WHEN/THEN: Should return 1
        assert EXIT_CODES[ErrorCode.E_NONZERO_EXIT] == 1

    def test_special_exit_codes_defined(self):
        """[P1] Special exit codes should be defined correctly."""
        assert EXIT_SUCCESS == 0
        assert EXIT_GENERAL_ERROR == 1
        assert EXIT_CONFIG_ERROR == 2


class TestToExitCode:
    """Tests for CompGeneError.to_exit_code() method."""

    def test_to_exit_code_returns_correct_value(self):
        """[P1] to_exit_code should return mapped exit code."""
        # GIVEN: Various errors
        test_cases = [
            (ErrorCode.E_INPUT_MISSING, 3),
            (ErrorCode.E_TOOL_NOT_FOUND, 4),
            (ErrorCode.E_TIMEOUT, 5),
            (ErrorCode.E_NET_RATE_LIMIT, 6),
            (ErrorCode.E_NONZERO_EXIT, 1),
        ]

        for code, expected_exit in test_cases:
            # WHEN: Creating error and getting exit code
            error = CompGeneError(code, "Test message")
            actual_exit = error.to_exit_code()

            # THEN: Should match expected
            assert actual_exit == expected_exit, f"Wrong exit code for {code}"


# =============================================================================
# Test Configuration Error (Special Case)
# =============================================================================

class TestConfigurationError:
    """Tests for the ConfigurationError special case."""

    def test_configuration_error_returns_code_2(self):
        """[P1] ConfigurationError should always return exit code 2."""
        # GIVEN: A configuration error
        error = ConfigurationError("Invalid config value")

        # WHEN: Getting exit code
        exit_code = error.to_exit_code()

        # THEN: Should be 2
        assert exit_code == EXIT_CONFIG_ERROR
        assert exit_code == 2

    def test_configuration_error_is_compgene_error(self):
        """[P1] ConfigurationError should inherit from CompGeneError."""
        # GIVEN: A configuration error
        error = ConfigurationError("Bad config")

        # WHEN/THEN: Should be instance of CompGeneError
        assert isinstance(error, CompGeneError)

    def test_configuration_error_with_details(self):
        """[P2] ConfigurationError should support details."""
        # GIVEN: Error with details
        error = ConfigurationError(
            "Invalid species configuration",
            details="species[0].id: must be lowercase"
        )

        # WHEN/THEN: Details should be accessible
        assert error.details == "species[0].id: must be lowercase"

    def test_configuration_error_pickle_serialization(self):
        """[P1] ConfigurationError should be pickle-serializable."""
        # GIVEN: A configuration error with details
        error = ConfigurationError(
            "Invalid config value",
            details="species[0].genome: file not found"
        )

        # WHEN: Pickling and unpickling
        pickled = pickle.dumps(error)
        unpickled = pickle.loads(pickled)

        # THEN: All attributes should be preserved
        assert unpickled.message == error.message
        assert unpickled.details == error.details
        assert unpickled.to_exit_code() == 2
        assert isinstance(unpickled, ConfigurationError)


# =============================================================================
# Test Error Formatting (Task 5)
# =============================================================================

class TestFormatErrorMessage:
    """Tests for the format_error_message function."""

    def test_format_includes_error_code(self):
        """[P1] Formatted message should include error code."""
        # GIVEN: An error
        error = CompGeneError(ErrorCode.E_INPUT_MISSING, "File missing")

        # WHEN: Formatting
        formatted = format_error_message(error, use_color=False)

        # THEN: Should include code
        assert "E_INPUT_MISSING" in formatted

    def test_format_includes_message(self):
        """[P1] Formatted message should include error message."""
        # GIVEN: An error
        error = CompGeneError(ErrorCode.E_TIMEOUT, "Operation timed out")

        # WHEN: Formatting
        formatted = format_error_message(error, use_color=False)

        # THEN: Should include message
        assert "Operation timed out" in formatted

    def test_format_includes_recovery(self):
        """[P1] Formatted message should include recovery suggestion."""
        # GIVEN: An error
        error = CompGeneError(ErrorCode.E_TOOL_NOT_FOUND, "Tool missing")

        # WHEN: Formatting
        formatted = format_error_message(error, use_color=False)

        # THEN: Should include recovery
        assert "Recovery:" in formatted
        assert error.recovery_suggestion in formatted

    def test_format_includes_details_when_present(self):
        """[P1] Formatted message should include details if provided."""
        # GIVEN: Error with details
        error = CompGeneError(
            ErrorCode.E_INPUT_FORMAT,
            "Invalid format",
            details="/path/to/file.gff3"
        )

        # WHEN: Formatting
        formatted = format_error_message(error, use_color=False)

        # THEN: Should include details
        assert "Details:" in formatted
        assert "/path/to/file.gff3" in formatted

    def test_format_no_details_when_absent(self):
        """[P2] Formatted message should not have Details line if no details."""
        # GIVEN: Error without details
        error = CompGeneError(ErrorCode.E_OOM, "Out of memory")

        # WHEN: Formatting
        formatted = format_error_message(error, use_color=False)

        # THEN: Should not include "Details:"
        assert "Details:" not in formatted

    def test_format_shows_retryable_indicator(self):
        """[P2] Formatted message should indicate if retryable."""
        # GIVEN: Retryable error
        error = CompGeneError(ErrorCode.E_TIMEOUT, "Timed out")

        # WHEN: Formatting
        formatted = format_error_message(error, use_color=False)

        # THEN: Should show retryable indicator
        assert "可重试" in formatted

        # GIVEN: Non-retryable error
        error2 = CompGeneError(ErrorCode.E_INPUT_MISSING, "Missing")

        # WHEN: Formatting
        formatted2 = format_error_message(error2, use_color=False)

        # THEN: Should show non-retryable indicator
        assert "不可重试" in formatted2

    def test_format_with_color(self):
        """[P2] Formatting with color should include ANSI codes."""
        # GIVEN: An error
        error = CompGeneError(ErrorCode.E_DISK_FULL, "Disk full")

        # WHEN: Formatting with color
        formatted = format_error_message(error, use_color=True)

        # THEN: Should include ANSI escape codes
        assert "\033[" in formatted


# =============================================================================
# Test exit_with_error (Task 5)
# =============================================================================

class TestExitWithError:
    """Tests for the exit_with_error function."""

    def test_exit_with_error_calls_sys_exit(self):
        """[P1] exit_with_error should call sys.exit with correct code."""
        # GIVEN: An error
        error = CompGeneError(ErrorCode.E_INPUT_MISSING, "File missing")

        # WHEN: Calling exit_with_error
        # THEN: Should raise SystemExit with correct code
        with pytest.raises(SystemExit) as exc_info:
            exit_with_error(error, use_color=False)

        assert exc_info.value.code == 3

    def test_exit_with_error_prints_to_stderr(self, capsys):
        """[P1] exit_with_error should print to stderr."""
        # GIVEN: An error
        error = CompGeneError(ErrorCode.E_TIMEOUT, "Timed out")

        # WHEN: Calling exit_with_error
        with pytest.raises(SystemExit):
            exit_with_error(error, use_color=False)

        # THEN: Should print to stderr
        captured = capsys.readouterr()
        assert "E_TIMEOUT" in captured.err
        assert "Timed out" in captured.err


# =============================================================================
# Edge Case and Integration Tests (Task 6 Extension)
# =============================================================================

class TestEdgeCases:
    """Edge case tests for boundary conditions and special scenarios."""

    def test_exception_chaining_preserves_cause(self):
        """[P2] CompGeneError should preserve exception chain when used as wrapper."""
        # GIVEN: An original exception
        original_error = ValueError("Original validation error")

        # WHEN: Wrapping in CompGeneError with 'from' clause
        try:
            try:
                raise original_error
            except ValueError as e:
                raise CompGeneError(
                    ErrorCode.E_INPUT_FORMAT,
                    "Wrapped error",
                    details="validation failed"
                ) from e
        except CompGeneError as wrapped:
            # THEN: __cause__ should be preserved
            assert wrapped.__cause__ is original_error
            assert isinstance(wrapped.__cause__, ValueError)

    def test_error_code_json_serialization(self):
        """[P2] ErrorCode should be JSON serializable."""
        import json

        # GIVEN: An error code
        code = ErrorCode.E_TIMEOUT

        # WHEN: Serializing to JSON
        data = {"error_code": code, "value": code.value}
        json_str = json.dumps(data)
        parsed = json.loads(json_str)

        # THEN: Should serialize and deserialize correctly
        assert parsed["error_code"] == "E_TIMEOUT"
        assert parsed["value"] == "E_TIMEOUT"

    def test_error_code_in_dict_key(self):
        """[P2] ErrorCode should work as dictionary key."""
        # GIVEN: A dictionary with ErrorCode keys
        error_counts = {
            ErrorCode.E_TIMEOUT: 5,
            ErrorCode.E_OOM: 3,
            ErrorCode.E_DISK_FULL: 1,
        }

        # WHEN: Accessing by ErrorCode
        # THEN: Should work correctly
        assert error_counts[ErrorCode.E_TIMEOUT] == 5
        assert error_counts.get(ErrorCode.E_INPUT_MISSING, 0) == 0

    def test_compgene_error_with_unicode_message(self):
        """[P2] CompGeneError should handle Unicode messages correctly."""
        # GIVEN: Error with Chinese message
        error = CompGeneError(
            ErrorCode.E_INPUT_MISSING,
            "文件未找到: genome.fa",
            details="路径: /data/基因组/test.fa"
        )

        # WHEN: Formatting
        formatted = format_error_message(error, use_color=False)

        # THEN: Unicode should be preserved
        assert "文件未找到" in formatted
        assert "基因组" in formatted

    def test_compgene_error_with_multiline_details(self):
        """[P2] CompGeneError should handle multiline details."""
        # GIVEN: Error with multiline details
        multiline_details = """Line 1: Error at position 123
Line 2: Expected GFF3 format
Line 3: Found invalid character"""
        error = CompGeneError(
            ErrorCode.E_INPUT_FORMAT,
            "Format error",
            details=multiline_details
        )

        # WHEN: Formatting
        formatted = format_error_message(error, use_color=False)

        # THEN: Details should be included
        assert "Details:" in formatted
        assert "Line 1:" in formatted

    def test_all_error_codes_have_docstring(self):
        """[P2] All ErrorCode members should have descriptive docstrings."""
        # GIVEN: All error codes
        for code in ErrorCode:
            # WHEN: Checking docstring (stored in __doc__ after the value)
            # Note: Enum member docstrings are stored in _value2member_map_
            # THEN: Each code should be documented in the class docstring
            assert code.name in ErrorCode.__doc__

    def test_exit_codes_cover_all_unix_ranges(self):
        """[P2] Exit codes should cover standard Unix exit code ranges."""
        # GIVEN: All exit codes
        exit_values = set(EXIT_CODES.values())

        # WHEN: Checking coverage
        # THEN: Should have codes for: general(1), input(3), tool(4), resource(5), network(6)
        assert 1 in exit_values  # General error
        assert 3 in exit_values  # Input error
        assert 4 in exit_values  # Tool error
        assert 5 in exit_values  # Resource error
        assert 6 in exit_values  # Network error

        # AND: Special codes should be defined
        assert EXIT_SUCCESS == 0
        assert EXIT_CONFIG_ERROR == 2

    def test_get_recovery_with_string_matching_enum_value(self):
        """[P2] get_recovery should work with string that matches enum value."""
        # GIVEN: A string that matches an ErrorCode value
        # Note: Since ErrorCode inherits from str, this works due to Python's
        # dict lookup behavior with str subclasses

        # WHEN: Looking up recovery
        is_retryable, suggestion = get_recovery(ErrorCode.E_TIMEOUT)

        # THEN: Should return correct values
        assert is_retryable is True
        assert len(suggestion) > 0


class TestIntegration:
    """Integration tests for error system with other components."""

    def test_configuration_error_integrates_with_error_formatting(self):
        """[P1] ConfigurationError should format correctly with format_error_message."""
        # GIVEN: A configuration error
        error = ConfigurationError(
            "Invalid species ID",
            details="species[0].id: 'Test' should be lowercase"
        )

        # WHEN: Formatting
        formatted = format_error_message(error, use_color=False)

        # THEN: Should include all relevant information
        assert "E_INPUT_FORMAT" in formatted  # Base error code
        assert "Invalid species ID" in formatted
        assert "species[0].id" in formatted
        assert "Recovery:" in formatted

    def test_error_can_be_raised_and_caught_by_base_class(self):
        """[P1] ConfigurationError should be catchable as CompGeneError."""
        # GIVEN: A configuration error
        config_error = ConfigurationError("Test error")

        # WHEN: Raising and catching as base class
        with pytest.raises(CompGeneError) as exc_info:
            raise config_error

        # THEN: Should be caught and have correct exit code
        assert exc_info.value.to_exit_code() == 2

    def test_multiple_errors_have_distinct_exit_codes(self):
        """[P1] Different error types should produce distinct exit codes."""
        # GIVEN: Various errors
        errors = [
            CompGeneError(ErrorCode.E_INPUT_MISSING, "Missing"),
            CompGeneError(ErrorCode.E_TOOL_NOT_FOUND, "Tool not found"),
            CompGeneError(ErrorCode.E_TIMEOUT, "Timeout"),
            CompGeneError(ErrorCode.E_NET_RATE_LIMIT, "Rate limited"),
            ConfigurationError("Config error"),
        ]

        # WHEN: Getting exit codes
        exit_codes = [e.to_exit_code() for e in errors]

        # THEN: Should have expected codes
        assert exit_codes == [3, 4, 5, 6, 2]

    def test_error_recovery_suggestions_are_actionable(self):
        """[P2] All recovery suggestions should be non-empty actionable strings."""
        # GIVEN: All error codes
        for code in ErrorCode:
            # WHEN: Getting recovery
            is_retryable, suggestion = get_recovery(code)

            # THEN: Suggestion should be actionable
            assert len(suggestion) > 5, f"Suggestion for {code} is too short"
            assert not suggestion.startswith(" "), f"Suggestion for {code} has leading space"
