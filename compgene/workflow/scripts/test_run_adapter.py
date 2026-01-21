"""
Tests for the run_adapter.py CLI script.

Tests cover:
- Argument parsing
- JSON argument validation
- Adapter loading
- Error handling
"""

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from workflow.scripts.run_adapter import (
    parse_args,
    load_adapter,
    parse_json_arg,
    paths_from_dict,
    main,
)
from workflow.lib.errors import ErrorCode, CompGeneError


# =============================================================================
# Test parse_json_arg
# =============================================================================

class TestParseJsonArg:
    """Tests for parse_json_arg function."""

    def test_valid_json_object(self):
        """Test parsing valid JSON object."""
        result = parse_json_arg('{"key": "value"}', "test")
        assert result == {"key": "value"}

    def test_empty_json_object(self):
        """Test parsing empty JSON object."""
        result = parse_json_arg('{}', "test")
        assert result == {}

    def test_nested_json_object(self):
        """Test parsing nested JSON object."""
        result = parse_json_arg('{"outer": {"inner": "value"}}', "test")
        assert result == {"outer": {"inner": "value"}}

    def test_json_array_raises_error(self):
        """Test that JSON array raises CompGeneError."""
        with pytest.raises(CompGeneError) as exc_info:
            parse_json_arg('[1, 2, 3]', "test")
        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT
        assert "JSON object" in str(exc_info.value)

    def test_invalid_json_raises_error(self):
        """Test that invalid JSON raises CompGeneError."""
        with pytest.raises(CompGeneError) as exc_info:
            parse_json_arg('not valid json', "test")
        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT
        assert "Invalid JSON" in str(exc_info.value)

    def test_json_string_raises_error(self):
        """Test that JSON string raises CompGeneError."""
        with pytest.raises(CompGeneError) as exc_info:
            parse_json_arg('"just a string"', "test")
        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT


# =============================================================================
# Test paths_from_dict
# =============================================================================

class TestPathsFromDict:
    """Tests for paths_from_dict function."""

    def test_converts_strings_to_paths(self):
        """Test conversion of string values to Path objects."""
        result = paths_from_dict({"a": "/path/to/a", "b": "relative/b"})
        assert result == {"a": Path("/path/to/a"), "b": Path("relative/b")}

    def test_empty_dict(self):
        """Test empty dict returns empty dict."""
        assert paths_from_dict({}) == {}


# =============================================================================
# Test load_adapter
# =============================================================================

class TestLoadAdapter:
    """Tests for load_adapter function."""

    def test_module_not_found(self):
        """Test that missing module raises CompGeneError."""
        with pytest.raises(CompGeneError) as exc_info:
            load_adapter("nonexistent_adapter")
        assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND
        assert "module not found" in str(exc_info.value).lower()

    def test_class_not_found(self):
        """Test that missing class raises CompGeneError."""
        # Create a mock module without the expected class
        with patch.dict('sys.modules', {'workflow.adapters.fake': MagicMock(spec=[])}):
            with pytest.raises(CompGeneError) as exc_info:
                load_adapter("fake")
            assert exc_info.value.error_code == ErrorCode.E_TOOL_NOT_FOUND
            assert "class not found" in str(exc_info.value).lower()


# =============================================================================
# Test parse_args
# =============================================================================

class TestParseArgs:
    """Tests for parse_args function."""

    def test_minimal_args(self):
        """Test parsing minimal required arguments."""
        test_args = [
            "run_adapter.py",
            "orthofinder",
            "--inputs", '{"proteins": "data/proteins/"}',
            "--outputs", '{"result": "out.tsv"}',
        ]
        with patch.object(sys, 'argv', test_args):
            args = parse_args()
            assert args.adapter_name == "orthofinder"
            assert args.inputs == '{"proteins": "data/proteins/"}'
            assert args.outputs == '{"result": "out.tsv"}'
            assert args.config == "{}"  # default
            assert args.wildcards == "{}"  # default
            assert args.threads == 1  # default
            assert args.max_retries == 3  # default
            assert args.grace_period == 10.0  # default
            assert args.log_level == "INFO"  # default
            assert args.no_color is False  # default

    def test_all_args(self):
        """Test parsing all arguments."""
        test_args = [
            "run_adapter.py",
            "busco",
            "--inputs", '{"genome": "genome.fa"}',
            "--outputs", '{"result": "busco.txt"}',
            "--config", '{"lineage": "bacteria"}',
            "--wildcards", '{"species": "ecoli"}',
            "--threads", "8",
            "--meta-dir", "custom/meta",
            "--log-dir", "custom/logs",
            "--max-retries", "5",
            "--grace-period", "20.0",
            "--log-level", "DEBUG",
            "--no-color",
        ]
        with patch.object(sys, 'argv', test_args):
            args = parse_args()
            assert args.adapter_name == "busco"
            assert args.config == '{"lineage": "bacteria"}'
            assert args.wildcards == '{"species": "ecoli"}'
            assert args.threads == 8
            assert args.meta_dir == "custom/meta"
            assert args.log_dir == "custom/logs"
            assert args.max_retries == 5
            assert args.grace_period == 20.0
            assert args.log_level == "DEBUG"
            assert args.no_color is True


# =============================================================================
# Test main function
# =============================================================================

class TestMain:
    """Tests for main function."""

    def test_invalid_inputs_json(self):
        """Test main returns error for invalid inputs JSON."""
        test_args = [
            "run_adapter.py",
            "orthofinder",
            "--inputs", "not valid json",
            "--outputs", '{"result": "out.tsv"}',
        ]
        with patch.object(sys, 'argv', test_args):
            with patch('workflow.scripts.run_adapter.exit_with_error') as mock_exit:
                main()
                mock_exit.assert_called_once()
                error = mock_exit.call_args[0][0]
                assert error.error_code == ErrorCode.E_INPUT_FORMAT

    def test_invalid_outputs_json(self):
        """Test main returns error for invalid outputs JSON."""
        test_args = [
            "run_adapter.py",
            "orthofinder",
            "--inputs", '{"proteins": "data/"}',
            "--outputs", "not valid json",
        ]
        with patch.object(sys, 'argv', test_args):
            with patch('workflow.scripts.run_adapter.exit_with_error') as mock_exit:
                main()
                mock_exit.assert_called_once()
                error = mock_exit.call_args[0][0]
                assert error.error_code == ErrorCode.E_INPUT_FORMAT

    def test_adapter_not_found(self):
        """Test main returns error for missing adapter."""
        test_args = [
            "run_adapter.py",
            "nonexistent_adapter",
            "--inputs", '{"input": "file.txt"}',
            "--outputs", '{"output": "out.txt"}',
        ]
        with patch.object(sys, 'argv', test_args):
            with patch('workflow.scripts.run_adapter.exit_with_error') as mock_exit:
                main()
                mock_exit.assert_called_once()
                error = mock_exit.call_args[0][0]
                assert error.error_code == ErrorCode.E_TOOL_NOT_FOUND

    def test_successful_execution(self):
        """Test successful adapter execution."""
        test_args = [
            "run_adapter.py",
            "mock",
            "--inputs", '{"input": "file.txt"}',
            "--outputs", '{"output": "out.txt"}',
        ]

        # Create mock adapter and runner
        mock_adapter = MagicMock()
        mock_adapter.spec.name = "mock"
        mock_result = MagicMock()
        mock_result.outputs = {"output": Path("out.txt")}
        mock_result.summary = "success"

        with patch.object(sys, 'argv', test_args):
            with patch('workflow.scripts.run_adapter.load_adapter', return_value=mock_adapter):
                with patch('workflow.scripts.run_adapter.create_logger') as mock_logger:
                    with patch('workflow.scripts.run_adapter.AdapterRunner') as MockRunner:
                        mock_runner_instance = MagicMock()
                        mock_runner_instance.run.return_value = mock_result
                        MockRunner.return_value = mock_runner_instance

                        result = main()

                        assert result == 0
                        mock_runner_instance.run.assert_called_once()
