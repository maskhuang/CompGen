"""
Unit tests for the NCBI download module.

Tests cover:
- Accession parsing (GCF/GCA formats)
- Rate limiter behavior
- Retry logic
- Cache mechanism
- Error handling

Source: Story 2.3 - NCBI Annotation Download
"""

import json
import os
import tempfile
import threading
import time
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch
import urllib.error

import pytest

from workflow.lib.errors import CompGeneError, ErrorCode
from workflow.lib.ncbi import (
    # Constants
    DEFAULT_RATE_LIMIT,
    API_KEY_RATE_LIMIT,
    CACHE_MANIFEST,
    # Classes
    AccessionType,
    NCBIAccession,
    RateLimiter,
    DownloadResult,
    NCBIDownloadResult,
    NCBIClient,
    # Functions
    parse_ncbi_accession,
    is_ncbi_accession,
    get_ncbi_api_key,
    build_url_with_api_key,
    with_retry,
    is_cached,
    verify_cache,
    write_download_manifest,
)


# =============================================================================
# Accession Parsing Tests
# =============================================================================


class TestParseNCBIAccession:
    """Tests for parse_ncbi_accession function."""

    def test_parse_gcf_with_prefix_and_version(self):
        """Parse GCF accession with NCBI: prefix and version."""
        acc = parse_ncbi_accession("NCBI:GCF_000165445.1")

        assert acc.prefix == "GCF"
        assert acc.numeric == "000165445"
        assert acc.version == "1"
        assert acc.accession_type == AccessionType.REFSEQ
        assert acc.full_accession == "GCF_000165445.1"

    def test_parse_gcf_without_prefix(self):
        """Parse GCF accession without NCBI: prefix."""
        acc = parse_ncbi_accession("GCF_000165445.1")

        assert acc.prefix == "GCF"
        assert acc.numeric == "000165445"
        assert acc.version == "1"
        assert acc.full_accession == "GCF_000165445.1"

    def test_parse_gca_genbank(self):
        """Parse GCA (GenBank) accession."""
        acc = parse_ncbi_accession("GCA_000001405.15")

        assert acc.prefix == "GCA"
        assert acc.numeric == "000001405"
        assert acc.version == "15"
        assert acc.accession_type == AccessionType.GENBANK

    def test_parse_without_version(self):
        """Parse accession without version number."""
        acc = parse_ncbi_accession("GCF_000165445")

        assert acc.prefix == "GCF"
        assert acc.numeric == "000165445"
        assert acc.version is None
        assert acc.full_accession == "GCF_000165445"

    def test_parse_case_insensitive(self):
        """Verify case-insensitive parsing."""
        acc_lower = parse_ncbi_accession("gcf_000165445.1")
        acc_upper = parse_ncbi_accession("GCF_000165445.1")

        assert acc_lower.prefix == "GCF"
        assert acc_lower.full_accession == acc_upper.full_accession

    def test_ftp_path_parts(self):
        """Test FTP path part generation."""
        acc = parse_ncbi_accession("GCF_000165445.1")

        assert acc.ftp_path_parts == ("000", "165", "445")

    def test_ftp_base_path(self):
        """Test FTP base path generation."""
        acc = parse_ncbi_accession("GCF_000165445.1")

        assert acc.ftp_base_path == "GCF/000/165/445"

    def test_base_accession(self):
        """Test base accession without version."""
        acc = parse_ncbi_accession("GCF_000165445.1")

        assert acc.base_accession == "GCF_000165445"

    def test_invalid_format_raises_error(self):
        """Invalid format raises CompGeneError."""
        with pytest.raises(CompGeneError) as exc_info:
            parse_ncbi_accession("invalid_accession")

        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT

    def test_empty_string_raises_error(self):
        """Empty string raises CompGeneError."""
        with pytest.raises(CompGeneError) as exc_info:
            parse_ncbi_accession("")

        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT

    def test_wrong_prefix_raises_error(self):
        """Wrong prefix (not GCF/GCA) raises error."""
        with pytest.raises(CompGeneError) as exc_info:
            parse_ncbi_accession("GCX_000165445.1")

        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT

    def test_invalid_numeric_length(self):
        """Numeric part with wrong length raises error."""
        with pytest.raises(CompGeneError) as exc_info:
            parse_ncbi_accession("GCF_12345678.1")  # 8 digits instead of 9

        assert exc_info.value.error_code == ErrorCode.E_INPUT_FORMAT


class TestIsNCBIAccession:
    """Tests for is_ncbi_accession function."""

    def test_valid_gcf_returns_true(self):
        """Valid GCF accession returns True."""
        assert is_ncbi_accession("GCF_000165445.1") is True

    def test_valid_gca_returns_true(self):
        """Valid GCA accession returns True."""
        assert is_ncbi_accession("GCA_000001405.15") is True

    def test_with_ncbi_prefix_returns_true(self):
        """Accession with NCBI: prefix returns True."""
        assert is_ncbi_accession("NCBI:GCF_000165445.1") is True

    def test_invalid_returns_false(self):
        """Invalid string returns False."""
        assert is_ncbi_accession("not_an_accession") is False

    def test_empty_returns_false(self):
        """Empty string returns False."""
        assert is_ncbi_accession("") is False

    def test_none_returns_false(self):
        """None returns False."""
        assert is_ncbi_accession(None) is False

    def test_local_path_returns_false(self):
        """Local file path returns False."""
        assert is_ncbi_accession("/path/to/file.gff") is False


# =============================================================================
# Rate Limiter Tests
# =============================================================================


class TestRateLimiter:
    """Tests for RateLimiter class."""

    def test_initial_tokens(self):
        """Rate limiter starts with full tokens."""
        limiter = RateLimiter(3.0)

        # Should not wait for first few requests
        wait1 = limiter.acquire()
        assert wait1 == 0.0

        wait2 = limiter.acquire()
        assert wait2 == 0.0

        wait3 = limiter.acquire()
        assert wait3 == 0.0

    def test_rate_limiting_triggers(self):
        """Rate limiter enforces waiting after tokens exhausted."""
        limiter = RateLimiter(10.0)  # 10 req/s for faster test

        # Exhaust initial tokens
        for _ in range(10):
            limiter.acquire()

        # Next request should wait
        start = time.monotonic()
        limiter.acquire()
        elapsed = time.monotonic() - start

        # Should have waited ~0.1 seconds (1/10 req/s)
        assert elapsed >= 0.05  # Allow some tolerance

    def test_tokens_refill(self):
        """Tokens refill over time."""
        limiter = RateLimiter(10.0)

        # Exhaust tokens
        for _ in range(10):
            limiter.acquire()

        # Wait for refill
        time.sleep(0.3)  # Should refill ~3 tokens

        # Should be able to acquire without much wait
        wait = limiter.acquire()
        assert wait < 0.1

    def test_set_rate(self):
        """Rate can be updated."""
        limiter = RateLimiter(3.0)

        limiter.set_rate(10.0)

        assert limiter.rate == 10.0

    def test_thread_safety(self):
        """Rate limiter is thread-safe."""
        limiter = RateLimiter(100.0)  # High rate for fast test
        results = []

        def acquire_token():
            limiter.acquire()
            results.append(True)

        threads = [
            threading.Thread(target=acquire_token)
            for _ in range(50)
        ]

        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert len(results) == 50


# =============================================================================
# Retry Logic Tests
# =============================================================================


class TestWithRetry:
    """Tests for with_retry function."""

    def test_success_on_first_try(self):
        """Successful function returns immediately."""
        mock_func = Mock(return_value="success")

        result = with_retry(mock_func, max_retries=3)

        assert result == "success"
        assert mock_func.call_count == 1

    def test_retries_on_failure(self):
        """Retries on retryable exceptions."""
        mock_func = Mock(
            side_effect=[
                urllib.error.URLError("error 1"),
                urllib.error.URLError("error 2"),
                "success"
            ]
        )

        result = with_retry(
            mock_func,
            max_retries=3,
            base_delay=0.01,  # Fast for tests
            max_delay=0.1
        )

        assert result == "success"
        assert mock_func.call_count == 3

    def test_raises_after_max_retries(self):
        """Raises exception after max retries exhausted."""
        mock_func = Mock(
            side_effect=urllib.error.URLError("persistent error")
        )

        with pytest.raises(urllib.error.URLError):
            with_retry(
                mock_func,
                max_retries=2,
                base_delay=0.01,
                max_delay=0.1
            )

        assert mock_func.call_count == 3  # Initial + 2 retries

    def test_exponential_backoff(self):
        """Delay increases exponentially."""
        delays = []

        def capturing_sleep(duration):
            delays.append(duration)

        mock_func = Mock(
            side_effect=urllib.error.URLError("error")
        )

        with patch("time.sleep", capturing_sleep):
            with pytest.raises(urllib.error.URLError):
                with_retry(
                    mock_func,
                    max_retries=3,
                    base_delay=1.0,
                    max_delay=100.0
                )

        # Delays should be: 1.0, 2.0, 4.0
        assert delays == [1.0, 2.0, 4.0]

    def test_max_delay_cap(self):
        """Delay is capped at max_delay."""
        delays = []

        def capturing_sleep(duration):
            delays.append(duration)

        mock_func = Mock(
            side_effect=urllib.error.URLError("error")
        )

        with patch("time.sleep", capturing_sleep):
            with pytest.raises(urllib.error.URLError):
                with_retry(
                    mock_func,
                    max_retries=5,
                    base_delay=1.0,
                    max_delay=5.0
                )

        # Delays should be: 1.0, 2.0, 4.0, 5.0, 5.0 (capped)
        assert delays == [1.0, 2.0, 4.0, 5.0, 5.0]

    def test_non_retryable_exception_not_retried(self):
        """Non-retryable exceptions are not retried."""
        mock_func = Mock(side_effect=ValueError("not retryable"))

        with pytest.raises(ValueError):
            with_retry(
                mock_func,
                max_retries=3,
                retryable_exceptions=(urllib.error.URLError,)
            )

        assert mock_func.call_count == 1


# =============================================================================
# API Key Tests
# =============================================================================


class TestAPIKey:
    """Tests for API key handling."""

    def test_get_api_key_from_env(self):
        """API key is read from environment."""
        with patch.dict(os.environ, {"NCBI_API_KEY": "test_key_123"}):
            key = get_ncbi_api_key()
            assert key == "test_key_123"

    def test_get_api_key_when_not_set(self):
        """Returns None when env var not set."""
        with patch.dict(os.environ, {}, clear=True):
            # Remove NCBI_API_KEY if present
            os.environ.pop("NCBI_API_KEY", None)
            key = get_ncbi_api_key()
            assert key is None

    def test_build_url_with_api_key(self):
        """API key is appended to URL."""
        url = build_url_with_api_key(
            "https://api.ncbi.nlm.nih.gov/test",
            api_key="my_key"
        )

        assert "api_key=my_key" in url

    def test_build_url_without_api_key(self):
        """URL unchanged when no API key."""
        with patch.dict(os.environ, {}, clear=True):
            os.environ.pop("NCBI_API_KEY", None)
            url = build_url_with_api_key("https://api.ncbi.nlm.nih.gov/test")

            assert "api_key" not in url

    def test_build_url_existing_query_params(self):
        """API key uses & for URL with existing params."""
        url = build_url_with_api_key(
            "https://api.ncbi.nlm.nih.gov/test?param=value",
            api_key="my_key"
        )

        assert "&api_key=my_key" in url


# =============================================================================
# Cache Tests
# =============================================================================


class TestCacheManagement:
    """Tests for cache management functions."""

    def test_is_cached_returns_false_when_no_manifest(self):
        """is_cached returns False when manifest doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            assert is_cached("GCF_000165445.1", Path(tmpdir)) is False

    def test_is_cached_returns_false_when_files_missing(self):
        """is_cached returns False when cached files are missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            accession_dir = cache_dir / "GCF_000165445.1"
            accession_dir.mkdir()

            # Create manifest with missing file reference
            manifest = {
                "files": [
                    {"path": str(accession_dir / "missing.gz")}
                ]
            }

            with open(accession_dir / CACHE_MANIFEST, "w") as f:
                json.dump(manifest, f)

            assert is_cached("GCF_000165445.1", cache_dir) is False

    def test_is_cached_returns_true_when_all_files_exist(self):
        """is_cached returns True when all files exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            accession_dir = cache_dir / "GCF_000165445.1"
            accession_dir.mkdir()

            # Create test file
            test_file = accession_dir / "genome.fna.gz"
            test_file.write_text("test content")

            # Create manifest
            manifest = {
                "files": [
                    {"path": str(test_file)}
                ]
            }

            with open(accession_dir / CACHE_MANIFEST, "w") as f:
                json.dump(manifest, f)

            assert is_cached("GCF_000165445.1", cache_dir) is True

    def test_verify_cache_with_valid_checksum(self):
        """verify_cache returns True for valid checksums."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            accession_dir = cache_dir / "GCF_000165445.1"
            accession_dir.mkdir()

            # Create test file
            test_file = accession_dir / "genome.fna.gz"
            content = b"test content for checksum"
            test_file.write_bytes(content)

            # Calculate correct checksum
            import hashlib
            expected_sha256 = hashlib.sha256(content).hexdigest()

            # Create manifest with correct checksum
            manifest = {
                "files": [
                    {
                        "path": str(test_file),
                        "sha256": expected_sha256
                    }
                ]
            }

            with open(accession_dir / CACHE_MANIFEST, "w") as f:
                json.dump(manifest, f)

            assert verify_cache("GCF_000165445.1", cache_dir) is True

    def test_verify_cache_with_invalid_checksum(self):
        """verify_cache returns False for invalid checksums."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            accession_dir = cache_dir / "GCF_000165445.1"
            accession_dir.mkdir()

            # Create test file
            test_file = accession_dir / "genome.fna.gz"
            test_file.write_bytes(b"test content")

            # Create manifest with wrong checksum
            manifest = {
                "files": [
                    {
                        "path": str(test_file),
                        "sha256": "wrong_checksum_value"
                    }
                ]
            }

            with open(accession_dir / CACHE_MANIFEST, "w") as f:
                json.dump(manifest, f)

            assert verify_cache("GCF_000165445.1", cache_dir) is False

    def test_write_download_manifest(self):
        """write_download_manifest creates valid manifest."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)

            files = [
                {
                    "type": "genome",
                    "path": str(output_dir / "genome.fna.gz"),
                    "sha256": "abc123"
                }
            ]

            manifest_path = write_download_manifest(
                "GCF_000165445.1",
                output_dir,
                files,
                metadata={"source": "test"}
            )

            assert manifest_path.exists()

            with open(manifest_path) as f:
                manifest = json.load(f)

            assert manifest["accession"] == "GCF_000165445.1"
            assert len(manifest["files"]) == 1
            assert manifest["metadata"]["source"] == "test"
            assert "downloaded_at" in manifest


# =============================================================================
# NCBIClient Tests
# =============================================================================


class TestNCBIClient:
    """Tests for NCBIClient class."""

    def test_client_initialization(self):
        """Client initializes with cache directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            client = NCBIClient(Path(tmpdir))

            assert client.cache_dir == Path(tmpdir)
            assert client.rate_limiter is not None

    def test_client_with_api_key(self):
        """Client with API key uses higher rate limit."""
        with tempfile.TemporaryDirectory() as tmpdir:
            client = NCBIClient(Path(tmpdir), api_key="test_key")

            assert client.api_key == "test_key"
            assert client.rate_limiter.rate == API_KEY_RATE_LIMIT

    def test_client_without_api_key(self):
        """Client without API key uses default rate limit."""
        with tempfile.TemporaryDirectory() as tmpdir:
            with patch.dict(os.environ, {}, clear=True):
                os.environ.pop("NCBI_API_KEY", None)
                client = NCBIClient(Path(tmpdir))

                assert client.rate_limiter.rate == DEFAULT_RATE_LIMIT

    def test_client_is_cached(self):
        """Client.is_cached checks cache correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            client = NCBIClient(Path(tmpdir))

            # Not cached initially
            assert client.is_cached("GCF_000165445.1") is False

    def test_client_clear_cache(self):
        """Client.clear_cache removes cached data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            client = NCBIClient(Path(tmpdir))
            accession_dir = client.cache_dir / "GCF_000165445.1"
            accession_dir.mkdir()

            # Create a file in cache
            (accession_dir / "test.txt").write_text("test")

            assert accession_dir.exists()

            result = client.clear_cache("GCF_000165445.1")

            assert result is True
            assert not accession_dir.exists()

    def test_client_clear_cache_nonexistent(self):
        """Client.clear_cache returns False for non-existent cache."""
        with tempfile.TemporaryDirectory() as tmpdir:
            client = NCBIClient(Path(tmpdir))

            result = client.clear_cache("GCF_000165445.1")

            assert result is False


# =============================================================================
# DownloadResult Tests
# =============================================================================


class TestDownloadResult:
    """Tests for DownloadResult dataclass."""

    def test_download_result_fields(self):
        """DownloadResult has expected fields."""
        result = DownloadResult(
            url="https://example.com/file.gz",
            path="/path/to/file.gz",
            size_bytes=12345,
            sha256="abc123def456"
        )

        assert result.url == "https://example.com/file.gz"
        assert result.path == "/path/to/file.gz"
        assert result.size_bytes == 12345
        assert result.sha256 == "abc123def456"


class TestNCBIDownloadResult:
    """Tests for NCBIDownloadResult dataclass."""

    def test_ncbi_download_result_defaults(self):
        """NCBIDownloadResult has sensible defaults."""
        result = NCBIDownloadResult(
            accession="GCF_000165445.1",
            output_dir=Path("/tmp")
        )

        assert result.genome_path is None
        assert result.annotation_path is None
        assert result.protein_path is None
        assert result.from_cache is False


# =============================================================================
# Integration Test (Mocked Network)
# =============================================================================


class TestDownloadNCBIDatasetMocked:
    """Integration tests with mocked network calls."""

    @patch("workflow.lib.ncbi.fetch_url_content")
    @patch("workflow.lib.ncbi.download_file")
    def test_download_uses_cache_when_available(
        self,
        mock_download,
        mock_fetch,
    ):
        """Download uses cache when data is already cached."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            accession_dir = cache_dir / "GCF_000165445.1"
            accession_dir.mkdir()

            # Create cached files
            genome_file = accession_dir / "genome.fna.gz"
            genome_file.write_text("cached genome")

            annotation_file = accession_dir / "annotation.gff.gz"
            annotation_file.write_text("cached annotation")

            # Create manifest
            manifest = {
                "accession": "GCF_000165445.1",
                "downloaded_at": "2024-01-01T00:00:00Z",
                "files": [
                    {"type": "genome", "path": str(genome_file)},
                    {"type": "annotation", "path": str(annotation_file)},
                ]
            }

            with open(accession_dir / CACHE_MANIFEST, "w") as f:
                json.dump(manifest, f)

            # Import here to ensure fresh state
            from workflow.lib.ncbi import download_ncbi_dataset

            result = download_ncbi_dataset(
                "GCF_000165445.1",
                cache_dir,
                force=False
            )

            # Should use cache, no network calls
            mock_fetch.assert_not_called()
            mock_download.assert_not_called()

            assert result.from_cache is True
            assert result.genome_path == genome_file
            assert result.annotation_path == annotation_file

    @patch("workflow.lib.ncbi.fetch_url_content")
    @patch("workflow.lib.ncbi.download_file")
    def test_force_download_ignores_cache(
        self,
        mock_download,
        mock_fetch,
    ):
        """Force download ignores existing cache."""
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_dir = Path(tmpdir)
            accession_dir = cache_dir / "GCF_000165445.1"
            accession_dir.mkdir()

            # Create cached files
            genome_file = accession_dir / "genome.fna.gz"
            genome_file.write_text("cached genome")

            # Create manifest
            manifest = {
                "accession": "GCF_000165445.1",
                "downloaded_at": "2024-01-01T00:00:00Z",
                "files": [
                    {"type": "genome", "path": str(genome_file)},
                ]
            }

            with open(accession_dir / CACHE_MANIFEST, "w") as f:
                json.dump(manifest, f)

            # Mock FTP directory listing
            mock_fetch.return_value = b'''
            <html>
            <a href="GCF_000165445.1_Test/">GCF_000165445.1_Test/</a>
            </html>
            '''

            # Mock file download
            mock_download.return_value = MagicMock(
                url="https://example.com/genome.fna.gz",
                path=str(genome_file),
                size_bytes=1000,
                sha256="abc123"
            )

            from workflow.lib.ncbi import download_ncbi_dataset

            # Force should trigger new downloads even with cache
            result = download_ncbi_dataset(
                "GCF_000165445.1",
                cache_dir,
                force=True
            )

            # Should have made network calls despite cache
            assert mock_fetch.called


# =============================================================================
# Constants Tests
# =============================================================================


class TestConstants:
    """Tests for module constants."""

    def test_default_rate_limit(self):
        """Default rate limit is 3 requests per second."""
        assert DEFAULT_RATE_LIMIT == 3.0

    def test_api_key_rate_limit(self):
        """API key rate limit is 10 requests per second."""
        assert API_KEY_RATE_LIMIT == 10.0

    def test_ncbi_ftp_base(self):
        """NCBI FTP base URL is correct."""
        from workflow.lib.ncbi import NCBI_FTP_BASE
        assert "ftp.ncbi.nlm.nih.gov" in NCBI_FTP_BASE
