"""Tests for globsim.download.RDA resilient download behaviour."""
import os
import time
import unittest
from unittest.mock import MagicMock, patch, call

import requests

from globsim.download.RDA import (
    CONNECT_TIMEOUT,
    MAX_RETRIES,
    READ_TIMEOUT,
    RETRY_BACKOFF,
    Rdams,
)


def _make_rdams():
    """Return an Rdams instance without needing a real auth file."""
    obj = object.__new__(Rdams)
    obj.token = "fake-token"
    obj.DEFAULT_AUTH_FILE = "/dev/null"
    return obj


class TestModuleConstants(unittest.TestCase):
    """Sanity-check the module-level timeout/retry constants."""

    def test_connect_timeout_positive(self):
        self.assertGreater(CONNECT_TIMEOUT, 0)

    def test_read_timeout_positive(self):
        self.assertGreater(READ_TIMEOUT, 0)

    def test_max_retries_positive(self):
        self.assertGreater(MAX_RETRIES, 0)

    def test_retry_backoff_positive(self):
        self.assertGreater(RETRY_BACKOFF, 0)


class TestDownloadFilesRetry(unittest.TestCase):
    """download_files should retry on failure and succeed if a later attempt works."""

    def setUp(self):
        self.rdams = _make_rdams()

    @patch("globsim.download.RDA.time.sleep")
    def test_retries_on_exception_then_succeeds(self, mock_sleep):
        """If the first attempt raises, second attempt should succeed without further retries."""
        call_count = {"n": 0}

        def side_effect(url, out_dir):
            call_count["n"] += 1
            if call_count["n"] == 1:
                raise requests.Timeout("stall")

        with patch.object(self.rdams, "_download_file", side_effect=side_effect):
            self.rdams.download_files(["http://example.com/file.nc"], out_dir="/tmp/")

        self.assertEqual(call_count["n"], 2)
        # One sleep between attempt 1 and 2
        mock_sleep.assert_called_once()

    @patch("globsim.download.RDA.time.sleep")
    def test_gives_up_after_max_retries(self, mock_sleep):
        """download_files should stop after MAX_RETRIES attempts and not re-raise."""
        with patch.object(
            self.rdams, "_download_file", side_effect=requests.Timeout("stall")
        ):
            # Should not raise even though all attempts fail
            self.rdams.download_files(
                ["http://example.com/file.nc"], out_dir="/tmp/", retries=3
            )

        self.assertEqual(mock_sleep.call_count, 2)  # sleep between attempt 1-2 and 2-3

    def test_no_retry_on_success(self):
        """If first attempt succeeds, _download_file is called exactly once."""
        with patch.object(self.rdams, "_download_file") as mock_dl:
            self.rdams.download_files(["http://example.com/file.nc"], out_dir="/tmp/")

        mock_dl.assert_called_once()

    @patch("globsim.download.RDA.time.sleep")
    def test_backoff_grows_exponentially(self, mock_sleep):
        """Sleep durations should be 5, 10, 20, … (RETRY_BACKOFF * 2^attempt)."""
        with patch.object(
            self.rdams, "_download_file", side_effect=Exception("err")
        ):
            self.rdams.download_files(
                ["http://example.com/file.nc"], out_dir="/tmp/", retries=4
            )

        expected_delays = [
            RETRY_BACKOFF * (2 ** i) for i in range(3)  # retries-1 sleeps
        ]
        actual_delays = [c.args[0] for c in mock_sleep.call_args_list]
        self.assertEqual(actual_delays, expected_delays)


class TestDownloadFileTimeouts(unittest.TestCase):
    """_download_file should pass (CONNECT_TIMEOUT, READ_TIMEOUT) to requests calls."""

    def setUp(self):
        self.rdams = _make_rdams()

    def _make_head_response(self, content_length=1000, accept_ranges="bytes"):
        resp = MagicMock()
        resp.status_code = 200
        resp.headers = {
            "Content-Length": str(content_length),
            "Accept-Ranges": accept_ranges,
        }
        resp.raise_for_status = MagicMock()
        return resp

    def _make_get_response(self, chunks=(b"x" * 1000,)):
        resp = MagicMock()
        resp.status_code = 200
        resp.raise_for_status = MagicMock()
        resp.iter_content = MagicMock(return_value=iter(chunks))
        resp.__enter__ = MagicMock(return_value=resp)
        resp.__exit__ = MagicMock(return_value=False)
        return resp

    @patch("globsim.download.RDA.requests.get")
    @patch("globsim.download.RDA.requests.head")
    def test_head_uses_timeout(self, mock_head, mock_get):
        mock_head.return_value = self._make_head_response()
        mock_get.return_value = self._make_get_response()

        with patch("builtins.open", unittest.mock.mock_open()):
            with patch("os.path.exists", return_value=False):
                with patch("os.stat"):
                    self.rdams._download_file(
                        "http://example.com/file.nc", "/tmp/"
                    )

        _, head_kwargs = mock_head.call_args
        self.assertEqual(head_kwargs["timeout"], (CONNECT_TIMEOUT, READ_TIMEOUT))

    @patch("globsim.download.RDA.requests.get")
    @patch("globsim.download.RDA.requests.head")
    def test_get_uses_timeout(self, mock_head, mock_get):
        mock_head.return_value = self._make_head_response()
        mock_get.return_value = self._make_get_response()

        with patch("builtins.open", unittest.mock.mock_open()):
            with patch("os.path.exists", return_value=False):
                with patch("os.stat"):
                    self.rdams._download_file(
                        "http://example.com/file.nc", "/tmp/"
                    )

        _, get_kwargs = mock_get.call_args
        self.assertEqual(get_kwargs["timeout"], (CONNECT_TIMEOUT, READ_TIMEOUT))


class TestDownloadFileSkipAndResume(unittest.TestCase):
    """_download_file should skip complete files and resume partial ones."""

    def setUp(self):
        self.rdams = _make_rdams()

    def _make_head_response(self, content_length=1000, accept_ranges="bytes"):
        resp = MagicMock()
        resp.headers = {
            "Content-Length": str(content_length),
            "Accept-Ranges": accept_ranges,
        }
        resp.raise_for_status = MagicMock()
        return resp

    @patch("globsim.download.RDA.requests.get")
    @patch("globsim.download.RDA.requests.head")
    def test_skips_complete_file(self, mock_head, mock_get):
        """No GET request should be made when local file size == Content-Length."""
        mock_head.return_value = self._make_head_response(content_length=500)

        stat_result = MagicMock()
        stat_result.st_size = 500

        with patch("os.path.exists", return_value=True):
            with patch("os.stat", return_value=stat_result):
                self.rdams._download_file("http://example.com/file.nc", "/tmp/")

        mock_get.assert_not_called()

    @patch("globsim.download.RDA.requests.get")
    @patch("globsim.download.RDA.requests.head")
    def test_resumes_partial_file(self, mock_head, mock_get):
        """GET request should include a Range header for a partial local file."""
        mock_head.return_value = self._make_head_response(
            content_length=1000, accept_ranges="bytes"
        )

        get_resp = MagicMock()
        get_resp.raise_for_status = MagicMock()
        get_resp.iter_content = MagicMock(return_value=iter([b"x" * 500]))
        mock_get.return_value = get_resp

        stat_result = MagicMock()
        stat_result.st_size = 500

        with patch("builtins.open", unittest.mock.mock_open()):
            with patch("os.path.exists", return_value=True):
                with patch("os.stat", return_value=stat_result):
                    with patch.object(self.rdams, "check_file_status"):
                        self.rdams._download_file(
                            "http://example.com/file.nc", "/tmp/"
                        )

        _, get_kwargs = mock_get.call_args
        self.assertIn("Range", get_kwargs["headers"])
        self.assertEqual(get_kwargs["headers"]["Range"], "bytes=500-")

    @patch("globsim.download.RDA.requests.get")
    @patch("globsim.download.RDA.requests.head")
    def test_no_resume_when_no_range_support(self, mock_head, mock_get):
        """When Accept-Ranges is 'none', download should restart from byte 0 (mode 'wb')."""
        mock_head.return_value = self._make_head_response(
            content_length=1000, accept_ranges="none"
        )

        get_resp = MagicMock()
        get_resp.raise_for_status = MagicMock()
        get_resp.iter_content = MagicMock(return_value=iter([b"x" * 1000]))
        mock_get.return_value = get_resp

        stat_result = MagicMock()
        stat_result.st_size = 500

        mock_open = unittest.mock.mock_open()
        with patch("builtins.open", mock_open):
            with patch("os.path.exists", return_value=True):
                with patch("os.stat", return_value=stat_result):
                    with patch.object(self.rdams, "check_file_status"):
                        self.rdams._download_file(
                            "http://example.com/file.nc", "/tmp/"
                        )

        # open should have been called with 'wb' (not 'ab')
        open_call_args = mock_open.call_args
        self.assertEqual(open_call_args[0][1], "wb")

        _, get_kwargs = mock_get.call_args
        self.assertNotIn("Range", get_kwargs.get("headers", {}))


if __name__ == "__main__":
    unittest.main()
