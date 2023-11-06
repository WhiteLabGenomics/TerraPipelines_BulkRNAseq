from io import StringIO
import sys
import pytest

from src.hello import print_name


class TestPrintName:
    """Object to collect the relevant tests for the function `print_name` in
    the hello module."""

    def test_print_empty(self):
        """Test for the empty case and the right exception"""
        with pytest.raises(ValueError):
            print_name("")

    def test_print_name_type(self):
        """Test for the usual case type"""
        assert print_name("Oscar") is None

    def test_print_name(self):
        """Test the usual case"""
        expected = "The name is test\n"
        capturedOutput = StringIO()  # Create StringIO object
        sys.stdout = capturedOutput  #  and redirect stdout.
        print_name("test")  # Call unchanged function.
        sys.stdout = sys.__stdout__  # Reset redirect.
        assert expected == capturedOutput.getvalue()  # Now works as before.

    def test_print_name_error(self):
        """Test whether it raises the right error"""

        with pytest.raises(
            ValueError, match=r"Name given is empty, please enter a valid name."
        ):
            print_name("")
