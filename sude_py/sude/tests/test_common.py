import unittest

from sude import SUDE, __version__, sude


class TestPublicAPI(unittest.TestCase):
    def test_public_exports(self):
        self.assertTrue(callable(sude))
        self.assertTrue(callable(SUDE))
        self.assertRegex(__version__, r"^\d+\.\d+\.\d+$")


if __name__ == "__main__":
    unittest.main()
