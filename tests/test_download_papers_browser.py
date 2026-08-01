import unittest
from pathlib import Path


class BrowserDownloadTests(unittest.TestCase):
    def test_browser_script_requests_detached_chrome(self):
        source = Path("scripts/download_papers_browser.py").read_text()
        self.assertIn('add_experimental_option("detach", True)', source)
        self.assertNotIn("driver.quit()", source.split("if __name__")[0].split("# Intentionally")[0])

if __name__ == "__main__":
    unittest.main()
