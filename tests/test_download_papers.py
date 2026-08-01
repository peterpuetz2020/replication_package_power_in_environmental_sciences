import csv
import tempfile
import unittest
from pathlib import Path

from scripts.download_papers import (
    PROVIDED_LINKS,
    extract_doi,
    extract_link,
    output_name,
    read_link_text,
    read_links,
)


class DownloadPapersTests(unittest.TestCase):
    def test_extracts_requested_scopus_prefix(self):
        link = "https://www.scopus.com/inward/record.uri?eid=2-s2.0-85028274653&doi=10.1016/j.foo.2017.01.001&partnerID=40"
        self.assertEqual(
            extract_link(link),
            "https://www.scopus.com/inward/record.uri?eid=2-s2.0-85028274653&doi=10.1016",
        )

    def test_rejects_missing_or_short_doi_prefix(self):
        self.assertIsNone(extract_link("https://example.test/no-doi"))
        self.assertIsNone(extract_link("https://example.test/?doi=10.123"))

    def test_extracts_complete_url_encoded_doi_for_resolution(self):
        link = "https://www.scopus.com/inward/record.uri?eid=one&doi=10.1016%2fj.envpol.2020.114833&partnerID=40"
        self.assertEqual(extract_doi(link), "10.1016/j.envpol.2020.114833")

    def test_rejects_truncated_doi_for_resolution(self):
        self.assertIsNone(extract_doi("https://x.test/?doi=10.1016"))

    def test_reads_case_insensitive_link_column_and_skips_invalid_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "links.csv"
            with path.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=["id", "Link"])
                writer.writeheader()
                writer.writerow({"id": "1", "Link": "https://x.test/?doi=10.1234/rest"})
                writer.writerow({"id": "2", "Link": "not a link"})
            self.assertEqual(read_links(path), ["https://x.test/?doi=10.1234"])

    def test_filename_is_portable(self):
        name = output_name("https://x.test/path?doi=10.1234")
        self.assertEqual(name, "https_x.test_path_doi_10.1234.pdf")

    def test_reads_urls_pasted_as_plain_text(self):
        text = """Notes before the URLs
https://www.scopus.com/inward/record.uri?eid=one&doi=10.1016%2fj.foo
https://www.scopus.com/inward/record.uri?eid=two&doi=10.3390%2fbar
"""
        self.assertEqual(
            read_link_text(text),
            [
                "https://www.scopus.com/inward/record.uri?eid=one&doi=10.1016",
                "https://www.scopus.com/inward/record.uri?eid=two&doi=10.3390",
            ],
        )

    def test_reads_plain_text_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "provided_links.txt"
            path.write_text("https://x.test/?doi=10.1289%2fEHP5312\n")
            self.assertEqual(read_links(path), ["https://x.test/?doi=10.1289"])

    def test_bundled_links_are_valid_and_unique(self):
        links = read_link_text(PROVIDED_LINKS)
        self.assertEqual(len(links), 40)
        self.assertEqual(len(links), len(set(links)))


if __name__ == "__main__":
    unittest.main()
