#!/usr/bin/env python3
"""Download PDFs referenced by a spreadsheet's ``link`` column.

Only the portion of each link ending at ``doi=10.NNNN`` is used.  The script
follows redirects and, when the resulting response is HTML, follows the PDF
URL advertised by common citation metadata or an ordinary PDF link.
"""

from __future__ import annotations

import argparse
import csv
import html.parser
import re
import sys
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from xml.etree import ElementTree


LINK_RE = re.compile(r"^(.*?doi=10\.\d{4})", re.IGNORECASE)
URL_RE = re.compile(r"https?://\S+", re.IGNORECASE)
SPREADSHEET_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
USER_AGENT = "Mozilla/5.0 (compatible; replication-package-paper-downloader/1.0)"

# Links supplied for this replication package.  Keeping these in the executable
# makes the requested download reproducible without recreating a spreadsheet.
# The part after the registrant is intentionally retained here so the source
# list remains useful, although extract_link() applies the requested truncation.
PROVIDED_LINKS = """
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85086505103&doi=10.1016%2fj.jclepro.2020.122449&partnerID=40&md5=60eeb0459d6a152fcb4674f6b5bf74b8
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85086447955&doi=10.1016%2fj.envpol.2020.114833&partnerID=40&md5=f4236422dd6c4d5fdee7304d448116f3
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85084478888&doi=10.1016%2fj.envpol.2020.114751&partnerID=40&md5=68958b0ebe9cead3222a5ec06b9c19f8
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85085100834&doi=10.1016%2fj.chemosphere.2020.126984&partnerID=40&md5=bb8cac4f9e6803ea7801835b9bc3a40f
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85086457872&doi=10.1016%2fj.envint.2020.105831&partnerID=40&md5=8fc4c86cc200d18e676c1fdd1bb4f042
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85084319206&doi=10.1016%2fj.ecolind.2020.106471&partnerID=40&md5=30dc1f617069538f858f831594e862fc
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85084841045&doi=10.1007%2fs10841-020-00247-x&partnerID=40&md5=1a28070f37939790c415ed01a5202e5e
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85082869204&doi=10.1016%2fj.chemosphere.2020.126668&partnerID=40&md5=f142805c8326edf4bf64c3e1eca98e4f
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85082077073&doi=10.1111%2fina.12657&partnerID=40&md5=40014c0051423b0efd7fe9d39d03e64f
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85085060663&doi=10.1111%2fgcb.15127&partnerID=40&md5=a7cdbff1b1566b4508b9cd83f051f361
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85081027069&doi=10.1016%2fj.scitotenv.2020.137419&partnerID=40&md5=0309073281c5e85d0ae49992cdae853a
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85079840128&doi=10.1002%2feap.2074&partnerID=40&md5=9a895b0a92d1d5deda9a4b5e88a7ec7e
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85087198905&doi=10.1289%2fEHP5312&partnerID=40&md5=ebb4f16df0ee48c865e961cd1e686ac1
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85079119026&doi=10.1016%2fj.envpol.2020.114091&partnerID=40&md5=1afb371a5b9f17e25207c2a9bd5015ee
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85079873360&doi=10.1016%2fj.scitotenv.2020.137349&partnerID=40&md5=bf6153b289e4822cc62f4e392f36ffeb
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85081268314&doi=10.1016%2fj.apr.2020.02.006&partnerID=40&md5=1393646d412e0b021d3b23568a38282a
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85083046805&doi=10.1111%2fgcb.15071&partnerID=40&md5=07497a906e4a55578c73c3fae6e25f94
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85084266824&doi=10.1016%2fj.actao.2020.103577&partnerID=40&md5=386033884f36ed8f023832c74a6544e4
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85078199430&doi=10.1016%2fj.scitotenv.2020.136571&partnerID=40&md5=d0a7938f3c73dca3af5eebc43607302b
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85077951166&doi=10.1016%2fj.scitotenv.2020.136635&partnerID=40&md5=8670c16a722cdc98d57f45e2812e1ff4
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85083996711&doi=10.3390%2fijerph17082949&partnerID=40&md5=fe418487a320a4f53a5d84dc493578d6
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85083407608&doi=10.3390%2fijerph17082623&partnerID=40&md5=2b7b73bcb7b0d44958f2550ced622b10
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85081624922&doi=10.1007%2fs11356-020-08134-3&partnerID=40&md5=1f934d70fdf04ecabe1a187ed29c7792
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85083023056&doi=10.3390%2fijerph17072474&partnerID=40&md5=3b873c0f76b52cfeffdc41a412aadbcc
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85081893535&doi=10.1007%2fs00477-020-01786-0&partnerID=40&md5=143165bb12de85791e13992ed0ace56f
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85083631929&doi=10.2166%2fwh.2020.088&partnerID=40&md5=e037f5ddc29a82a2996eb4ef11b5e757
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85079052785&doi=10.1016%2fj.envint.2019.105408&partnerID=40&md5=662e311c18ebed895445d09409872df8
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85078737203&doi=10.1016%2fj.envint.2020.105500&partnerID=40&md5=5930dfa3c438d82243d8ce2163b3af59
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85082380659&doi=10.3390%2fijerph17062130&partnerID=40&md5=16d8589ee1375f28c0f8abc4890a416f
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85077340136&doi=10.1016%2fj.envpol.2019.113476&partnerID=40&md5=ee05012c3c606ca647bd34464a5a4a31
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85076250974&doi=10.1016%2fj.envpol.2019.113589&partnerID=40&md5=19d14ccc313811516dc18744c35d089c
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85077045597&doi=10.1016%2fj.ecoenv.2019.110089&partnerID=40&md5=709159718951583f40dab31759658cbb
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85082849112&doi=10.3390%2frecycling5010006&partnerID=40&md5=eeee7a6fd8baaa3f5eb389b09d76be34
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85077658560&doi=10.1016%2fj.envres.2019.109087&partnerID=40&md5=291e02b3e630008a3549963d34fd4034
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85075748165&doi=10.1002%2fldr.3478&partnerID=40&md5=4599cc3a9a010c1a55bd8a1ae3fbfb5b
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85077875182&doi=10.1002%2feap.2044&partnerID=40&md5=0b5f135f1150990d372ec5a38668ff1e
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85081031130&doi=10.1002%2fece3.6066&partnerID=40&md5=2f20d5771889ac7f29ae80ce66238a97
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85074883195&doi=10.1016%2fj.scitotenv.2019.134989&partnerID=40&md5=46f69509b6e878da12b47a9d4bc58018
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85074689873&doi=10.1016%2fj.scitotenv.2019.134985&partnerID=40&md5=0a8e013032a26fc3ea1a61711ded9ba9
https://www.scopus.com/inward/record.uri?eid=2-s2.0-85076914838&doi=10.1016%2fj.psep.2019.12.004&partnerID=40&md5=cfbf386ebfc630e7af707e7e92aca95a
"""


def extract_link(value: object) -> str | None:
    """Return a link truncated immediately after the four DOI digits."""
    match = LINK_RE.match(str(value).strip()) if value is not None else None
    return match.group(1) if match else None


def output_name(link: str) -> str:
    """Turn the extracted link into a portable, recognisable PDF filename."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", link).strip("._") + ".pdf"


def _xlsx_rows(path: Path) -> list[dict[str, str]]:
    """Read the first worksheet of an XLSX file using the standard library."""
    with zipfile.ZipFile(path) as archive:
        try:
            strings_xml = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            strings = [
                "".join(node.text or "" for node in item.iter(f"{SPREADSHEET_NS}t"))
                for item in strings_xml.iter(f"{SPREADSHEET_NS}si")
            ]
        except KeyError:
            strings = []

        workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
        relationships = ElementTree.fromstring(
            archive.read("xl/_rels/workbook.xml.rels")
        )
        relationship_id = next(workbook.iter(f"{SPREADSHEET_NS}sheet")).attrib[
            "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
        ]
        target = next(
            node.attrib["Target"]
            for node in relationships
            if node.attrib["Id"] == relationship_id
        )
        worksheet_path = target.lstrip("/")
        if not worksheet_path.startswith("xl/"):
            worksheet_path = f"xl/{worksheet_path}"
        worksheet = ElementTree.fromstring(archive.read(worksheet_path))

    table: list[list[str]] = []
    for row in worksheet.iter(f"{SPREADSHEET_NS}row"):
        cells: dict[int, str] = {}
        for cell in row.findall(f"{SPREADSHEET_NS}c"):
            column_letters = re.match(r"[A-Z]+", cell.attrib["r"])
            if not column_letters:
                continue
            column = 0
            for letter in column_letters.group():
                column = column * 26 + ord(letter) - ord("A") + 1
            value_node = cell.find(f"{SPREADSHEET_NS}v")
            inline_node = cell.find(f".//{SPREADSHEET_NS}t")
            value = "" if value_node is None else (value_node.text or "")
            if cell.attrib.get("t") == "s" and value:
                value = strings[int(value)]
            elif cell.attrib.get("t") == "inlineStr" and inline_node is not None:
                value = inline_node.text or ""
            cells[column - 1] = value
        if cells:
            table.append([cells.get(index, "") for index in range(max(cells) + 1)])

    if not table:
        return []
    headers = table[0]
    return [dict(zip(headers, row)) for row in table[1:]]


def read_links(path: Path) -> list[str]:
    if path.suffix.lower() in {".txt", ".list"}:
        return read_link_text(path.read_text(encoding="utf-8"))
    if path.suffix.lower() == ".csv":
        with path.open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream))
    elif path.suffix.lower() == ".xlsx":
        rows = _xlsx_rows(path)
    else:
        raise ValueError("input must be a .csv, .xlsx, .txt, or .list file")

    if not rows:
        return []
    link_column = next((name for name in rows[0] if name.lower() == "link"), None)
    if link_column is None:
        raise ValueError("input does not contain a column named 'link'")
    return [link for row in rows if (link := extract_link(row.get(link_column)))]


def read_link_text(text: str) -> list[str]:
    """Extract and truncate URLs pasted directly into a plain-text file/stdin."""
    return [link for value in URL_RE.findall(text) if (link := extract_link(value))]


class _PdfLinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.metadata_url: str | None = None
        self.link_url: str | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = {key.lower(): value for key, value in attrs}
        if tag.lower() == "meta" and attributes.get("name", "").lower() in {
            "citation_pdf_url",
            "wkhealth_pdf_url",
        }:
            self.metadata_url = attributes.get("content")
        if tag.lower() == "a" and self.link_url is None:
            href = attributes.get("href")
            if href and re.search(r"(?:\.pdf(?:[?#]|$)|/pdf(?:[/?#]|$))", href, re.I):
                self.link_url = href


def _request(url: str, timeout: float):
    return urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": USER_AGENT}), timeout=timeout
    )


def download_pdf(link: str, destination: Path, timeout: float = 30) -> None:
    """Follow *link*, locate a PDF if necessary, and atomically save it."""
    with _request(link, timeout) as response:
        final_url = response.geturl()
        body = response.read()

    if not body.startswith(b"%PDF-"):
        parser = _PdfLinkParser()
        parser.feed(body.decode("utf-8", errors="replace"))
        candidate = parser.metadata_url or parser.link_url
        if not candidate:
            raise RuntimeError(f"no PDF link was found on {final_url}")
        pdf_url = urllib.parse.urljoin(final_url, candidate)
        with _request(pdf_url, timeout) as response:
            body = response.read()
        if not body.startswith(b"%PDF-"):
            raise RuntimeError(f"the PDF link did not return a PDF: {pdf_url}")

    temporary = destination.with_suffix(destination.suffix + ".part")
    temporary.write_bytes(body)
    temporary.replace(destination)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        help="CSV/XLSX link table or plain-text URL list; omit to use bundled links",
    )
    parser.add_argument("--output", type=Path, default=Path("data/papers"))
    parser.add_argument("--timeout", type=float, default=30)
    args = parser.parse_args(argv)

    args.output.mkdir(parents=True, exist_ok=True)
    try:
        if args.input:
            links = read_links(args.input)
        else:
            pasted = sys.stdin.read() if not sys.stdin.isatty() else ""
            links = read_link_text(pasted or PROVIDED_LINKS)
    except (OSError, ValueError, zipfile.BadZipFile) as error:
        parser.error(str(error))

    failures = 0
    for link in dict.fromkeys(links):
        destination = args.output / output_name(link)
        if destination.exists():
            print(f"SKIP {destination}")
            continue
        try:
            download_pdf(link, destination, args.timeout)
            print(f"SAVED {destination}")
        except Exception as error:  # Keep processing the other spreadsheet rows.
            failures += 1
            print(f"FAILED {link}: {error}", file=sys.stderr)
    print(f"Processed {len(set(links))} unique links; {failures} failed.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
