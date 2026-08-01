#!/usr/bin/env python3
"""Open every supplied paper in Chrome, click “View PDF”, and keep tabs open.

Install the browser-only dependency first with
``python -m pip install -r requirements-download-papers.txt``. Chrome's detached
mode deliberately leaves the browser and all article tabs open when this script
finishes.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

from selenium import webdriver
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as expected
from selenium.webdriver.support.ui import WebDriverWait

from download_papers import PROVIDED_LINKS, URL_RE, extract_link, output_name


VIEW_PDF_XPATH = (
    "//*[self::a or self::button]["
    "contains(translate(normalize-space(.), 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', "
    "'abcdefghijklmnopqrstuvwxyz'), 'view pdf') or "
    "contains(translate(@aria-label, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', "
    "'abcdefghijklmnopqrstuvwxyz'), 'view pdf')]"
)


def wait_for_download(directory: Path, previous: set[Path], timeout: float) -> Path:
    """Return a newly completed PDF from Chrome's download directory."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        current = set(directory.iterdir())
        partial = {path for path in current if path.suffix.lower() == ".crdownload"}
        pdfs = [path for path in current - previous if path.suffix.lower() == ".pdf"]
        if pdfs and not partial:
            return max(pdfs, key=lambda path: path.stat().st_mtime)
        time.sleep(0.5)
    raise TimeoutError("the PDF download did not finish before the timeout")


def browser_options(download_directory: Path) -> webdriver.ChromeOptions:
    options = webdriver.ChromeOptions()
    options.add_experimental_option("detach", True)
    options.add_experimental_option(
        "prefs",
        {
            "download.default_directory": str(download_directory),
            "download.prompt_for_download": False,
            "download.directory_upgrade": True,
            "plugins.always_open_pdf_externally": True,
        },
    )
    return options


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/papers"))
    parser.add_argument("--page-timeout", type=float, default=60)
    parser.add_argument("--download-timeout", type=float, default=120)
    args = parser.parse_args()

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    driver = webdriver.Chrome(options=browser_options(output))
    links = list(dict.fromkeys(URL_RE.findall(PROVIDED_LINKS)))

    failures = 0
    for index, source_url in enumerate(links):
        extracted = extract_link(source_url)
        if extracted is None:
            continue
        destination = output / output_name(extracted)
        if index:
            driver.execute_script("window.open(arguments[0], '_blank');", source_url)
            driver.switch_to.window(driver.window_handles[-1])
        else:
            driver.get(source_url)

        if destination.exists():
            print(f"SKIP {destination}")
            continue
        before = set(output.iterdir())
        try:
            view_pdf = WebDriverWait(driver, args.page_timeout).until(
                expected.element_to_be_clickable((By.XPATH, VIEW_PDF_XPATH))
            )
            driver.execute_script("arguments[0].scrollIntoView({block: 'center'});", view_pdf)
            view_pdf.click()
            downloaded = wait_for_download(output, before, args.download_timeout)
            downloaded.replace(destination)
            print(f"SAVED {destination}")
        except (TimeoutException, TimeoutError, OSError) as error:
            failures += 1
            print(f"FAILED {source_url}: {error}")

    print(f"Processed {len(links)} links; {failures} failed. Browser tabs remain open.")
    # Intentionally do not call driver.quit(): Chrome was started with detach=True.
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
