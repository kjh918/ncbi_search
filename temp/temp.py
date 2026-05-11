import requests
from bs4 import BeautifulSoup
from typing import Dict, Any

# [MODIFIED] Added function to strictly fetch the raw HTML from the ClinVar URL
def fetch_rcv_webpage(url: str) -> str:
    # [Rule 0] Strict parameter validation
    if not url:
        raise ValueError("The 'url' parameter is strictly required.")
    if not url.startswith("http"):
        raise ValueError(f"Invalid URL format provided: {url}")

    try:
        # Explicit timeout is enforced. No default assumptions.
        response = requests.get(url, timeout=15.0)
        response.raise_for_status()
        return response.text
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Failed to fetch the ClinVar RCV webpage: {e}")

# [MODIFIED] Added function to parse specific DOM elements (Summary data) from the HTML
def parse_rcv_summary_data(html_content: str) -> Dict[str, Any]:
    # [Rule 0] Strict parameter validation
    if not html_content:
        raise ValueError("The 'html_content' parameter is strictly required.")

    soup = BeautifulSoup(html_content, "html.parser")
    parsed_data = {}

    # 1. Parse Title
    title_tag = soup.find("h1")
    if not title_tag:
        raise RuntimeError("Failed to parse the page title. The ClinVar DOM structure might have changed.")
    parsed_data["title"] = title_tag.text.strip()

    # 2. Parse Clinical Significance
    # ClinVar usually structures this as <dt>Clinical significance:</dt> followed by a <dd>
    clin_sig_dt = soup.find("dt", string=lambda text: text and "Clinical significance:" in text)
    if clin_sig_dt:
        clin_sig_dd = clin_sig_dt.find_next_sibling("dd")
        parsed_data["clinical_significance"] = clin_sig_dd.text.strip() if clin_sig_dd else "Not found"
    else:
        parsed_data["clinical_significance"] = "Not found"

    # 3. Parse Condition
    condition_dt = soup.find("dt", string=lambda text: text and "Condition" in text)
    if condition_dt:
        condition_dd = condition_dt.find_next_sibling("dd")
        parsed_data["condition"] = condition_dd.text.strip() if condition_dd else "Not found"
    else:
        parsed_data["condition"] = "Not found"

    # 4. Parse Review Status
    review_dt = soup.find("dt", string=lambda text: text and "Review status:" in text)
    if review_dt:
        review_dd = review_dt.find_next_sibling("dd")
        parsed_data["review_status"] = review_dd.text.strip() if review_dd else "Not found"
    else:
        parsed_data["review_status"] = "Not found"

    return parsed_data

target_url = "https://www.ncbi.nlm.nih.gov/clinvar/RCV005861489.1/"

# 1. Fetch HTML
raw_html = fetch_rcv_webpage(url=target_url)

# 2. Parse Data
rcv_data = parse_rcv_summary_data(html_content=raw_html)

# 3. Output
print(raw_html)
# {'title': 'NM_000140.4(PAH):c.1066-11G>A AND Phenylketonuria', 'clinical_significance': 'Pathogenic', ...}