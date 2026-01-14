# script.py
from __future__ import annotations
from pathlib import Path

FILES: dict[str, str] = {}

FILES["pyproject.toml"] = """\
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "ncbi-search"
version = "0.2.0"
description = "NCBI Entrez search utilities (ClinVar + PubMed)"
readme = "README.md"
requires-python = ">=3.10"
dependencies = [
  "biopython>=1.81",
]

[project.scripts]
ncbi_search = "ncbi_search.cli:main"  # ✅ CHANGED: underscore command

[tool.setuptools]
package-dir = {"" = "src"}

[tool.setuptools.packages.find]
where = ["src"]
"""

FILES["README.md"] = """\
# ncbi-search

Reusable NCBI Entrez search utilities.
- ClinVar: disease-based search
- PubMed: keyword/topic search

## Install (python3.11)
pip3.11 install -e .

## CLI
ncbi_search clinvar --email you@example.com --disease "Noonan syndrome" --max-results 5
ncbi_search pubmed  --email you@example.com --query "Noonan syndrome RASopathy" --max-results 5

## Python API
from ncbi_search import ClinVarClient, PubMedClient
"""

FILES["src/ncbi_search/__init__.py"] = """\
# ✅ CHANGED: PubMedClient export
from .clinvar import ClinVarClient
from .pubmed import PubMedClient

__all__ = ["ClinVarClient", "PubMedClient"]
"""

FILES["src/ncbi_search/models.py"] = """\
from __future__ import annotations
from dataclasses import dataclass

@dataclass(frozen=True)
class ClinVarRecord:
    title: str
    gene_symbol: str | None
    clinical_significance: str | None
    variation_id: str | None

    @property
    def url(self) -> str | None:
        if not self.variation_id:
            return None
        return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{self.variation_id}/"


# ✅ CHANGED: PubMed article model 추가
@dataclass(frozen=True)
class PubMedArticle:
    pmid: str
    title: str | None
    journal: str | None
    pub_date: str | None

    @property
    def url(self) -> str:
        return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"
"""

FILES["src/ncbi_search/utils.py"] = """\
from __future__ import annotations
import time
import random

def normalize_query(q: str) -> str:
    return " ".join((q or "").strip().split())

def polite_sleep(base: float = 0.34, jitter: float = 0.15) -> None:
    # NCBI E-utilities는 과도한 호출을 피하는 것이 권장됩니다.
    time.sleep(max(0.0, base + random.uniform(-jitter, jitter)))
"""

FILES["src/ncbi_search/clinvar.py"] = """\
from __future__ import annotations
from typing import Any, Iterable
from Bio import Entrez

from .models import ClinVarRecord
from .utils import normalize_query, polite_sleep


class ClinVarClient:
    def __init__(
        self,
        *,
        email: str,
        api_key: str | None = None,
        tool: str = "ncbi-search",
        max_retries: int = 3,
    ) -> None:
        if not email or "@" not in email:
            raise ValueError("Entrez.email 은 필수이며 유효한 이메일이어야 합니다.")

        Entrez.email = email
        Entrez.tool = tool
        if api_key:
            Entrez.api_key = api_key

        self.max_retries = max_retries

    def search_by_disease(
        self,
        disease_name: str,
        *,
        max_results: int = 10,
    ) -> list[ClinVarRecord]:
        disease_name = normalize_query(disease_name)
        if not disease_name:
            return []

        ids = self._esearch_ids(
            term=f"{disease_name}[disease]",
            retmax=max_results,
        )
        if not ids:
            return []

        docs = self._esummary(ids)
        return [self._to_record(d) for d in docs]

    def _esearch_ids(self, *, term: str, retmax: int) -> list[str]:
        last_err: Exception | None = None
        for attempt in range(1, self.max_retries + 1):
            try:
                polite_sleep()
                with Entrez.esearch(db="clinvar", term=term, retmax=retmax) as h:
                    res = Entrez.read(h)
                return list(res.get("IdList", []))
            except Exception as e:
                last_err = e
                polite_sleep(base=0.8 * attempt)
        raise RuntimeError("ClinVar ESearch 실패") from last_err

    def _esummary(self, ids: Iterable[str]) -> list[dict[str, Any]]:
        last_err: Exception | None = None
        ids_str = ",".join(ids)
        for attempt in range(1, self.max_retries + 1):
            try:
                polite_sleep()
                with Entrez.esummary(db="clinvar", id=ids_str) as h:
                    res = Entrez.read(h)
                return list(res["DocumentSummarySet"]["DocumentSummary"])
            except Exception as e:
                last_err = e
                polite_sleep(base=0.8 * attempt)
        raise RuntimeError("ClinVar ESummary 실패") from last_err

    def _to_record(self, doc: dict[str, Any]) -> ClinVarRecord:
        title = str(doc.get("title") or "")

        gene_symbol = None
        genes = doc.get("genes") or []
        if genes and isinstance(genes, list) and isinstance(genes[0], dict):
            gene_symbol = genes[0].get("symbol")

        clin_sig = None
        cs = doc.get("clinical_significance") or {}
        if isinstance(cs, dict):
            clin_sig = cs.get("description")

        variation_id = None
        vs = doc.get("variation_set") or []
        if vs and isinstance(vs, list) and isinstance(vs[0], dict):
            variation_id = str(vs[0].get("variation_id") or "") or None

        return ClinVarRecord(
            title=title,
            gene_symbol=gene_symbol,
            clinical_significance=clin_sig,
            variation_id=variation_id,
        )
"""

FILES["src/ncbi_search/pubmed.py"] = """\
# ✅ CHANGED: new module - PubMed keyword search
from __future__ import annotations
from typing import Any, Iterable

from Bio import Entrez

from .models import PubMedArticle
from .utils import normalize_query, polite_sleep


class PubMedClient:
    def __init__(
        self,
        *,
        email: str,
        api_key: str | None = None,
        tool: str = "ncbi-search",
        max_retries: int = 3,
    ) -> None:
        if not email or "@" not in email:
            raise ValueError("Entrez.email 은 필수이며 유효한 이메일이어야 합니다.")

        Entrez.email = email
        Entrez.tool = tool
        if api_key:
            Entrez.api_key = api_key

        self.max_retries = max_retries

    def search(
        self,
        query: str,
        *,
        max_results: int = 10,
        sort: str = "relevance",
    ) -> list[PubMedArticle]:
        query = normalize_query(query)
        if not query:
            return []

        pmids = self._esearch_pmids(term=query, retmax=max_results, sort=sort)
        if not pmids:
            return []

        docs = self._esummary(pmids)
        return [self._to_article(d) for d in docs]

    def _esearch_pmids(self, *, term: str, retmax: int, sort: str) -> list[str]:
        last_err: Exception | None = None
        for attempt in range(1, self.max_retries + 1):
            try:
                polite_sleep()
                with Entrez.esearch(db="pubmed", term=term, retmax=retmax, sort=sort) as h:
                    res = Entrez.read(h)
                return list(res.get("IdList", []))
            except Exception as e:
                last_err = e
                polite_sleep(base=0.8 * attempt)
        raise RuntimeError("PubMed ESearch 실패") from last_err

    def _esummary(self, pmids: Iterable[str]) -> list[dict[str, Any]]:
        last_err: Exception | None = None
        ids_str = ",".join(pmids)
        for attempt in range(1, self.max_retries + 1):
            try:
                polite_sleep()
                with Entrez.esummary(db="pubmed", id=ids_str, retmode="xml") as h:
                    res = Entrez.read(h)
                # Entrez.read(esummary pubmed)는 보통 list[dict] 형태로 반환
                return list(res)
            except Exception as e:
                last_err = e
                polite_sleep(base=0.8 * attempt)
        raise RuntimeError("PubMed ESummary 실패") from last_err

    def _to_article(self, doc: dict[str, Any]) -> PubMedArticle:
        # PubMed ESummary key들은 케이스가 섞일 수 있어 안전하게 처리
        pmid = str(doc.get("Id") or doc.get("uid") or "")

        title = doc.get("Title") or doc.get("title")
        journal = doc.get("FullJournalName") or doc.get("fulljournalname") or doc.get("Source") or doc.get("source")

        pub_date = (
            doc.get("PubDate")
            or doc.get("pubdate")
            or doc.get("EPubDate")
            or doc.get("epubdate")
        )

        return PubMedArticle(
            pmid=pmid,
            title=title,
            journal=journal,
            pub_date=pub_date,
        )
"""

FILES["src/ncbi_search/cli.py"] = """\
from __future__ import annotations
import argparse
import json

from .clinvar import ClinVarClient
from .pubmed import PubMedClient  # ✅ CHANGED


def main() -> None:
    # ✅ CHANGED: subcommands 구조로 변경 (clinvar / pubmed)
    p = argparse.ArgumentParser("ncbi_search")
    sub = p.add_subparsers(dest="cmd", required=True)

    # ---- ClinVar ----
    p_cv = sub.add_parser("clinvar", help="Search ClinVar by disease")
    p_cv.add_argument("--email", required=True)
    p_cv.add_argument("--api-key", default=None)
    p_cv.add_argument("--disease", required=True)
    p_cv.add_argument("--max-results", type=int, default=10)
    p_cv.add_argument("--json", action="store_true")

    # ---- PubMed ----
    p_pm = sub.add_parser("pubmed", help="Search PubMed by keyword/topic")
    p_pm.add_argument("--email", required=True)
    p_pm.add_argument("--api-key", default=None)
    p_pm.add_argument("--query", required=True, help="keyword/topic query (e.g. 'Noonan syndrome RASopathy')")
    p_pm.add_argument("--max-results", type=int, default=10)
    p_pm.add_argument("--sort", default="relevance", choices=["relevance", "pub+date"], help="PubMed sort")
    p_pm.add_argument("--json", action="store_true")

    args = p.parse_args()

    if args.cmd == "clinvar":
        client = ClinVarClient(email=args.email, api_key=args.api_key)
        records = client.search_by_disease(args.disease, max_results=args.max_results)

        if args.json:
            print(json.dumps([r.__dict__ | {"url": r.url} for r in records], ensure_ascii=False, indent=2))
            return

        for i, r in enumerate(records, start=1):
            print("-" * 30)
            print(f"결과 {i}")
            print(f"Title: {r.title}")
            print(f"Gene: {r.gene_symbol or 'N/A'}")
            print(f"Clinical Significance: {r.clinical_significance or 'N/A'}")
            print(f"URL: {r.url or 'N/A'}")
        return

    if args.cmd == "pubmed":
        client = PubMedClient(email=args.email, api_key=args.api_key)
        articles = client.search(args.query, max_results=args.max_results, sort=args.sort)

        if args.json:
            print(json.dumps([a.__dict__ | {"url": a.url} for a in articles], ensure_ascii=False, indent=2))
            return

        for i, a in enumerate(articles, start=1):
            print("-" * 30)
            print(f"결과 {i}")
            print(f"PMID: {a.pmid}")
            print(f"Title: {a.title or 'N/A'}")
            print(f"Journal: {a.journal or 'N/A'}")
            print(f"PubDate: {a.pub_date or 'N/A'}")
            print(f"URL: {a.url}")
        return
"""

FILES["tests/test_smoke.py"] = """\
def test_import():
    import ncbi_search  # noqa: F401
"""

def main() -> None:
    root = Path("ncbi_search_pkg")
    for rel, content in FILES.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(content, encoding="utf-8")
    print(f"✅ ncbi_search_pkg 생성 완료: {root.resolve()}")

if __name__ == "__main__":
    main()
