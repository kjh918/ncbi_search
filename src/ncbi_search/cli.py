from __future__ import annotations

import argparse
import json
from pathlib import Path

from .clinvar import ClinVarClient
from .pubmed import PubMedClient
from .clinvar_fetch import ClinVarFetchClient  # ✅ CHANGED
from .export import write_json, write_tsv       # ✅ CHANGED


def _add_common_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--email", required=True, help="NCBI Entrez email (required)")
    p.add_argument("--api-key", default=None, help="NCBI API key (optional)")
    p.add_argument("--debug", action="store_true", help="Print debug logs")


def _cmd_clinvar(ns: argparse.Namespace) -> int:
    client = ClinVarClient(email=ns.email, api_key=ns.api_key)
    max_results = None if ns.max_results == 0 else ns.max_results
    records = client.search_by_disease(ns.disease, max_results=max_results, debug=ns.debug)

    if ns.json:
        print(json.dumps([r.__dict__ | {"url": r.url} for r in records], ensure_ascii=False, indent=2))
        return 0

    for i, r in enumerate(records, start=1):
        print("-" * 30)
        print(f"결과 {i}")
        print(f"Title: {r.title}")
        print(f"Gene: {r.gene_symbol or 'N/A'}")
        print(f"Clinical Significance: {r.clinical_significance or 'N/A'}")
        print(f"URL: {r.url or 'N/A'}")
    return 0


def _cmd_pubmed(ns: argparse.Namespace) -> int:
    client = PubMedClient(email=ns.email, api_key=ns.api_key)
    articles = client.search(ns.query, max_results=ns.max_results, sort=ns.sort)

    if ns.json:
        print(json.dumps([a.__dict__ | {"url": a.url} for a in articles], ensure_ascii=False, indent=2))
        return 0

    for i, a in enumerate(articles, start=1):
        print("-" * 30)
        print(f"결과 {i}")
        print(f"PMID: {a.pmid}")
        print(f"Title: {a.title or 'N/A'}")
        print(f"Journal: {a.journal or 'N/A'}")
        print(f"PubDate: {a.pub_date or 'N/A'}")
        print(f"URL: {a.url}")
    return 0


# ✅ CHANGED: URL 리스트를 받아서 JSON/TSV로 export
def _cmd_clinvar_export(ns: argparse.Namespace) -> int:
    fetcher = ClinVarFetchClient(email=ns.email, api_key=ns.api_key)

    # 입력: 파일 또는 커맨드라인 리스트
    items: list[str] = []
    if ns.input:
        p = Path(ns.input)
        raw = p.read_text(encoding="utf-8").splitlines()
        items.extend([x.strip() for x in raw if x.strip()])
    if ns.urls:
        items.extend(ns.urls)

    details = fetcher.fetch_details_from_urls(items, debug=ns.debug, batch_size=ns.batch_size)

    out_json = Path(ns.out_json) if ns.out_json else None
    out_tsv = Path(ns.out_tsv) if ns.out_tsv else None

    if out_json:
        write_json(details, out_json)
    if out_tsv:
        write_tsv(details, out_tsv)

    # 기본: stdout 요약
    print(f"exported: {len(details)} variants")
    if out_json:
        print(f"json: {out_json}")
    if out_tsv:
        print(f"tsv: {out_tsv}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="ncbi_search")
    sub = p.add_subparsers(dest="cmd", required=True)

    # clinvar
    p_cv = sub.add_parser("clinvar", help="Search ClinVar by disease")
    _add_common_args(p_cv)
    p_cv.add_argument("--disease", required=True)
    p_cv.add_argument("--max-results", type=int, default=10, help="0 means fetch as much as possible (clinvar search)")
    p_cv.add_argument("--json", action="store_true")
    p_cv.set_defaults(func=_cmd_clinvar)

    # pubmed
    p_pm = sub.add_parser("pubmed", help="Search PubMed by keyword/topic")
    _add_common_args(p_pm)
    p_pm.add_argument("--query", required=True)
    p_pm.add_argument("--sort", default="relevance", choices=["relevance", "pub+date"])
    p_pm.add_argument("--max-results", type=int, default=10)
    p_pm.add_argument("--json", action="store_true")
    p_pm.set_defaults(func=_cmd_pubmed)

    # ✅ CHANGED: clinvar-export
    p_ce = sub.add_parser("clinvar-export", help="Fetch ClinVar Variant Details from variation URLs/IDs and export to JSON/TSV")
    _add_common_args(p_ce)
    p_ce.add_argument("--input", default=None, help="Text file containing ClinVar variation URLs/IDs (one per line)")
    p_ce.add_argument("--urls", nargs="*", default=None, help="ClinVar variation URLs/IDs as arguments")
    p_ce.add_argument("--batch-size", type=int, default=50, help="efetch batch size")
    p_ce.add_argument("--out-json", default="clinvar_variants.json", help="Output JSON path")
    p_ce.add_argument("--out-tsv", default="clinvar_variants.tsv", help="Output TSV path")
    p_ce.set_defaults(func=_cmd_clinvar_export)

    return p


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    ns = parser.parse_args(argv)
    rc = ns.func(ns)
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
