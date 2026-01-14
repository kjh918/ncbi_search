#!/usr/bin/env python3
# scripts/clinvar_search_export.py
from __future__ import annotations

import argparse
from pathlib import Path

from ncbi_search.clinvar import ClinVarClient
from ncbi_search.clinvar_fetch import ClinVarFetchClient
from ncbi_search.export import write_json, write_tsv


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="clinvar_search_export.py",
        description="ClinVar disease search -> fetch variant details -> export JSON/TSV",
    )
    p.add_argument("--email", required=True, help="NCBI Entrez email (required)")
    p.add_argument("--api-key", default=None, help="NCBI API key (optional)")

    p.add_argument("--disease", required=True, help="Disease/condition query for ClinVar search")
    p.add_argument("--search-max", type=int, default=50, help="How many variation results to collect from search")

    p.add_argument("--out-dir", default=".", help="Output directory")
    p.add_argument("--prefix", default="clinvar", help="Output prefix name (json/tsv)")

    p.add_argument("--batch-size", type=int, default=50, help="efetch batch size")
    p.add_argument("--debug", action="store_true", help="Debug logs")

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) ClinVar 검색 (disease-only API를 쓰고 있으면 여기서 그대로 호출)
    #    - 네가 disease/condition union을 원하면 해당 메서드로 변경하면 됨.
    search_client = ClinVarClient(email=args.email, api_key=args.api_key)

    # ✅ CHANGED: 지금까지 대화 흐름상 "disease" 검색이 메인이라 그대로 사용
    # 필요하면 네 코드에 맞춰 search_by_disease_or_condition(...)로 바꿔도 됨.
    records = search_client.search_by_disease_only(
        args.disease,
        max_results=args.search_max,
        debug=args.debug,
    )

    # 검색 결과가 0이어도 에러 없이 종료
    if not records:
        if args.debug:
            print("[PIPELINE] No search results.")
        # 빈 파일도 만들고 싶으면 아래 주석 해제
        # write_json([], out_dir / f"{args.prefix}.json")
        # write_tsv([], out_dir / f"{args.prefix}.tsv")
        return 0

    # 2) 검색 결과에서 URL/ID 수집
    urls_or_ids: list[str] = []
    for r in records:
        if r.url:
            urls_or_ids.append(r.url)
        elif r.variation_id:
            urls_or_ids.append(str(r.variation_id))

    if args.debug:
        print(f"[PIPELINE] collected urls_or_ids={len(urls_or_ids)}")

    # 3) Variant Details fetch + parse
    fetch_client = ClinVarFetchClient(email=args.email, api_key=args.api_key)
    details = fetch_client.fetch_details_from_urls(
        urls_or_ids,
        debug=args.debug,
        batch_size=args.batch_size,
    )

    # 4) Export
    out_json = out_dir / f"{args.prefix}.json"
    out_tsv = out_dir / f"{args.prefix}.tsv"

    write_json(details, out_json)
    write_tsv(details, out_tsv)

    print(f"exported: {len(details)} variants")
    print(f"json: {out_json}")
    print(f"tsv: {out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
