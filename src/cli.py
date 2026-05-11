from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import os
sys.path.append(str(Path(__file__).parent))  # 도메인 모듈 import를 위해 경로 추가

# 도메인 모듈 (프로젝트 구조에 맞게 import 경로 조정 필요)
from clinvar import ClinVarClient
from pubmed import PubMedClient
from clinvar_fetch import ClinVarFetchClient
from export import write_json, write_tsv
# [MODIFIED] 이전 단계에서 작성한 RAG 파이프라인 모듈 import 추가
from rag_pipeline.pipelines.article_pipeline import RAGDataPipeline 

def _add_common_args(p: argparse.ArgumentParser) -> None:
    # [MODIFIED] API Key는 선택사항이되 default=None 명시를 제거하여 argparse 기본 동작에 맡김
    p.add_argument("--email", required=True, help="NCBI Entrez email (required)")
    p.add_argument("--api-key", help="NCBI API key (optional)")
    p.add_argument("--debug", action="store_true", help="Print debug logs")


def _cmd_clinvar(ns: argparse.Namespace) -> int:
    # [MODIFIED] Rule 0: max_results 필수값 검증
    if ns.max_results is None or ns.max_results < 0:
        raise ValueError("max_results 값이 누락되었거나 음수입니다. 무제한 검색시 0을 명시적으로 입력하세요.")

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
    # [MODIFIED] Rule 0: sort 및 max_results 강제 예외처리
    if not ns.sort:
        raise ValueError("sort 옵션이 누락되었습니다. (예: 'relevance' 또는 'pub+date')")
    if ns.max_results is None or ns.max_results <= 0:
        raise ValueError("max_results 옵션은 1 이상이어야 합니다.")

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


def _cmd_clinvar_export(ns: argparse.Namespace) -> int:
    # [MODIFIED] Rule 0: 출력 대상, batch_size 유효성 강제 예외처리
    if not ns.out_json and not ns.out_tsv:
        raise ValueError("출력 파일 경로(--out-json 또는 --out-tsv) 중 하나 이상은 반드시 제공되어야 합니다.")
    if ns.batch_size is None or ns.batch_size <= 0:
        raise ValueError("batch_size는 1 이상의 정수여야 합니다.")

    fetcher = ClinVarFetchClient(email=ns.email, api_key=ns.api_key)

    items: list[str] = []
    if ns.input:
        p = Path(ns.input)
        raw = p.read_text(encoding="utf-8").splitlines()
        items.extend([x.strip() for x in raw if x.strip()])
    if ns.urls:
        items.extend(ns.urls)

    if not items:
        raise ValueError("--input 또는 --urls를 통해 최소 1개 이상의 변이 정보(ID/URL)를 제공해야 합니다.")

    details = fetcher.fetch_details_from_urls(items, debug=ns.debug, batch_size=ns.batch_size)

    if ns.out_json:
        out_json = Path(ns.out_json)
        write_json(details, out_json)
        print(f"json: {out_json}")
        
    if ns.out_tsv:
        out_tsv = Path(ns.out_tsv)
        write_tsv(details, out_tsv)
        print(f"tsv: {out_tsv}")

    print(f"exported: {len(details)} variants")
    return 0

# [NEW] 이전 단계에서 만든 PMC RAG 파이프라인 연동 커맨드
def _cmd_pubmed_rag(ns: argparse.Namespace) -> int:
    if ns.max_results is None or ns.max_results <= 0:
        raise ValueError("max_results 옵션은 1 이상이어야 합니다.")
        
    pipeline = RAGDataPipeline(email=ns.email, api_key=ns.api_key or "")
    
    results = pipeline.run(query=ns.query, max_results=ns.max_results)
    
    if ns.out_json:
        out_path = Path(ns.out_json)
        # Dataclass를 JSON 직렬화 가능하도록 변환
        dump_data = [
            {
                "metadata": {k: getattr(res["metadata"], k) for k in ["pmid", "pmcid", "title", "journal", "pub_date"]},
                "sections": [{"title": s.title, "content": s.content} for s in res["full_text"].sections]
            }
            for res in results
        ]
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(dump_data, f, ensure_ascii=False, indent=2)
        print(f"RAG Data exported to: {out_path}")
        
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="ncbi_search")
    sub = p.add_subparsers(dest="cmd", required=True)

    # clinvar
    p_cv = sub.add_parser("clinvar", help="Search ClinVar by disease")
    _add_common_args(p_cv)
    p_cv.add_argument("--disease", required=True)
    # [MODIFIED] default 제거, required=True 로 강제
    p_cv.add_argument("--max-results", type=int, required=True, help="0 means fetch as much as possible")
    p_cv.add_argument("--json", action="store_true")
    p_cv.set_defaults(func=_cmd_clinvar)

    # pubmed
    p_pm = sub.add_parser("pubmed", help="Search PubMed by keyword/topic")
    _add_common_args(p_pm)
    p_pm.add_argument("--query", required=True)
    # [MODIFIED] default 제거, required=True 로 강제
    p_pm.add_argument("--sort", required=True, choices=["relevance", "pub+date"])
    p_pm.add_argument("--max-results", type=int, required=True)
    p_pm.add_argument("--json", action="store_true")
    p_pm.set_defaults(func=_cmd_pubmed)

    # clinvar-export
    p_ce = sub.add_parser("clinvar-export", help="Fetch ClinVar Variant Details")
    _add_common_args(p_ce)
    p_ce.add_argument("--input", help="Text file containing ClinVar variation URLs/IDs")
    p_ce.add_argument("--urls", nargs="*", help="ClinVar variation URLs/IDs as arguments")
    # [MODIFIED] default 제거
    p_ce.add_argument("--batch-size", type=int, required=True, help="efetch batch size (e.g., 50)")
    p_ce.add_argument("--out-json", help="Output JSON path")
    p_ce.add_argument("--out-tsv", help="Output TSV path")
    p_ce.set_defaults(func=_cmd_clinvar_export)

    # [NEW] pubmed-rag (Full-text XML 파이프라인)
    p_pr = sub.add_parser("pubmed-rag", help="Fetch PubMed Open Access XMLs and parse into RAG chunks")
    _add_common_args(p_pr)
    p_pr.add_argument("--query", required=True)
    p_pr.add_argument("--max-results", type=int, required=True)
    p_pr.add_argument("--out-json", required=True, help="Output JSON path for parsed sections")
    p_pr.set_defaults(func=_cmd_pubmed_rag)

    return p


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    ns = parser.parse_args(argv)
    
    try:
        rc = ns.func(ns)
        raise SystemExit(rc)
    except Exception as e:
        # [MODIFIED] Pipeline에서 올라오는 에러를 일관되게 캐치하여 출력
        print(f"[Error] {str(e)}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()