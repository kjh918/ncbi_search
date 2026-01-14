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

    # -------------------------------------------------
    # ✅ CHANGED: disease 단독(=dis/disease) 대비 결과가 줄지 않게
    # - tag별 Count 확인 → 페이지네이션으로 충분히 가져와 union 구성
    # - 마지막에 max_results로 자름
    # -------------------------------------------------
    def search_by_disease_or_condition(
        self,
        query: str,
        *,
        max_results: int = 10000,
        debug: bool = False,
        tags: tuple[str | None, ...] = ("dis", "disease", "phenotype", None),
        per_tag_fetch: int | None = None,  # ✅ CHANGED: 태그별로 얼마나 넉넉히 가져올지
        page_size: int = 500,              # ✅ CHANGED
        hard_cap_per_tag: int = 5000,      # ✅ CHANGED: 태그별 최대 수집 상한
    ) -> list[ClinVarRecord]:
        """
        tag별 top N만 가져오면 중복 때문에 union이 '줄어드는 것처럼' 보일 수 있음.
        그래서 tag별로 넉넉히(page) 가져와 union을 만들고, 마지막에 max_results로 자름.
        """
        query = normalize_query(query)
        if not query:
            return []

        # ✅ CHANGED: per_tag_fetch 기본값 = max_results * 5 (겹침 대비 넉넉히)
        if per_tag_fetch is None:
            per_tag_fetch = max(max_results * 5, max_results)

        # ✅ CHANGED: id 수집(발견 순서 유지)
        id_order: list[str] = []
        id_seen: set[str] = set()

        for tag in tags:
            term = self._build_term(query, tag)

            # ✅ CHANGED: 먼저 Count 확인 (retmax=0)
            _, meta0 = self._esearch_ids_with_meta(term=term, retmax=0, retstart=0)
            total = int(meta0.get("Count") or 0)

            if debug:
                self._print_count_debug(tag, term, meta0)

            if total <= 0:
                continue

            # ✅ CHANGED: 이 태그에서 가져올 목표량(넉넉히, 단 상한)
            target = min(total, per_tag_fetch, hard_cap_per_tag)

            # ✅ CHANGED: 페이지네이션으로 target까지 ids 확보
            fetched = 0
            retstart = 0
            while fetched < target:
                need = min(page_size, target - fetched)
                ids_page, meta = self._esearch_ids_with_meta(term=term, retmax=need, retstart=retstart)

                if debug and retstart == 0:
                    # 첫 페이지에서만 translation/ids 미리보기
                    self._print_search_debug(tag, term, ids_page, meta)

                if not ids_page:
                    break

                for _id in ids_page:
                    if _id not in id_seen:
                        id_seen.add(_id)
                        id_order.append(_id)

                fetched += len(ids_page)
                retstart += len(ids_page)

                if debug:
                    print(f"[ClinVar DEBUG] tag={tag} fetched={fetched}/{target} union={len(id_order)} retstart={retstart}")

        # ✅ CHANGED: union 결과가 없으면 종료(에러 X)
        if not id_order:
            return []

        # ✅ CHANGED: 최종 출력은 max_results로 자름
        ids_final = id_order[:max_results]

        docs = self._esummary_safe(ids_final, debug=debug)
        if docs:
            return [self._to_record(d) for d in docs]

        # esummary가 깨져도 최소한 추적 가능하도록 variation_id 기반 반환
        return [
            ClinVarRecord(
                title=f"(matched by union search; query={query})",
                gene_symbol=None,
                clinical_significance=None,
                variation_id=_id,
            )
            for _id in ids_final
        ]

    # -------------------------------------------------
    # ✅ CHANGED: disease-only도 같은 방식으로 비교 가능하게 제공
    # -------------------------------------------------
    def search_by_disease_only(
        self,
        disease: str,
        *,
        max_results: int = 10,
        debug: bool = False,
        tag: str = "dis",
    ) -> list[ClinVarRecord]:
        disease = normalize_query(disease)
        if not disease:
            return []

        term = self._build_term(disease, tag)
        ids, meta = self._esearch_ids_with_meta(term=term, retmax=max_results, retstart=0)

        if debug:
            self._print_search_debug(tag, term, ids, meta)

        if not ids:
            return []

        docs = self._esummary_safe(ids, debug=debug)
        if docs:
            return [self._to_record(d) for d in docs]

        return [
            ClinVarRecord(
                title=f"(matched by disease-only; term={term})",
                gene_symbol=None,
                clinical_significance=None,
                variation_id=_id,
            )
            for _id in ids
        ]

    def _build_term(self, q: str, tag: str | None) -> str:
        q_quoted = f"\"{q}\""
        if tag:
            return f"({q_quoted}[{tag}] OR {q}[{tag}])"
        return f"({q_quoted} OR {q})"

    # -----------------------------
    # ESearch
    # -----------------------------
    def _esearch_ids_with_meta(self, *, term: str, retmax: int, retstart: int) -> tuple[list[str], dict[str, Any]]:
        last_err: Exception | None = None
        for attempt in range(1, self.max_retries + 1):
            try:
                polite_sleep()
                with Entrez.esearch(
                    db="clinvar",
                    term=term,
                    retmax=retmax,
                    retstart=retstart,
                    usehistory="y",
                ) as h:
                    res = Entrez.read(h)

                ids = list(res.get("IdList", []))
                meta = {
                    "Count": res.get("Count"),
                    "RetMax": res.get("RetMax"),
                    "RetStart": res.get("RetStart"),
                    "QueryTranslation": res.get("QueryTranslation"),
                }
                return ids, meta
            except Exception as e:
                last_err = e
                polite_sleep(base=0.8 * attempt)
        raise RuntimeError("ClinVar ESearch 실패") from last_err

    def _print_count_debug(self, tag: str | None, term: str, meta: dict[str, Any]) -> None:
        print("=" * 60)
        print(f"[ClinVar DEBUG] tag={tag} term={term}")
        print(f"[ClinVar DEBUG] Count={meta.get('Count')} (retmax=0 preflight)")
        print(f"[ClinVar DEBUG] QueryTranslation={meta.get('QueryTranslation')}")
        print("=" * 60)

    def _print_search_debug(self, tag: str | None, term: str, ids: list[str], meta: dict[str, Any]) -> None:
        print("=" * 60)
        print(f"[ClinVar DEBUG] tag={tag} term={term}")
        print(f"[ClinVar DEBUG] Count={meta.get('Count')} RetMax={meta.get('RetMax')} RetStart={meta.get('RetStart')}")
        print(f"[ClinVar DEBUG] QueryTranslation={meta.get('QueryTranslation')}")
        print(f"[ClinVar DEBUG] IdList({len(ids)})={ids[:10]}{'...' if len(ids) > 10 else ''}")
        print("=" * 60)

    # -----------------------------
    # esummary “안전 버전”
    # -----------------------------
    def _esummary_safe(self, ids: Iterable[str], *, debug: bool = False) -> list[dict[str, Any]]:
        ids_str = ",".join(ids)
        last_err: Exception | None = None

        for attempt in range(1, self.max_retries + 1):
            try:
                polite_sleep()
                with Entrez.esummary(db="clinvar", id=ids_str) as h:
                    res = Entrez.read(h)

                if isinstance(res, dict):
                    dss = res.get("DocumentSummarySet")
                    if isinstance(dss, dict):
                        ds = dss.get("DocumentSummary")
                        if isinstance(ds, list):
                            return ds
                    if debug:
                        print(f"[ClinVar DEBUG] esummary dict keys={list(res.keys())}")

                if debug:
                    print(f"[ClinVar DEBUG] esummary unexpected type={type(res)}")
                return []

            except Exception as e:
                last_err = e
                if debug:
                    print(f"[ClinVar DEBUG] esummary failed attempt={attempt}: {e}")
                polite_sleep(base=0.8 * attempt)

        if debug and last_err:
            print(f"[ClinVar DEBUG] esummary giving up: {last_err}")
        return []

    # -----------------------------
    # summary dict -> record 변환
    # -----------------------------
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
