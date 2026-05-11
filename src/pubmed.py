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
