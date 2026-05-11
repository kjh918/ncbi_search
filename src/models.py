from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import Any, Optional


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


@dataclass(frozen=True)
class PubMedArticle:
    pmid: str
    title: str | None
    journal: str | None
    pub_date: str | None

    @property
    def url(self) -> str:
        return f"https://pubmed.ncbi.nlm.nih.gov/{self.pmid}/"

@dataclass(frozen=True)
class ArticleSection:
    title: str
    content: str

    def __post_init__(self):
        if not self.content:
            raise ValueError("ArticleSection: 본문 내용(content)은 필수입니다.")

# [NEW] XML Document 전체 모델
@dataclass(frozen=True)
class PMCXMLDocument:
    pmcid: str
    sections: list[ArticleSection]

    def __post_init__(self):
        if not self.pmcid:
            raise ValueError("PMCXMLDocument: pmcid는 필수입니다.")
        if not self.sections:
            raise ValueError("PMCXMLDocument: 파싱된 sections 데이터가 없습니다.")

    @property
    def url(self) -> str:
        return f"https://www.ncbi.nlm.nih.gov/pmc/articles/{self.pmcid}/"
    
# ✅ NEW: ClinVar variation details 정규화 모델 (export용)
@dataclass(frozen=True)
class ClinVarVariantDetails:
    variation_id: str
    accession: str | None
    accession_version: str | None
    variant_type: str | None
    title: str | None

    gene_symbol: str | None
    gene_id: str | None

    hgvs_genomic: str | None
    hgvs_cdna: str | None
    hgvs_protein: str | None

    condition_name: str | None
    clinical_significance: str | None
    review_status: str | None
    last_evaluated: str | None

    canonical_spdi: str | None

    url: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
