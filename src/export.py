from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable

from .models import ClinVarVariantDetails


TSV_COLUMNS = [
    "variation_id",
    "accession",
    "accession_version",
    "variant_type",
    "title",
    "gene_symbol",
    "gene_id",
    "hgvs_genomic",
    "hgvs_cdna",
    "hgvs_protein",
    "condition_name",
    "clinical_significance",
    "review_status",
    "last_evaluated",
    "canonical_spdi",
    "url",
]


def write_json(items: Iterable[ClinVarVariantDetails], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    data = [it.to_dict() for it in items]
    out_path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def write_tsv(items: Iterable[ClinVarVariantDetails], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=TSV_COLUMNS, delimiter="\t")
        w.writeheader()
        for it in items:
            d = it.to_dict()
            w.writerow({k: d.get(k) for k in TSV_COLUMNS})
