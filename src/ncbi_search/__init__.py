# ✅ CHANGED: PubMedClient export
from .clinvar import ClinVarClient
from .pubmed import PubMedClient
from .clinvar_fetch import ClinVarFetchClient

__all__ = ["ClinVarClient", "PubMedClient", "ClinVarFetchClient"]
