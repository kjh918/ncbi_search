import requests
import xml.etree.ElementTree as ET
from typing import Dict, List
from core.exceptions import ParsingError

class ClinVarAPIClient:
    """
    NCBI E-utilities를 통해 ClinVar에서 특정 변이의 ACMG 분류 결과 및 증거를 추출하는 클라이언트
    """
    def __init__(self, email: str):
        if not email or not isinstance(email, str):
            raise ValueError("NCBI API 사용을 위해 email은 필수 입력값입니다.")
        self.email = email
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def fetch_variant_clinical_significance(self, gene: str, variant: str) -> Dict[str, str]:
        """
        유전자와 변이명을 입력받아 ClinVar의 요약된 임상적 의의를 반환합니다.
        
        Args:
            gene (str): 유전자 심볼 (예: 'FGFR3')
            variant (str): 변이명 (예: 'c.1138G>A')
            
        Returns:
            Dict: Variant ID, Clinical Significance 등이 포함된 딕셔너리
        """
        if not gene.strip() or not variant.strip():
            raise ValueError("gene과 variant는 필수값입니다.")

        # 1. ESearch: 변이 검색을 통해 ClinVar ID 획득
        search_term = f"{gene}[Gene] AND {variant}[Variant Name]"
        search_url = f"{self.base_url}/esearch.fcgi"
        search_params = {
            "db": "clinvar",
            "term": search_term,
            "retmode": "json",
            "email": self.email
        }
        
        try:
            search_res = requests.get(search_url, params=search_params, timeout=10)
            search_res.raise_for_status()
            search_data = search_res.json()
        except requests.RequestException as e:
            raise ParsingError(f"ClinVar 검색 API 호출 중 오류 발생: {str(e)}")

        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            raise ParsingError(f"ClinVar에서 해당 변이({gene} {variant})를 찾을 수 없습니다.")

        variant_id = id_list[0] # 가장 관련성 높은 첫 번째 ID 사용

        # 2. ESummary: 해당 ID의 임상 데이터 추출
        summary_url = f"{self.base_url}/esummary.fcgi"
        summary_params = {
            "db": "clinvar",
            "id": variant_id,
            "retmode": "json",
            "email": self.email
        }

        try:
            summary_res = requests.get(summary_url, params=summary_params, timeout=10)
            summary_res.raise_for_status()
            summary_data = summary_res.json()
        except requests.RequestException as e:
            raise ParsingError(f"ClinVar Summary API 호출 중 오류 발생: {str(e)}")

        try:
            # JSON 응답에서 핵심 임상 판정 결과(Pathogenic 등) 추출
            uid_data = summary_data["result"][variant_id]
            clinical_significance = uid_data.get("clinical_significance", {}).get("description", "Unknown")
            
            return {
                "clinvar_id": variant_id,
                "gene": gene,
                "variant": variant,
                "clinical_significance": clinical_significance,
                "source_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variant_id}/"
            }
        except KeyError as e:
            raise ParsingError(f"ClinVar 데이터 구조 분석 실패: {str(e)}")