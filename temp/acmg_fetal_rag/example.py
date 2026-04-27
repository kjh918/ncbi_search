import requests
import json
import sys

# ---------------------------------------------------------
# 1. 예외 처리 클래스
# ---------------------------------------------------------
class ParsingError(Exception):
    pass

# ---------------------------------------------------------
# 2. ClinVar API 클라이언트 (이전 설계 적용)
# ---------------------------------------------------------
class ClinVarAPIClient:
    def __init__(self, email: str):
        if not email or not isinstance(email, str):
            raise ValueError("NCBI API 사용을 위해 email은 필수 입력값입니다.")
        self.email = email
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def fetch_variant_clinical_significance(self, gene: str, variant: str) -> dict:
        if not gene.strip() or not variant.strip():
            raise ValueError("gene과 variant는 필수값입니다.")

        # ESearch: ClinVar ID 검색
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

        variant_id = id_list[0]

        # ESummary: 임상 데이터 추출
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

# ---------------------------------------------------------
# 3. RAG 프롬프트 조립 실행부 (Mock Pipeline)
# ---------------------------------------------------------
def generate_rag_prompt(gene: str, variant: str, user_email: str) -> str:
    if not gene or not variant or not user_email:
        raise ValueError("gene, variant, user_email 값은 필수입니다.")

    print(f"🔄 외부 API 연동 중: ClinVar에서 {gene} {variant} 데이터 조회...")
    
    # API 클라이언트 초기화 및 호출
    clinvar_client = ClinVarAPIClient(email=user_email)
    try:
        clinvar_data = clinvar_client.fetch_variant_clinical_significance(gene=gene, variant=variant)
        print("✅ ClinVar 데이터 조회 성공\n")
    except Exception as e:
        print(f"❌ 오류 발생: {e}")
        sys.exit(1)

    # (가정) 로컬 벡터 DB에서 검색된 문헌 데이터
    mock_retrieved_literature = (
        "[PMID: 29165669]\n"
        "Mutations in the FGFR3 gene, specifically the c.1138G>A transition resulting in a Gly380Arg substitution, "
        "are responsible for more than 98% of achondroplasia cases. Functional assays indicate this leads to a gain-of-function."
    )

    # 최종 LLM 주입용 프롬프트 조립
    system_prompt = (
        "<|begin_of_text|><|start_header_id|>system<|end_header_id|>\n"
        "You are an expert clinical molecular geneticist. Evaluate the variant based on the provided literature and database facts.\n"
    )
    
    context_section = (
        f"--- Literature Context ---\n{mock_retrieved_literature}\n\n"
        f"--- ClinVar Database Fact ---\n"
        f"ClinVar ID: {clinvar_data['clinvar_id']}\n"
        f"Current Consensus: {clinvar_data['clinical_significance']}\n"
        f"Reference URL: {clinvar_data['source_url']}\n"
    )
    
    user_prompt = (
        "<|start_header_id|>user<|end_header_id|>\n"
        f"Based on the context, what is the clinical significance of {gene} {variant} and why?<|eot_id|>\n"
        "<|start_header_id|>assistant<|end_header_id|>\n"
    )

    final_prompt = system_prompt + context_section + user_prompt
    return final_prompt

# ---------------------------------------------------------
# 4. 스크립트 실행
# ---------------------------------------------------------
if __name__ == "__main__":
    # 실행 시 본인의 이메일로 변경하여 사용하십시오 (NCBI API 정책)
    TEST_EMAIL = "your.email@example.com" 
    TARGET_GENE = "FGFR3"
    TARGET_VARIANT = "c.1138G>A"
    
    try:
        prompt_result = generate_rag_prompt(
            gene=TARGET_GENE, 
            variant=TARGET_VARIANT, 
            user_email=TEST_EMAIL
        )
        
        print("="*60)
        print("🧠 [LLM에 전달될 최종 프롬프트]")
        print("="*60)
        print(prompt_result)
        
    except ValueError as ve:
        print(f"입력값 오류: {ve}")