import sys
import time
import json
import argparse
import urllib.parse
import xml.etree.ElementTree as ET
from typing import List, Dict, Any
from Bio import Entrez
import xmltodict

def convert_element_to_dict(doc: ET.Element) -> dict:
    if doc is None:
        raise ValueError("The 'doc' element must be explicitly provided. None is not allowed.")
        
    try:
        # Convert Element object to string
        doc_xml_str = ET.tostring(doc, encoding='unicode')
    except Exception as e:
        raise RuntimeError(f"Failed to convert XML element to string: {e}")

    if not doc_xml_str.strip():
        raise ValueError("The converted XML string is empty.")

    try:
        # [MODIFIED] Convert the raw XML string directly into a Python dictionary
        doc_dict = xmltodict.parse(doc_xml_str)
        return doc_dict
    except Exception as e:
        raise RuntimeError(f"Failed to parse XML string into dictionary using xmltodict: {e}")

def convert_element_to_dict(doc: ET.Element) -> dict:
    if doc is None:
        raise ValueError("The 'doc' element must be explicitly provided. None is not allowed.")
        
    try:
        # Convert Element object to string
        doc_xml_str = ET.tostring(doc, encoding='unicode')
    except Exception as e:
        raise RuntimeError(f"Failed to convert XML element to string: {e}")

    if not doc_xml_str.strip():
        raise ValueError("The converted XML string is empty.")

    try:
        # [MODIFIED] Convert the raw XML string directly into a Python dictionary
        doc_dict = xmltodict.parse(doc_xml_str)
        return doc_dict
    except Exception as e:
        raise RuntimeError(f"Failed to parse XML string into dictionary using xmltodict: {e}")
    
# [MODIFIED] 과거 코드의 안전한 디코딩 방식 차용
def _decode_xml(payload: bytes) -> str:
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError:
        return payload.decode("latin-1", errors="replace")

def fetch_clinvar_evidence_by_disease(email: str, disease_name: str, max_results: int) -> List[Dict[str, Any]]:
    # [Rule 0] 필수 파라미터 강제 검증
    if not email: raise ValueError("email 파라미터가 누락되었습니다.")
    if not disease_name: raise ValueError("disease_name 파라미터가 누락되었습니다.")
    if max_results is None or max_results <= 0: raise ValueError("max_results는 1 이상이어야 합니다.")

    Entrez.email = email

    clean_disease_name = urllib.parse.unquote(disease_name).strip()
    if not clean_disease_name:
        raise ValueError("디코딩된 disease_name이 비어있습니다.")

    search_term = f'"{clean_disease_name}"[Condition] OR "{clean_disease_name}"'
    
    print(f"  -> [Debug] NCBI ESearch 쿼리: {search_term}")

    # 1. ESearch (단순 ID 리스트 반환이므로 Entrez.read 사용 안전)
    try:
        with Entrez.esearch(db="clinvar", term=search_term, retmax=max_results) as handle:
            record = Entrez.read(handle)
        id_list = record.get("IdList", [])
    except Exception as e:
        raise RuntimeError(f"ClinVar esearch 호출 실패: {e}")

    print(f"  -> [Debug] ESearch 확보된 Variation ID 개수: {len(id_list)}개")

    if not id_list:
        raise RuntimeError(f"검색 결과가 0건입니다. 쿼리명: '{search_term}'")

    time.sleep(0.35)

    # 2. ESummary (Raw XML 직접 추출 및 파싱 - DTD 에러 원천 차단)
    try:
        with Entrez.esummary(db="clinvar", id=",".join(id_list)) as handle:
            # [MODIFIED] Entrez.read()를 버리고 과거 코드처럼 h.read() + ET.fromstring 사용
            raw_xml_bytes = handle.read()
            xml_text = _decode_xml(raw_xml_bytes)
            root = ET.fromstring(xml_text)
    except Exception as e:
        raise RuntimeError(f"ClinVar esummary 호출 및 XML 파싱 실패: {e}")

    evidences = []
    for doc in root.findall(".//DocumentSummary"):
        # [MODIFIED] doc 객체의 전체 XML 구조를 문자열로 변환하여 출력 (디버깅용)
        # encoding='unicode'를 주어 bytes가 아닌 읽기 편한 str 형태로 출력합니다.
        doc_xml_str = ET.tostring(doc, encoding='unicode')
        print("========== [DEBUG: DocumentSummary XML] ==========")
        #print(doc_xml_str)
        dict_data = convert_element_to_dict(doc)
        print(dict_data)
        print("==================================================")
        #exit()
        var_id = doc.get("uid", "")
        
        title_elem = doc.find("title")
        title = title_elem.text if title_elem is not None else ""
        
        clin_sig = "Not Provided"
        review_status = "Not Provided"
        clin_sig_elem = doc.find("clinical_significance")
        if clin_sig_elem is not None:
            desc_elem = clin_sig_elem.find("description")
            if desc_elem is not None and desc_elem.text:
                clin_sig = desc_elem.text
            rev_elem = clin_sig_elem.find("review_status")
            if rev_elem is not None and rev_elem.text:
                review_status = rev_elem.text

        genes = []
        for gene_elem in doc.findall(".//genes/gene"):
            sym = gene_elem.get("symbol")
            if sym:
                genes.append(sym)
        
        if var_id:
            evidences.append({
                "variation_id": var_id,
                "title": title,
                "genes": genes,
                "clinical_significance": clin_sig,
                "review_status": review_status,
                "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{var_id}/"
            })


    return evidences

def main():
    parser = argparse.ArgumentParser(description="Fetch ClinVar evidence by Disease Name (Raw XML Parsing)")
    # [Rule 0] default 배제, required=True 강제
    parser.add_argument("--email", required=True, help="NCBI API 이메일")
    parser.add_argument("--disease", required=True, help="질환명 (예: Achondroplasia)")
    parser.add_argument("--max-results", type=int, required=True, help="최대 검색 건수")

    args = parser.parse_args()

    if args.max_results <= 0:
        print("[Error] max-results는 1 이상이어야 합니다.")
        sys.exit(1)

    try:
        print(f"▶ '{args.disease}' 질환에 대한 ClinVar 근거 검색 시작...")
        results = fetch_clinvar_evidence_by_disease(
            email=args.email, 
            disease_name=args.disease, 
            max_results=args.max_results
        )
        
        print(f"\n▶ 검색 완료: 총 {len(results)}건 확보")
        print(json.dumps(results, ensure_ascii=False, indent=2))
        
    except Exception as e:
        print(f"\n[Fatal Error] {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()