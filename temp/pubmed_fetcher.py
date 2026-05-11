import sys
import time
import json
import argparse
import re
import xml.etree.ElementTree as ET
from typing import List, Dict, Any
from Bio import Entrez

def _decode_xml(payload: bytes) -> str:
    try:
        return payload.decode("utf-8")
    except UnicodeDecodeError:
        return payload.decode("latin-1", errors="replace")

# [MODIFIED] 신규 추가: URL이나 버전(.1)이 포함된 문자열에서 순수 RCV Accession 추출
def parse_rcv_accession(raw_input: str) -> str:
    if not raw_input: 
        raise ValueError("RCV 입력값이 누락되었습니다.")
    
    # 정규식을 통해 RCV로 시작하는 문자열 파싱 (예: https://.../RCV005861489.1/ -> RCV005861489.1)
    match = re.search(r'(RCV\d+(?:\.\d+)?)', raw_input)
    if not match:
        raise ValueError(f"유효한 RCV Accession을 찾을 수 없습니다. 입력값: {raw_input}")
    
    return match.group(1)

def fetch_pmids_from_rcv(email: str, raw_rcv_input: str) -> List[str]:
    # [Rule 0] 필수 파라미터 강제 검증
    if not email: raise ValueError("email 파라미터가 누락되었습니다.")
    if not raw_rcv_input: raise ValueError("raw_rcv_input 파라미터가 누락되었습니다.")

    Entrez.email = email
    
    # 1. 입력값 정제 및 Base Accession 분리
    target_rcv = parse_rcv_accession(raw_rcv_input)
    base_rcv = target_rcv.split('.')[0] # XML 노드 매칭을 위해 버전 제외 (예: RCV005861489)
    print(f"  -> [Debug] Target RCV: {target_rcv} (Base: {base_rcv})")

    # 2. RCV 번호로 ClinVar 내부 Variation ID 검색 (NCBI efetch 규칙 대응)
    try:
        with Entrez.esearch(db="clinvar", term=f"{base_rcv}") as handle:
            record = Entrez.read(handle)
        var_ids = record.get("IdList", [])
        print(var_ids)
    except Exception as e:
        raise RuntimeError(f"ClinVar esearch 호출 실패 (RCV 검색): {e}")

    if not var_ids:
        raise RuntimeError(f"RCV Accession '{target_rcv}'에 맵핑되는 내부 Variation ID가 없습니다.")

    time.sleep(0.35)

    # 3. VCV 전체가 아닌 clinvarset(RCV 포함 세트) XML 다운로드
    try:
        with Entrez.efetch(db="clinvar", id=var_ids[0], rettype="clinvarset", retmode="xml") as handle:
            raw_xml_bytes = handle.read()
            
        xml_text = _decode_xml(raw_xml_bytes)
        
        # [Rule 0] Explicit exception handling for empty or malformed (HTML) NCBI responses
        if not xml_text or not xml_text.strip():
            raise RuntimeError("NCBI returned an empty response instead of valid XML.")
            
        if xml_text.strip().startswith("<!DOCTYPE html>") or "<html" in xml_text[:50].lower():
            raise RuntimeError(f"NCBI API returned an HTML error instead of XML. Preview: {xml_text[:200]}...")
            
        root = ET.fromstring(xml_text)
        
    except Exception as e:
        raise RuntimeError(f"Failed to fetch or parse clinvarset XML for Variation ID '{var_ids[0]}': {e}")

    # 4. 다운로드된 XML에서 정확히 타겟 RCV에 해당하는 노드만 찾아 논문 추출
    pmids = set()
    
    # ReferenceClinVarAssertion 또는 ClinVarAssertion 노드 순회
    for assertion in root.findall(".//ReferenceClinVarAssertion") + root.findall(".//ClinVarAssertion"):
        accession_node = assertion.find(".//ClinVarAccession")
        if accession_node is not None:
            acc = accession_node.get("Accession", "")
            
            # [MODIFIED] 해당 노드가 우리가 찾고자 하는 RCV인지 검증
            if acc == base_rcv:
                # 해당 RCV 노드 하위에 있는 Citation의 PubMed ID만 추출
                for citation in assertion.findall(".//Citation"):
                    for id_node in citation.findall(".//ID"):
                        if id_node.get("Source") == "PubMed" and id_node.text:
                            pmids.add(id_node.text.strip())

    return list(pmids)

def fetch_pubmed_details(email: str, pmids: List[str]) -> List[Dict[str, Any]]:
    # [Rule 0]
    if not email: raise ValueError("email 파라미터가 누락되었습니다.")
    if not pmids: raise ValueError("조회할 pmids 리스트가 누락되었거나 비어있습니다.")

    Entrez.email = email
    time.sleep(0.35)

    try:
        with Entrez.esummary(db="pubmed", id=",".join(pmids)) as handle:
            summaries = Entrez.read(handle)
    except Exception as e:
        raise RuntimeError(f"PubMed esummary 호출 실패: {e}")

    papers = []
    for summary in summaries:
        pmid = str(summary.get("Id", ""))
        title = summary.get("Title", "")
        journal = summary.get("Source", "")
        
        if pmid and title:
            papers.append({
                "pmid": pmid,
                "title": title,
                "journal": journal,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            })

    return papers

def get_evidence_papers(email: str, raw_rcv_input: str) -> List[Dict[str, Any]]:
    if not email: raise ValueError("email 필수")
    if not raw_rcv_input: raise ValueError("raw_rcv_input 필수")
    
    print(f"▶ 1. 입력된 RCV ('{raw_rcv_input}')에서 전용 PubMed ID 추출 중...")
    pmids = fetch_pmids_from_rcv(email=email, raw_rcv_input=raw_rcv_input)
    
    if not pmids:
        print("  -> 해당 RCV에 명시된 PubMed 논문이 존재하지 않거나, 추출하지 못했습니다.")
        return []
        
    print(f"▶ 2. {len(pmids)}개의 논문 상세 정보 조회 중...")
    return fetch_pubmed_details(email=email, pmids=pmids)

def main():
    parser = argparse.ArgumentParser(description="Extract PubMed Evidence strictly for a specific ClinVar RCV")
    # [Rule 0] Default 배제
    parser.add_argument("--email", required=True, help="NCBI API 이메일")
    # [MODIFIED] URL 입력도 허용하도록 help 텍스트 수정
    parser.add_argument("--rcv", required=True, help="ClinVar RCV Accession 또는 URL (예: https://www.ncbi.nlm.nih.gov/clinvar/RCV005861489.1/)")

    args = parser.parse_args()

    try:
        results = get_evidence_papers(email=args.email, raw_rcv_input=args.rcv)
        
        print(f"\n=== Evidence Papers for {args.rcv} ===")
        print(json.dumps(results, ensure_ascii=False, indent=2))
        
    except Exception as e:
        print(f"\n[Fatal Error] {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()