from src.ncbi_search.clinvar_fetch import ClinVarFetchClient

c = ClinVarFetchClient(email="your_email@example.com")
details = c.debug_fetch_one(
    "https://www.ncbi.nlm.nih.gov/clinvar/variation/12395/",
    dump_xml_path="debug_out/one_12395.xml",      # 선택
    print_xml_head_lines=60,                      # 선택
)

print("parsed:", len(details))
if details:
    print(details[0].to_dict())
