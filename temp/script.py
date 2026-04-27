# 파일명: script.py
import os
from pathlib import Path

def setup_acmg_rag_structure(base_path: Path):
    if not isinstance(base_path, Path):
        raise TypeError("base_path는 pathlib.Path 객체여야 합니다.")

    dirs = [
        "core", "components", "tasks", "pipelines", "data/acmg_raw"
    ]

    files = {
        "core/base_vector_db.py": "# 벡터 DB 인터페이스 및 교차 검색 로직\n",
        "components/schemas.py": "# ACMGRule 및 DocumentChunk 스키마\n",
        "tasks/acmg_loader.py": """# ClinGen 웹 데이터를 파싱하여 DB에 적재하는 태스크
import json
from components.schemas import ACMGRule

class ACMGDataLoader:
    def load_from_json(self, file_path: str):
        with open(file_path, 'r') as f:
            data = json.load(f)
        # 기본값 없이 예외처리로 구성
        return [ACMGRule(**item) for item in data]
""",
        "pipelines/acmg_pipeline.py": "# 문헌 검색 + ACMG 규칙 교차 검색 파이프라인\n"
    }

    try:
        for d in dirs:
            (base_path / d).mkdir(parents=True, exist_ok=True)
        
        for f_path, content in files.items():
            with open(base_path / f_path, "w", encoding="utf-8") as f:
                f.write(content)
        
        print(f"✅ ACMG RAG 프로젝트 구조 생성 완료: {base_path}")
    except Exception as e:
        print(f"❌ 오류 발생: {e}")

if __name__ == "__main__":
    setup_acmg_rag_structure(Path("acmg_fetal_rag"))