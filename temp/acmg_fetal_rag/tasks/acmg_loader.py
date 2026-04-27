# ClinGen 웹 데이터를 파싱하여 DB에 적재하는 태스크
import json
from components.schemas import ACMGRule

class ACMGDataLoader:
    def load_from_json(self, file_path: str):
        with open(file_path, 'r') as f:
            data = json.load(f)
        # 기본값 없이 예외처리로 구성
        return [ACMGRule(**item) for item in data]
