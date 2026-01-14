from __future__ import annotations

import random
import re
import time
from urllib.parse import urlparse


def normalize_query(q: str) -> str:
    return " ".join((q or "").strip().split())


def polite_sleep(base: float = 0.34, jitter: float = 0.15) -> None:
    time.sleep(max(0.0, base + random.uniform(-jitter, jitter)))


# ✅ NEW: variation URL/ID에서 variation_id 추출
_VAR_RE = re.compile(r"/clinvar/variation/(\d+)(?:/|$)")


def extract_variation_id(s: str) -> str | None:
    """
    input:
      - "12395"
      - "https://www.ncbi.nlm.nih.gov/clinvar/variation/12395/"
    output: "12395" or None
    """
    s = (s or "").strip()
    if not s:
        return None

    # 숫자만이면 그대로
    if s.isdigit():
        return s

    # URL이면 path에서 추출
    try:
        u = urlparse(s)
        m = _VAR_RE.search(u.path)
        if m:
            return m.group(1)
    except Exception:
        pass

    # 혹시 전체 문자열에 포함된 경우도 추출
    m = _VAR_RE.search(s)
    if m:
        return m.group(1)

    return None
