from __future__ import annotations

from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree as ET

from Bio import Entrez

from .models import ClinVarVariantDetails
from .utils import extract_variation_id, polite_sleep


class ClinVarFetchClient:
    def __init__(
        self,
        *,
        email: str,
        api_key: str | None = None,
        tool: str = "ncbi-search",
        max_retries: int = 3,
    ) -> None:
        if not email or "@" not in email:
            raise ValueError("Entrez.email 은 필수이며 유효한 이메일이어야 합니다.")
        Entrez.email = email
        Entrez.tool = tool
        if api_key:
            Entrez.api_key = api_key
        self.max_retries = max_retries

    def fetch_details_from_urls(
        self,
        urls_or_ids: Iterable[str],
        *,
        debug: bool = False,
        batch_size: int = 50,
        dump_xml_dir: str | Path | None = None,
        print_xml_head_lines: int = 0,
    ) -> list[ClinVarVariantDetails]:
        # 1) URL/ID -> variation_id(숫자) 정규화
        ids: list[str] = []
        bad: list[str] = []
        for s in urls_or_ids:
            vid = extract_variation_id(s)
            if vid:
                ids.append(vid)
            else:
                bad.append(str(s))

        # 중복 제거(입력 순서 유지)
        seen: set[str] = set()
        uniq_ids: list[str] = []
        for vid in ids:
            if vid not in seen:
                seen.add(vid)
                uniq_ids.append(vid)

        if debug:
            print("=" * 80)
            print(f"[ClinVarFetch DEBUG] extracted={len(ids)} uniq={len(uniq_ids)} bad={len(bad)}")
            if bad:
                print("[ClinVarFetch DEBUG] bad examples:")
                for x in bad[:10]:
                    print("  -", x)
            print("[ClinVarFetch DEBUG] uniq id examples:", uniq_ids[:10])
            print("=" * 80)

        if not uniq_ids:
            return []

        dump_dir: Path | None = None
        if dump_xml_dir is not None:
            dump_dir = Path(dump_xml_dir)
            dump_dir.mkdir(parents=True, exist_ok=True)

        # ✅ FIXED: variation_id(숫자) -> VCV accession으로 매핑
        vcv_ids = self._map_variation_ids_to_vcv(uniq_ids, debug=debug)

        if debug:
            print(f"[ClinVarFetch DEBUG] mapped VCV ids={len(vcv_ids)} (from variation ids={len(uniq_ids)})")
            print("[ClinVarFetch DEBUG] VCV examples:", vcv_ids[:10])

        if not vcv_ids:
            return []

        # 2) VCV accession으로 efetch (rettype=vcv) + parse
        out: list[ClinVarVariantDetails] = []
        chunk_index = 0
        for i in range(0, len(vcv_ids), batch_size):
            chunk_index += 1
            chunk = vcv_ids[i : i + batch_size]
            details = self._efetch_vcv_and_parse(
                chunk,
                debug=debug,
                dump_dir=dump_dir,
                chunk_index=chunk_index,
                print_xml_head_lines=print_xml_head_lines,
            )
            out.extend(details)

        if debug:
            print(f"[ClinVarFetch DEBUG] total parsed={len(out)} from vcv_ids={len(vcv_ids)}")

        return out

    def debug_fetch_one(
        self,
        url_or_id: str,
        *,
        dump_xml_path: str | Path | None = None,
        print_xml_head_lines: int = 80,
    ) -> list[ClinVarVariantDetails]:
        vid = extract_variation_id(url_or_id)
        if not vid:
            raise ValueError(f"Cannot extract variation_id from: {url_or_id}")

        dump_dir = None
        dump_filename = None
        if dump_xml_path is not None:
            p = Path(dump_xml_path)
            p.parent.mkdir(parents=True, exist_ok=True)
            dump_dir = p.parent
            dump_filename = p.name

        # ✅ FIXED: 1개도 mapping -> vcv -> efetch(vcv)로 확인
        vcv_ids = self._map_variation_ids_to_vcv([vid], debug=True)
        if not vcv_ids:
            print("[ClinVarFetch DEBUG] mapping to VCV failed (0 results).")
            return []

        return self._efetch_vcv_and_parse(
            vcv_ids,
            debug=True,
            dump_dir=dump_dir,
            chunk_index=1,
            dump_filename=dump_filename,
            print_xml_head_lines=print_xml_head_lines,
        )

    # -------------------------------------------------
    # ✅ FIXED: bytes -> str decoding
    # -------------------------------------------------
    def _decode_xml(self, payload: str | bytes) -> str:
        if isinstance(payload, str):
            return payload
        try:
            return payload.decode("utf-8")
        except UnicodeDecodeError:
            return payload.decode("latin-1", errors="replace")

    # -------------------------------------------------
    # ✅ FIXED: variation_id -> VCV accession 매핑 (esummary 이용)
    # -------------------------------------------------
    def _map_variation_ids_to_vcv(self, variation_ids: list[str], *, debug: bool = False) -> list[str]:
        """
        variation_id(숫자)로 esummary(db=clinvar)를 호출해서
        DocumentSummary에서 VCV accession(예: VCV000012395.19 or VCV000012395)을 뽑는다.
        """
        # 여러개면 한 번에 가능(길면 batch로 잘라도 됨)
        ids_str = ",".join(variation_ids)
        last_err: Exception | None = None

        for attempt in range(1, self.max_retries + 1):
            try:
                polite_sleep()
                with Entrez.esummary(db="clinvar", id=ids_str) as h:
                    res = Entrez.read(h)

                # 안전 파싱
                docs = []
                if isinstance(res, dict):
                    dss = res.get("DocumentSummarySet")
                    if isinstance(dss, dict):
                        docs = dss.get("DocumentSummary") or []
                if not isinstance(docs, list):
                    docs = []

                vcv_ids: list[str] = []

                for d in docs:
                    if not isinstance(d, dict):
                        continue

                    # ✅ FIXED: 여러 형태 대비해서 가능한 필드들을 순서대로 시도
                    # (환경에 따라 키가 조금 다를 수 있음)
                    cand = (
                        d.get("accession")  # 흔히 VCV0000... 형태
                        or d.get("accession_version")  # 있으면 버전 포함
                        or d.get("Accession")
                        or d.get("AccessionVersion")
                    )

                    # accession이 없는 경우: variation_set 쪽에 variation_id만 있는 케이스도 있어서 skip
                    if not cand:
                        continue

                    s = str(cand).strip()

                    # VCV prefix가 아니면 제외(안전)
                    if not s.startswith("VCV"):
                        continue

                    # ClinVar efetch는 버전 유무 둘 다 받는 편이라 그대로 추가
                    vcv_ids.append(s)

                # 중복 제거(순서 유지)
                seen: set[str] = set()
                uniq: list[str] = []
                for x in vcv_ids:
                    if x not in seen:
                        seen.add(x)
                        uniq.append(x)

                if debug:
                    print(f"[ClinVarFetch DEBUG] esummary->VCV mapped={len(uniq)}")
                    if len(uniq) == 0:
                        # 진단용: 키 구조 확인
                        try:
                            keys = list(res.keys()) if hasattr(res, "keys") else None
                            print(f"[ClinVarFetch DEBUG] esummary top keys={keys}")
                        except Exception:
                            pass

                return uniq

            except Exception as e:
                last_err = e
                if debug:
                    print(f"[ClinVarFetch DEBUG] esummary mapping failed attempt={attempt}: {e}")
                polite_sleep(base=0.8 * attempt)

        if debug and last_err:
            print(f"[ClinVarFetch DEBUG] esummary mapping giving up: {last_err}")
        return []

    # -------------------------------------------------
    # ✅ FIXED: VCV efetch + parse
    # -------------------------------------------------
    def _efetch_vcv_and_parse(
        self,
        vcv_ids: list[str],
        *,
        debug: bool = False,
        dump_dir: Path | None = None,
        chunk_index: int = 1,
        dump_filename: str | None = None,
        print_xml_head_lines: int = 0,
    ) -> list[ClinVarVariantDetails]:
        ids_str = ",".join(vcv_ids)
        last_err: Exception | None = None

        for attempt in range(1, self.max_retries + 1):
            try:
                polite_sleep()
                with Entrez.efetch(db="clinvar", id=ids_str, rettype="vcv", retmode="xml") as h:
                    raw = h.read()

                xml_text = self._decode_xml(raw)

                if debug:
                    print("-" * 80)
                    print(
                        f"[ClinVarFetch DEBUG] efetch(vcv) chunk={chunk_index} attempt={attempt} "
                        f"ids={vcv_ids[:3]}{'...' if len(vcv_ids)>3 else ''} bytes={len(raw)}"
                    )

                if dump_dir is not None:
                    name = dump_filename or f"clinvar_vcv_chunk_{chunk_index}.xml"
                    out_path = dump_dir / name
                    out_path.write_text(xml_text, encoding="utf-8")
                    if debug:
                        print(f"[ClinVarFetch DEBUG] dumped xml: {out_path}")

                if debug and print_xml_head_lines > 0:
                    head = "\n".join(xml_text.splitlines()[:print_xml_head_lines])
                    print("[ClinVarFetch DEBUG] XML head:")
                    print(head)

                root = ET.fromstring(xml_text)
                if debug:
                    self._debug_xml_structure(root)

                parsed = self._parse_from_vcv_xml(root, debug=debug)

                if debug:
                    print(f"[ClinVarFetch DEBUG] parsed_count={len(parsed)} (chunk={chunk_index})")

                return parsed

            except Exception as e:
                last_err = e
                if debug:
                    print(f"[ClinVarFetch DEBUG] efetch(vcv)/parse failed attempt={attempt}: {e}")
                polite_sleep(base=0.8 * attempt)

        if debug and last_err:
            print(f"[ClinVarFetch DEBUG] efetch(vcv) giving up chunk={chunk_index}: {last_err}")

        return []

    def _debug_xml_structure(self, root: ET.Element) -> None:
        print(f"[ClinVarFetch DEBUG] root_tag={root.tag} root_attrib_keys={list(root.attrib.keys())}")

        xpaths = [
            ".//VariationArchive",
            ".//ClinVarResult",
            ".//TraitSet",
            ".//Interpretation",
            ".//MeasureSet",
        ]
        for xp in xpaths:
            cnt = len(root.findall(xp))
            print(f"[ClinVarFetch DEBUG] xpath_count {xp} = {cnt}")

        va_cnt = len(root.findall(".//VariationArchive"))
        if va_cnt == 0:
            sample_tags: list[str] = []
            for el in root.iter():
                sample_tags.append(el.tag)
                if len(sample_tags) >= 15:
                    break
            print("[ClinVarFetch DEBUG] VariationArchive=0. sample first tags:")
            for t in sample_tags:
                print("  -", t)

    def _parse_from_vcv_xml(self, root: ET.Element, *, debug: bool = False) -> list[ClinVarVariantDetails]:
        out: list[ClinVarVariantDetails] = []
        archives = root.findall(".//VariationArchive")

        if debug and not archives:
            print("[ClinVarFetch DEBUG] (vcv) No VariationArchive nodes found.")
            # Info 노드가 있으면 원인 출력
            info = root.findtext(".//Info")
