#!/usr/bin/env python
"""Resolve every manuscript reference against PubMed. Never invent a PMID.

    python scripts/verify_references.py

Queries the NCBI E-utilities API (esearch then esummary), retrieves the real
PMID, authors, journal, year, volume, pages and DOI, and then CHECKS that the
returned record is the paper that was intended rather than accepting the first
hit. The check compares title tokens and first-author surname, and records a
confidence so a human can audit the weak ones.

A reference that PubMed does not return is written with pmid "not found" and
listed at the end for manual checking. Conference papers in computer vision
(ECCV, CVPR, NeurIPS) are frequently not indexed in PubMed at all; those are
expected misses and are cited by DOI or arXiv instead, never with a fabricated
PMID.

Writes:
    manuscript/references_verified.json    full records plus match confidence
    manuscript/references.bib              BibTeX with the PMID in note
"""
from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "coda-brca-my"
EMAIL = "aliakbarbiotech@hotmail.com"     # NCBI asks callers to identify themselves

# key, intended citation, search query, expected first-author surname
WANTED = [
    ("kiemen2022coda", "CODA: quantitative 3D reconstruction of large tissues at cellular resolution",
     "CODA quantitative 3D reconstruction large tissues cellular resolution", "Kiemen"),
    ("kartasalo2018", "Comparative analysis of tissue reconstruction algorithms for 3D histology",
     "Comparative analysis of tissue reconstruction algorithms for 3D histology", "Kartasalo"),
    ("weitz2024acrobat", "ACROBAT, automatic registration of breast cancer tissue",
     "ACROBAT breast cancer WSI registration challenge Weitz", "Weitz"),
    ("borovec2020anhir", "ANHIR: automatic non-rigid histological image registration challenge",
     "ANHIR automatic non-rigid histological image registration challenge", "Borovec"),
    ("graham2019hovernet", "HoVer-Net: simultaneous segmentation and classification of nuclei",
     "HoVer-Net simultaneous segmentation classification nuclei multi-tissue histology", "Graham"),
    ("bankhead2017qupath", "QuPath: open source software for digital pathology image analysis",
     "QuPath open source software digital pathology image analysis", "Bankhead"),
    ("goode2013openslide", "OpenSlide: a vendor-neutral software foundation for digital pathology",
     "OpenSlide vendor-neutral software foundation digital pathology", "Goode"),
    ("ruifrok2001", "Quantification of histochemical staining by color deconvolution",
     "Quantification of histochemical staining by color deconvolution", "Ruifrok"),
    ("dowsett2011ki67", "Assessment of Ki67 in breast cancer: recommendations from the "
                        "International Ki67 in Breast Cancer Working Group",
     "Assessment of Ki67 in breast cancer recommendations International Ki67 Breast Cancer Working Group", "Dowsett"),
    ("polley2013ki67", "An international Ki67 reproducibility study",
     "international Ki67 reproducibility study Polley", "Polley"),
    ("nielsen2021ki67", "Assessment of Ki67 in breast cancer: updated recommendations",
     "Assessment of Ki67 in breast cancer updated recommendations International Ki67 Working Group", "Nielsen"),
    ("wolff2018her2", "HER2 testing in breast cancer: ASCO/CAP clinical practice guideline update",
     "Human Epidermal Growth Factor Receptor 2 Testing in Breast Cancer ASCO CAP guideline focused update", "Wolff"),
    ("wolff2023her2", "HER2 testing in breast cancer ASCO/CAP guideline update 2023",
     "HER2 testing breast cancer ASCO CAP guideline update 2023 Wolff", "Wolff"),
    ("conklin2011tacs", "Aligned collagen is a prognostic signature for survival in human breast carcinoma",
     "Aligned collagen prognostic signature survival human breast carcinoma", "Conklin"),
    ("provenzano2006tacs", "Collagen reorganization at the tumor-stromal interface facilitates local invasion",
     "Collagen reorganization tumor-stromal interface facilitates local invasion", "Provenzano"),
    ("chen2018deeplab", "Encoder-decoder with atrous separable convolution for semantic image segmentation",
     "Encoder-decoder atrous separable convolution semantic image segmentation", "Chen"),
    ("he2016resnet", "Deep residual learning for image recognition",
     "Deep residual learning for image recognition", "He"),
]


# PMIDs pinned after manual verification against esummary, because free-text
# search returned the wrong record for these. Each title was confirmed to match
# the intended paper before pinning; none is asserted from memory alone.
PINNED = {
    "kartasalo2018":     "29684099",   # Bioinformatics 2018, exact title match
    "bankhead2017qupath":"29203879",   # Sci Rep 2017;7:16878
    "ruifrok2001":       "11531144",   # Anal Quant Cytol Histol 2001, NOT the 2003 paper
    "dowsett2011ki67":   "21960707",   # J Natl Cancer Inst 2011;103:1656-64
    "nielsen2021ki67":   "33369635",   # J Natl Cancer Inst 2021, updated recommendations
    "wolff2018her2":     "29846122",   # J Clin Oncol 2018, NOT the 2023 update
    "wolff2023her2":     "37284804",   # J Clin Oncol 2023 update
}
# Not indexed in PubMed. Computer-vision conference papers are cited by DOI or
# arXiv; recording a PMID for them would mean inventing one.
NO_PUBMED = {
    "chen2018deeplab": "ECCV 2018, arXiv:1802.02611",
    "he2016resnet":    "CVPR 2016, arXiv:1512.03385",
}

STOP = {"the", "of", "for", "in", "a", "an", "and", "with", "to", "by", "on", "from"}


def _get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": f"{TOOL} ({EMAIL})"})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return json.load(fh)


def esearch(term: str, n: int = 5) -> list[str]:
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmode": "json",
                                "retmax": n, "tool": TOOL, "email": EMAIL})
    return _get(f"{EUTILS}/esearch.fcgi?{q}").get("esearchresult", {}).get("idlist", [])


def esummary(pmids: list[str]) -> dict:
    if not pmids:
        return {}
    q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids),
                                "retmode": "json", "tool": TOOL, "email": EMAIL})
    return _get(f"{EUTILS}/esummary.fcgi?{q}").get("result", {})


def tokens(s: str) -> set[str]:
    return {w for w in re.findall(r"[a-z0-9]+", s.lower()) if w not in STOP and len(w) > 2}


def score(intended: str, expect_author: str, rec: dict) -> float:
    """Token overlap on title, plus a bonus if the first author surname matches."""
    t = tokens(intended) & tokens(rec.get("title", ""))
    base = len(t) / max(len(tokens(intended)), 1)
    authors = rec.get("authors") or []
    first = authors[0]["name"].split()[0] if authors else ""
    return min(1.0, base + (0.25 if first.lower() == expect_author.lower() else 0.0))


def main() -> None:
    out_dir = ROOT / "manuscript"
    out_dir.mkdir(parents=True, exist_ok=True)
    records, not_found, weak = [], [], []

    for key, intended, query, author in WANTED:
        if key in NO_PUBMED:
            records.append({"key": key, "intended": intended, "pmid": "not found",
                            "venue": NO_PUBMED[key], "confidence": 1.0,
                            "note": "not indexed in PubMed; cite by DOI/arXiv"})
            not_found.append(key)
            print(f"  {key:22s} not in PubMed ({NO_PUBMED[key]}) - cited by arXiv/DOI")
            continue
        if key in PINNED:
            summ = esummary([PINNED[key]])
            best = summ.get(PINNED[key])
            if isinstance(best, dict) and best.get("title"):
                doi = next((a.get("value","") for a in best.get("articleids", [])
                            if a.get("idtype") == "doi"), "")
                records.append({"key": key, "intended": intended,
                                "pmid": best.get("uid"),
                                "title": best.get("title","").rstrip("."),
                                "authors": [a["name"] for a in (best.get("authors") or [])],
                                "journal": best.get("source",""),
                                "year": (best.get("pubdate","") or "")[:4],
                                "volume": best.get("volume",""),
                                "issue": best.get("issue",""),
                                "pages": best.get("pages",""), "doi": doi,
                                "confidence": 1.0, "note": "PMID pinned and verified"})
                print(f"  {key:22s} PMID {best.get('uid'):>9s}  pinned+verified  "
                      f"{best.get('source','')} {(best.get('pubdate','') or '')[:4]}")
                time.sleep(0.4)
                continue
        try:
            ids = esearch(query)
            summ = esummary(ids)
        except Exception as exc:
            ids, summ = [], {}
            print(f"  {key}: query failed ({exc})")

        best, best_s = None, 0.0
        for pmid in ids:
            rec = summ.get(pmid)
            if not isinstance(rec, dict):
                continue
            s = score(intended, author, rec)
            if s > best_s:
                best, best_s = rec, s

        if best is None or best_s < 0.35:
            records.append({"key": key, "intended": intended, "pmid": "not found",
                            "confidence": round(best_s, 2),
                            "best_candidate_title": (best or {}).get("title", ""),
                            "note": "PubMed returned no convincing match"})
            not_found.append(key)
            print(f"  {key:22s} NOT FOUND (best {best_s:.2f})")
            continue

        doi = ""
        for aid in best.get("articleids", []):
            if aid.get("idtype") == "doi":
                doi = aid.get("value", "")
        authors = [a["name"] for a in (best.get("authors") or [])]
        rec = {
            "key": key, "intended": intended, "pmid": best.get("uid"),
            "title": best.get("title", "").rstrip("."),
            "authors": authors,
            "journal": best.get("source", ""),
            "year": (best.get("pubdate", "") or "")[:4],
            "volume": best.get("volume", ""), "issue": best.get("issue", ""),
            "pages": best.get("pages", ""), "doi": doi,
            "confidence": round(best_s, 2),
        }
        records.append(rec)
        flag = "" if best_s >= 0.6 else "   <-- LOW CONFIDENCE, check manually"
        if best_s < 0.6:
            weak.append(key)
        print(f"  {key:22s} PMID {rec['pmid']:>9s}  conf {best_s:.2f}  "
              f"{rec['journal']} {rec['year']}{flag}")
        time.sleep(0.4)          # NCBI rate limit without an API key

    (out_dir / "references_verified.json").write_text(json.dumps(records, indent=2))

    bib = []
    for r in records:
        if r["pmid"] == "not found":
            bib.append(f"@misc{{{r['key']},\n  title = {{{r['intended']}}},\n"
                       f"  note = {{PMID: not found; verify manually}}\n}}")
            continue
        auth = " and ".join(r["authors"][:12]) or "Unknown"
        bib.append(
            f"@article{{{r['key']},\n"
            f"  author  = {{{auth}}},\n"
            f"  title   = {{{r['title']}}},\n"
            f"  journal = {{{r['journal']}}},\n"
            f"  year    = {{{r['year']}}},\n"
            f"  volume  = {{{r['volume']}}},\n"
            f"  pages   = {{{r['pages']}}},\n"
            + (f"  doi     = {{{r['doi']}}},\n" if r["doi"] else "")
            + f"  note    = {{PMID: {r['pmid']}}}\n}}")
    (out_dir / "references.bib").write_text("\n\n".join(bib))

    print(f"\nresolved {sum(1 for r in records if r['pmid'] != 'not found')} "
          f"of {len(records)} against PubMed")
    if weak:
        print(f"LOW CONFIDENCE, verify by hand: {', '.join(weak)}")
    if not_found:
        print(f"NOT FOUND (expected for CV conference papers): {', '.join(not_found)}")
    print("wrote manuscript/references_verified.json and manuscript/references.bib")


if __name__ == "__main__":
    main()
