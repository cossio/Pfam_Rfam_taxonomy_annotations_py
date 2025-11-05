#!/usr/bin/env python3
"""
Rfam → per-sequence taxonomy table (names or taxids)

Example:
  python rfam_taxa_table.py RF00162 \
    --email you@org.edu --api-key YOUR_KEY \
    --output-value name --lowercase --out RF00162.names.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Optional
from xml.etree import ElementTree as ET

import requests

RFAM_REGIONS_URL = "https://rfam.org/family/{fam}/regions"
EFETCH_TAXONOMY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
DESIRED_RANKS_DEFAULT = ("superkingdom","phylum","class","order","family","genus","species")

# ---------- utils ----------

def _chunked(it: Iterable[int], n: int) -> Iterable[List[int]]:
    buf: List[int] = []
    for x in it:
        buf.append(x)
        if len(buf) == n:
            yield buf
            buf = []
    if buf:
        yield buf

def _most_common(xs: List[int]) -> int:
    return Counter(xs).most_common(1)[0][0]

# ---- Handle top level ---

# Canonical domain anchors (NCBI TaxIDs)
DOMAIN_TAXA = {
    2: "Bacteria",
    2157: "Archaea",
    2759: "Eukaryota",
    10239: "Viruses",
    12884: "Viroids",
    12908: "Unclassified",  # unclassified sequences
}
DOMAIN_ORDER = [2, 2157, 2759, 10239, 12884, 12908]

def canonical_domain_from_lineage(lineage: list[dict]) -> dict:
    """
    Return a dict {rank:'superkingdom', name:str, taxid:int|None}
    based on hallmark TaxIDs present anywhere in the lineage.
    Falls back to any lineage node labeled superkingdom if anchors are absent.
    """
    tids = {n.get("taxid") for n in lineage if n.get("taxid") is not None}
    for tid in DOMAIN_ORDER:
        if tid in tids:
            return {"rank": "superkingdom", "name": DOMAIN_TAXA[tid], "taxid": tid}
    # fallback to whatever NCBI provided as superkingdom, if present
    for n in lineage:
        if n.get("rank") == "superkingdom":
            return {"rank": "superkingdom", "name": n.get("name", ""), "taxid": n.get("taxid")}
    return {"rank": "superkingdom", "name": "", "taxid": None}

# ---------- Rfam ----------

def fetch_rfam_regions(family: str, session: Optional[requests.Session] = None) -> List[dict]:
    s = session or requests.Session()
    r = s.get(RFAM_REGIONS_URL.format(fam=family), headers={"Accept":"text/plain"}, timeout=180)
    if r.status_code == 403:
        raise RuntimeError("Rfam regions too large (403). Try a smaller family.")
    r.raise_for_status()
    rows: List[dict] = []
    for ln in r.text.splitlines():
        if not ln or ln.startswith("#"):
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) < 3:
            continue
        try:
            rows.append({"accession": cols[0].strip(),
                         "species": cols[-2].strip(),
                         "taxid": int(cols[-1])})
        except Exception:
            continue
    return rows

# ---------- NCBI taxonomy ----------

def _parse_taxonomy_xml(xml_text: str) -> Dict[int, List[dict]]:
    out: Dict[int, List[dict]] = {}
    root = ET.fromstring(xml_text)
    for taxon in root.findall(".//Taxon"):
        tid_txt = taxon.findtext("TaxId")
        if not tid_txt:
            continue
        try:
            tid = int(tid_txt)
        except ValueError:
            continue
        name = (taxon.findtext("ScientificName") or "").strip()
        rank = (taxon.findtext("Rank") or "").strip()
        lineage: List[dict] = []
        lin_ex = taxon.find("LineageEx")
        if lin_ex is not None:
            for t in lin_ex.findall("Taxon"):
                lineage.append({
                    "taxid": int(t.findtext("TaxId") or "0"),
                    "rank": (t.findtext("Rank") or "").strip(),
                    "name": (t.findtext("ScientificName") or "").strip(),
                })
        lineage.append({"taxid": tid, "rank": rank, "name": name})
        out[tid] = lineage
    return out

def fetch_ncbi_lineages_resilient(
    taxids: Iterable[int],
    email: Optional[str],
    api_key: Optional[str],
    batch_size: int = 200,
    retries: int = 4,
    session: Optional[requests.Session] = None,
) -> Dict[int, List[dict]]:
    s = session or requests.Session()
    params_base = {"db":"taxonomy","retmode":"xml","tool":"rfam_taxa_table"}
    if email: params_base["email"] = email
    if api_key: params_base["api_key"] = api_key

    def fetch_batch(ids: List[int]) -> Dict[int, List[dict]]:
        if not ids:
            return {}
        data = dict(params_base)
        data["id"] = ",".join(str(t) for t in ids)
        backoff = 0.5
        for _ in range(retries):
            try:
                resp = s.post(EFETCH_TAXONOMY_URL, data=data, timeout=180)
                if resp.status_code == 429:
                    time.sleep(backoff); backoff = min(8.0, backoff*2); continue
                resp.raise_for_status()
                text = resp.text.strip()
                if not text or not text.startswith("<"):
                    raise ET.ParseError("non-XML or empty payload")
                return _parse_taxonomy_xml(text)
            except (requests.RequestException, ET.ParseError):
                time.sleep(backoff); backoff = min(8.0, backoff*2); continue
        if len(ids) == 1:
            return {}
        mid = len(ids)//2
        out = fetch_batch(ids[:mid]); out.update(fetch_batch(ids[mid:])); return out

    out: Dict[int, List[dict]] = {}
    unique = list(dict.fromkeys(int(t) for t in taxids if t))
    for chunk in _chunked(unique, batch_size):
        out.update(fetch_batch(chunk))
        time.sleep(0.12 if api_key else 0.34)
    return out

# ---------- lineage → outputs ----------

def lineage_rank_values(lineage: List[dict], ranks: List[str],
                        mode: str, lowercase: bool) -> List[str]:
    """
    mode: 'name' | 'taxid' | 'both'
    Returns values in the same order as `ranks`.
    Enforces canonical superkingdom/domain using hallmark TaxIDs.
    """
    by_rank = {n["rank"]: {"name": n.get("name"), "taxid": n.get("taxid")}
               for n in lineage if n.get("rank")}
    # enforce canonical domain
    dom = canonical_domain_from_lineage(lineage)
    by_rank["superkingdom"] = {"name": dom["name"], "taxid": dom["taxid"]}

    out: List[str] = []
    for r in ranks:
        entry = by_rank.get(r, {})
        nm = entry.get("name") or ""
        if lowercase and nm:
            nm = nm.lower()
        tid = entry.get("taxid")
        if mode == "name":
            out.append(nm)
        elif mode == "taxid":
            out.append("" if tid is None else str(tid))
        else:  # both
            out.append(nm)
            out.append("" if tid is None else str(tid))
    return out

# ---------- main ----------

def main():
    ap = argparse.ArgumentParser(description="Write per-sequence taxonomy for an Rfam family.")
    ap.add_argument("family", help="Rfam accession or ID, e.g., RF00162")
    ap.add_argument("--email", default=None)
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--batch-size", type=int, default=200)
    ap.add_argument("--ranks", default="superkingdom,phylum,class,order,family,genus,species")
    ap.add_argument("--output-value", choices=["name", "taxid", "both"], default="name", help="write names, taxids, or both")
    ap.add_argument("--lowercase", action="store_true", help="lowercase the names (useful if you want 'bacteria')")
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    ranks = [r.strip() for r in args.ranks.split(",") if r.strip()]

    with requests.Session() as s:
        regions = fetch_rfam_regions(args.family, session=s)
        if not regions:
            print("No regions parsed from Rfam.", file=sys.stderr); sys.exit(2)

        # collapse to one taxid per accession
        by_acc: Dict[str, List[int]] = defaultdict(list)
        for r in regions:
            by_acc[r["accession"]].append(r["taxid"])
        acc_taxid: Dict[str, int] = {acc: _most_common(tids) for acc, tids in by_acc.items()}

        lineages = fetch_ncbi_lineages_resilient(
            acc_taxid.values(), email=args.email, api_key=args.api_key,
            batch_size=args.batch_size, session=s,
        )

        # header
        if args.output_value == "name":
            header = ["accession"] + ranks
        elif args.output_value == "taxid":
            header = ["accession"] + [f"{r}_taxid" for r in ranks]
        else:  # both
            hdr = []
            for r in ranks:
                hdr.extend([f"{r}_name", f"{r}_taxid"])
            header = ["accession"] + hdr

        rows_out: List[List[str]] = []
        dropped = 0

        for acc, tid in acc_taxid.items():
            lin = lineages.get(tid)
            if not lin:
                rows_out.append([acc] + [""]* (len(header)-1))
                continue
            vals = lineage_rank_values(lin, ranks, args.output_value, args.lowercase)
            rows_out.append([acc] + vals)

        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            w.writerows(rows_out)

if __name__ == "__main__":
    main()
