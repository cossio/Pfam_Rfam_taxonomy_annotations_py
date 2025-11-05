# pfam_taxonomy.py

import sys
import time
import random
import argparse
from typing import Iterable, Optional, Dict, Any, List

import requests
import pandas as pd
from ete3 import NCBITaxa

_VERSION = "0.5"

# InterPro: list UniProt proteins for a Pfam entry
BASE = "https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/pfam/{pfam}/?page_size={page_size}"

UA = {
    "User-Agent": f"taxpy/{_VERSION} (+contact:you@example.org)",
    "Accept": "application/json",
    "Connection": "keep-alive",
}

# HTTP codes to retry
RETRIABLE = {408, 429, 500, 502, 503, 504}

# Taxonomic ranks to output
RANKS = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]

ncbi = NCBITaxa()  # downloads taxonomy DB on first use (~100 MB)


# ------------ Networking with retries and progress (logs to stderr) ------------

def _log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)

def _get_json_retry(session: requests.Session, url: str, timeout: int) -> Dict[str, Any]:
    """
    GET JSON with infinite retry on transient failures.
    Retries on RETRIABLE status and any RequestException.
    """
    attempt = 0
    start = time.time()
    while True:
        try:
            r = session.get(url, headers=UA, timeout=timeout)
            if r.status_code in RETRIABLE:
                delay = min(60.0, (2.0 ** attempt) * 0.75) + random.uniform(0.0, 0.5)
                _log(f"[interpro] {r.status_code} retry={attempt} sleep={delay:.1f}s elapsed={time.time()-start:.1f}s")
                time.sleep(delay)
                attempt += 1
                continue
            if r.status_code == 204:
                return {"results": [], "next": None}
            if r.status_code >= 400:
                # Fatal (400/404 etc.)
                r.raise_for_status()
            try:
                return r.json()
            except ValueError:
                delay = min(60.0, (2.0 ** attempt) * 0.75) + random.uniform(0.0, 0.5)
                _log(f"[interpro] non-JSON response retry={attempt} sleep={delay:.1f}s")
                time.sleep(delay)
                attempt += 1
        except requests.exceptions.RequestException as e:
            delay = min(60.0, (2.0 ** attempt) * 0.75) + random.uniform(0.0, 0.5)
            _log(f"[interpro] {e.__class__.__name__} retry={attempt} sleep={delay:.1f}s")
            time.sleep(delay)
            attempt += 1

def _paged_get(session: requests.Session, url: str, timeout: int, limit_pages: Optional[int]) -> Iterable[Dict[str, Any]]:
    """
    Iterate InterPro paginated JSON with robust retry.
    Prints page-level progress. Optional page limit for quick sanity checks.
    """
    pages = 0
    total = 0
    t0 = time.time()
    while url:
        data = _get_json_retry(session, url, timeout)
        results = data.get("results", []) or []
        pages += 1
        total += len(results)
        _log(f"[interpro] page={pages} +{len(results)} proteins total={total} elapsed={time.time()-t0:.1f}s")
        for item in results:
            yield item
        url = data.get("next")
        if limit_pages and pages >= limit_pages:
            _log(f"[interpro] reached limit_pages={limit_pages}")
            break


# ---------------------------- Data extraction ----------------------------

def _iter_pfam_proteins_with_taxa(pfam_id: str,
                                  reviewed: Optional[bool],
                                  page_size: int,
                                  timeout: int,
                                  limit_pages: Optional[int]) -> Iterable[Dict[str, Any]]:
    """
    Yield dicts with accession, tax_id, organism for all proteins in a Pfam family.
    reviewed=True -> Swiss‑Prot only; reviewed=False -> TrEMBL only; None -> both.
    """
    url = BASE.format(pfam=pfam_id, page_size=page_size)
    _log(f"[start] pfam={pfam_id} reviewed={reviewed} page_size={page_size} timeout={timeout}s version={_VERSION}")

    with requests.Session() as session:
        for res in _paged_get(session, url, timeout=timeout, limit_pages=limit_pages):
            meta = res.get("metadata", {}) or {}
            srcdb = meta.get("source_database")  # 'reviewed' or 'unreviewed'
            if reviewed is True and srcdb != "reviewed":
                continue
            if reviewed is False and srcdb != "unreviewed":
                continue

            acc = meta.get("accession")
            org = meta.get("source_organism") or meta.get("organism") or {}
            tax_id = org.get("taxId") or org.get("taxid")
            org_name = org.get("scientificName") or org.get("fullName")

            if acc and tax_id:
                try:
                    yield {"accession": acc, "tax_id": int(tax_id), "organism": org_name}
                except (TypeError, ValueError):
                    continue


def _ranks_for_taxids(tax_ids: list[int], ranks: list[str]) -> pd.DataFrame:
    """
    Map each tax_id to requested ranks.
    Special-case 'superkingdom' to handle taxonomy snapshots that mark it as 'domain'
    or leave it unranked. Fallback uses canonical top-level taxids.
    """
    # canonical top-level nodes
    TOP = {
        2: "Bacteria",
        2157: "Archaea",
        2759: "Eukaryota",
        10239: "Viruses",
        # add more if you care: 12884 "Viroids", etc.
    }

    uniq = sorted(set(tax_ids))
    rows = {}

    for tid in uniq:
        try:
            lin = ncbi.get_lineage(tid)                   # list of taxids (root -> tid)
            rmap = ncbi.get_rank(lin)                     # taxid -> rank str
            nmap = ncbi.get_taxid_translator(lin)         # taxid -> name

            # fill all requested ranks in order
            row = {r: "" for r in ranks}
            for t in lin:
                r = rmap.get(t)
                if r in row and not row[r]:
                    row[r] = nmap.get(t, "")

            # robust superkingdom:
            if "superkingdom" in row and not row["superkingdom"]:
                # try 'domain' label if present in this snapshot
                dom = next((nmap.get(t, "") for t in lin if rmap.get(t) == "domain"), "")
                if dom:
                    row["superkingdom"] = dom
                else:
                    # infer from canonical top-level nodes
                    for t in lin:
                        if t in TOP:
                            row["superkingdom"] = TOP[t]
                            break

            rows[tid] = row
        except Exception:
            rows[tid] = {r: "" for r in ranks}

    df = (
        pd.DataFrame.from_dict(rows, orient="index")
        .rename_axis("tax_id")
        .reset_index()
    )
    # ensure column order is exactly: tax_id + ranks
    return df[["tax_id", *ranks]]


def pfam_taxonomy_table(pfam_id: str,
                        reviewed: Optional[bool] = None,
                        ranks: Optional[List[str]] = None,
                        page_size: int = 150,
                        timeout: int = 45,
                        limit_pages: Optional[int] = None) -> pd.DataFrame:
    """
    Return one row per sequence with taxonomy ranks.
    Columns: accession, tax_id, organism, <ranks...>
    """
    if ranks is None:
        ranks = RANKS

    rows = list(_iter_pfam_proteins_with_taxa(
        pfam_id, reviewed=reviewed, page_size=page_size, timeout=timeout, limit_pages=limit_pages
    ))
    if not rows:
        raise ValueError(f"No proteins or TaxIDs found for {pfam_id}.")

    seq_df = (
        pd.DataFrame(rows)
        .drop_duplicates(subset=["accession"])
        .reset_index(drop=True)
    )
    _log(f"[interpro] sequences={len(seq_df)} unique_taxa={seq_df['tax_id'].nunique()}")

    tax_df = _ranks_for_taxids(seq_df["tax_id"].tolist(), ranks=ranks)

    out = seq_df.merge(tax_df, on="tax_id", how="left")
    # Fallback for species
    out["species"] = out["species"].mask(out["species"] == "", out["organism"])
    # Make sure superkingdom is 4th column
    final_cols = ["accession", "tax_id", "organism", *RANKS]
    out = out.reindex(columns=final_cols)

    return out


# ------------------------------- CLI -------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Export per-sequence taxonomy labels for a Pfam family.")
    p.add_argument("pfam", help="Pfam accession, e.g., PF00027")
    p.add_argument("subset", nargs="?", default="", help="'reviewed' or 'unreviewed'")
    p.add_argument("--page-size", type=int, default=150, help="InterPro page_size (smaller -> fewer timeouts)")
    p.add_argument("--timeout", type=int, default=45, help="Per-request timeout seconds")
    p.add_argument("--limit-pages", type=int, default=None, help="For a quick test, only fetch this many pages")
    p.add_argument("--ranks", default=",".join(RANKS), help="Comma-separated ranks to include")
    return p.parse_args()

if __name__ == "__main__":
    args = _parse_args()
    rev = None
    if args.subset.lower() in ("reviewed", "true", "swiss"):
        rev = True
    elif args.subset.lower() in ("unreviewed", "false", "trembl"):
        rev = False
    ranks = [r.strip() for r in args.ranks.split(",") if r.strip()]

    df = pfam_taxonomy_table(
        args.pfam,
        reviewed=rev,
        ranks=ranks,
        page_size=args.page_size,
        timeout=args.timeout,
        limit_pages=args.limit_pages,
    )
    out = f"{args.pfam}_seq_taxonomy.tsv"
    df.to_csv(out, sep="\t", index=False)
    _log(f"[done] wrote {out} rows={len(df)}")
