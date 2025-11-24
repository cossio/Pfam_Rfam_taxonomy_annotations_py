#!/usr/bin/env python3
import os, re, gzip, shutil, subprocess, numpy
import urllib.request
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

# Reference:
# https://www.nature.com/articles/s41467-024-54812-y
# See their supplementary data and Zenodo for the files used here.

# Configuration
ZENODO_TABLE = Path("data/garnet/genomic_coordinates_228_RNA.csv")       
OGT_TABLE = Path("data/garnet/Supplementary_Table_2_OGT_Comparison.csv") 
OUTPUT = Path("out/RF00162_with_OGT_garnet.tsv")
DOWNLOAD_DIR = Path("data/garnet/genomes")
DOWNLOAD_DIR.mkdir(exist_ok=True)

RFAM_TARGET = "RF00162"


def gtdb_to_ncbi(gtdb_acc):
    """
    Convert GTDB accession like RS_GCF_000005825.2 → GCF_000005825.2
    GTDB format:  {GB|RS}_{GCA|GCF}_XXXXXXXXX.X
    """
    parts = gtdb_acc.split("_", 1)
    if len(parts) != 2:
        raise ValueError(f"Bad GTDB accession: {gtdb_acc}")
    return parts[1]   # drop GB_ / RS_


def download_genome(ncbi_acc):
    """
    Download genome using: datasets download genome accession <ACC>
    Stores genome into DOWNLOAD_DIR/<ACC>/
    """
    outzip = DOWNLOAD_DIR / f"{ncbi_acc}.zip"
    outfolder = DOWNLOAD_DIR / ncbi_acc
    
    if outfolder.exists():
        return outfolder
    
    print(f"[DOWNLOAD] {ncbi_acc}")
    subprocess.run([
        "datasets", "download", "genome", "accession", ncbi_acc,
        "--include", "genome", "--filename", str(outzip)
    ], check=True)
    
    subprocess.run(["unzip", "-o", str(outzip), "-d", str(outfolder)], check=True)
    return outfolder


def download_genome_via_ftp(ncbi_acc, outdir=DOWNLOAD_DIR):
    """
    Download genomic FASTA for an NCBI assembly accession (GCA_xxx or GCF_xxx)
    using the NCBI genomes-all FTP.

    ncbi_acc: e.g. "GCA_024511045.1"
    outdir: directory to save downloaded files

    Returns: Path to the downloaded genomic FASTA (.gz)
    """
    outdir = Path(outdir)
    outdir.mkdir(exist_ok=True)

    # ---------------------------------------------------------------------
    # 1. Compute the FTP path prefix
    # ---------------------------------------------------------------------
    # Split into three-digit chunks: required by NCBI FTP structure
    prefix = ncbi_acc[:3]  # GCA or GCF
    digits = ncbi_acc.split("_")[1].split(".")[0]  # e.g. "024511045"
    p1, p2, p3 = digits[:3], digits[3:6], digits[6:9]

    ftp_base = (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/"
        f"{prefix}/{p1}/{p2}/{p3}/"
    )

    # ---------------------------------------------------------------------
    # 2. List the directory for the assembly (must find *the ASM folder*)
    # ---------------------------------------------------------------------
    # Example directory: GCA_024511045.1_ASM2451104v1/
    url = ftp_base
    with urllib.request.urlopen(url) as r:
        html = r.read().decode()

    # Find subdirectories beginning with "<acc>_"
    asm_dirs = re.findall(rf'href="({ncbi_acc}_.*?)\/"', html)
    if not asm_dirs:
        raise RuntimeError(f"Cannot locate ASM subdirectory for {ncbi_acc} in {url}")

    asm_dir = asm_dirs[0]  # only one is expected

    asm_url = ftp_base + asm_dir + "/"

    # ---------------------------------------------------------------------
    # 3. List files inside that ASM subdirectory, find genomic FASTA
    # ---------------------------------------------------------------------
    # with urllib.request.urlopen(asm_url) as r:
    #     html2 = r.read().decode()

    # # Typically: GCA_xxx_ASMxxxxxxv1_genomic.fna.gz
    # fna_files = re.findall(r'href="([^"]*_genomic\.fna\.gz)"', html2)
    # if not fna_files:
    #     raise RuntimeError(f"No genomic FASTA found in {asm_url}")

    #fna_gz_name = fna_files[0]
    fna_gz_name = f'{asm_dir}_genomic.fna.gz'
    download_url = asm_url + fna_gz_name

    # ---------------------------------------------------------------------
    # 4. Download to local directory
    # ---------------------------------------------------------------------
    gz_outpath = outdir / fna_gz_name

    if not gz_outpath.exists():
        print(f"[FTP DOWNLOAD] {download_url}")
        try:
            urllib.request.urlretrieve(download_url, gz_outpath)
        except urllib.error.HTTPError as e:
            print(f"Error downloading {download_url}: {e}. Skipping.")
            return None
    else:
        print(f"[CACHE] Using existing file: {gz_outpath}")

    fna_outpath = gz_outpath.with_suffix("")
    print(f"[EXTRACT] {gz_outpath} → {fna_outpath}")
    with gzip.open(gz_outpath, "rb") as f_in, open(fna_outpath, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    return fna_outpath



def find_genomic_fna(folder):
    """
    Locate genomic FASTA inside NCBI datasets folder.
    Typical path: ncbi_dataset/data/<ACC>/*.fna
    """
    for f in folder.rglob("*.fna"):
        return f
    raise FileNotFoundError(f"No .fna in {folder}")


def extract_sequence(genome_fna, contig, start, stop):
    """
    Extract the RNA sequence using 1-based inclusive coordinates.
    Reverse complement if start > stop.
    Returns DNA sequence (string).
    """
    # Load only matching contig
    for rec in SeqIO.parse(str(genome_fna), "fasta"):
        if rec.id == contig:
            seq = rec.seq
            break
    else:
        raise ValueError(f"Contig {contig} not found in {genome_fna}")

    if start < 1 or stop < 1:
        raise ValueError("Coordinates must be 1-based positive integers.")

    s = min(start, stop)
    e = max(start, stop)

    # Python slicing: seq[s-1:e] retrieves positions s..e inclusive
    subseq = seq[s-1:e]

    # If reversed orientation (start > stop), reverse complement
    if start > stop:
        subseq = subseq.reverse_complement()

    return str(subseq)


####################################################################
# MAIN
####################################################################

print("[LOAD] Zenodo genomic coordinates table…")
df = pd.read_csv(ZENODO_TABLE, sep=",", dtype=str)

# Clean column names in case Zenodo has # prefixes
df.columns = [c.lstrip("#") for c in df.columns]

print(f"[FILTER] Selecting {RFAM_TARGET} rows…")
rfam_df = df[df["Rfam_family"] == RFAM_TARGET].copy()

if rfam_df.empty:
    print(f"No {RFAM_TARGET} rows found. Check the file.")

# Limit for testing
# rfam_df = rfam_df.head(20)
# rfam_df = rfam_df[rfam_df.GTDB_accession == 'GB_GCA_900604495.1']

print(f"Found {len(rfam_df)} {RFAM_TARGET} entries.")
print("[LOAD] OGT table…")
ogt = pd.read_csv(OGT_TABLE, sep=",", dtype=str)
ogt = ogt.rename(columns={"GTDB_Genome_Accession": "gtdb_genome"})

################################################################
# Extract sequences for each SAM riboswitch row
################################################################
dna_seqs = []
rna_seqs = []

gtdb_downloaded = set()
for (idx, (i, row)) in enumerate(rfam_df.iterrows()):
    gtdb = row["GTDB_accession"]
    contig = row["contig"]
    start = int(row["start_pos"])
    stop = int(row["stop_pos"])

    ncbi_acc = gtdb_to_ncbi(gtdb)
    if gtdb not in gtdb_downloaded:
        genome_fna = download_genome_via_ftp(ncbi_acc)
        gtdb_downloaded.add(gtdb)
    else:
        genome_fna_files = list(DOWNLOAD_DIR.rglob(f"*{ncbi_acc}*.fna"))
        if genome_fna_files:
            genome_fna = genome_fna_files[0]
        else:
            genome_fna = None
    
    if genome_fna is None:
        print(f"[SKIP] No genome available for {ncbi_acc}")
        dna_seqs.append("")
        rna_seqs.append("")
    else:
        print("[EXTRACT] Sequence for", gtdb, contig, start, stop, " from ", genome_fna, " (idx =", idx + 1, " out of ", len(rfam_df), ")")

        # Extract DNA sequence
        dna = extract_sequence(genome_fna, contig, start, stop)
        rna = dna.replace("T", "U")

        dna_seqs.append(dna)
        rna_seqs.append(rna)

print("[MERGE] Attaching sequences to SAM table…")
rfam_df["dna_sequence"] = dna_seqs
rfam_df["rna_sequence"] = rna_seqs

print("[DEBUG] rfam_df columns:", rfam_df.columns.tolist())
print("[DEBUG] ogt columns:", ogt.columns.tolist())

rfam_df = rfam_df.rename(columns={"GTDB_accession": "gtdb_genome"})

print("[MERGE] Joining with OGT table…")
rfam_df["gtdb_genome"] = rfam_df["gtdb_genome"].str.strip().astype("category")
ogt["gtdb_genome"] = ogt["gtdb_genome"].str.strip().astype("category")

rfam_df = rfam_df.set_index("gtdb_genome")
ogt = ogt.set_index("gtdb_genome")

result = rfam_df.join(ogt, how="left").reset_index()

print(f"[WRITE] Saving {OUTPUT}")
result.to_csv(OUTPUT, sep="\t", index=False)

print("Done.")


