import http.client
from Bio import Entrez
from time import sleep
import urllib.error

Entrez.email = "jorge.fdcd@ipht.fr"
Entrez.tool = "seq-to-assembly-mapper"

CHUNK = 50        # NCBI allows 200 IDs per request


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


def fetch_with_retry(request_fn, max_retries=10, delay=1.0):
    """
    request_fn() must return a fresh handle each time.
    This wrapper retries on:
      - IncompleteRead
      - NotXMLError
      - URLError / HTTPError
      - empty response
    """
    for attempt in range(max_retries):
        try:
            handle = request_fn()
            data = Entrez.read(handle)
            return data
        except (http.client.IncompleteRead, Entrez.Parser.NotXMLError,
                urllib.error.URLError, urllib.error.HTTPError) as e:
            print(f"  Error {type(e).__name__}: retry {attempt+1}/{max_retries}")
            sleep(delay)
            continue
    raise RuntimeError("Maximum retries exceeded")


def seq_to_assembly_mapping(seq_ids):
    print(f"Total seq IDs: {len(seq_ids)}")

    # Step 1: seq -> UID
    seq_to_uid = {}
    batch_index = 1
    for batch in chunks(seq_ids, CHUNK):
        print(f"[Step 1] Processing batch {batch_index} "
              f"({len(batch)} seq IDs)...")

        data = fetch_with_retry(
            lambda: Entrez.elink(
                dbfrom="nuccore",
                db="assembly",
                id=batch,
                linkname="nuccore_assembly",
                retmax=0
            )
        )

        sleep(0.2)

        for record in data:
            seq_id = record["IdList"][0]
            linksets = record.get("LinkSetDb", [])
            if linksets:
                uid = linksets[0]["Link"][0]["Id"]
                seq_to_uid[seq_id] = uid
            else:
                seq_to_uid[seq_id] = None

        mapped = sum(1 for v in seq_to_uid.values() if v is not None)
        print(f"  ✓ Mapped so far: {mapped} / {len(seq_to_uid)}")

        batch_index += 1

    # Collect unique UIDs
    uids = list({u for u in seq_to_uid.values() if u is not None})
    print(f"Unique Assembly UIDs to fetch: {len(uids)}")

    # Step 2: UID -> GCA
    uid_to_gca = {}
    batch_index = 1
    for batch in chunks(uids, CHUNK):
        print(f"[Step 2] Fetching assembly summaries batch {batch_index} "
              f"({len(batch)} UIDs)...")

        data = fetch_with_retry(
            lambda: Entrez.esummary(db="assembly", id=batch, report="full")
        )

        sleep(0.2)

        for doc in data["DocumentSummarySet"]["DocumentSummary"]:
            uid = doc.attributes["uid"]
            uid_to_gca[uid] = doc["AssemblyAccession"]

        print(f"  ✓ Assemblies resolved so far: {len(uid_to_gca)}")

        batch_index += 1

    # Step 3: Combine
    print("[Step 3] Finalizing mapping...")
    out = {seq: uid_to_gca.get(uid) for seq, uid in seq_to_uid.items()}

    resolved = sum(1 for v in out.values() if v is not None)
    print(f"Resolved {resolved} / {len(out)} seq IDs to assemblies.")

    return out



file_path = '/Users/jfdcd/projects/2025/disentangle_taxonomy/seq_ids.txt'
with open(file_path) as f:
    seq_ids = [line.strip() for line in f if line.strip()]

mapping = seq_to_assembly_mapping(seq_ids)

with open("seq_to_assembly.tsv", "w") as f:
    for seq, asm in mapping.items():
        f.write(f"{seq}\t{asm}\n")
