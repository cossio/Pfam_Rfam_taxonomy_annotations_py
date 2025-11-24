import time
import csv
from Bio import Entrez

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------
Entrez.email = "your_email@example.com"  # <--- REPLACE WITH YOUR EMAIL
# Entrez.api_key = "YOUR_NCBI_API_KEY"   # Optional: Un-comment if you have a key
# -----------------------------------------------------------------------------

def chunker(seq, size):
    """Yield successive chunks from seq."""
    for i in range(0, len(seq), size):
        yield seq[i:i + size]

def fetch_assembly_ids(seq_ids):
    """
    Maps Nucleotide Accessions (e.g., CM000733.1) to Assembly Accessions 
    (e.g., GCA_002005145.1) using NCBI Entrez.
    """
    results = {}
    
    # Filter empty lines
    unique_ids = [x.strip() for x in seq_ids if x.strip()]
    total = len(unique_ids)
    print(f"Processing {total} IDs...")

    # Batch size 50 is safe for URL length limits
    BATCH_SIZE = 50
    
    for i, batch in enumerate(chunker(unique_ids, BATCH_SIZE)):
        print(f"  - Processing batch {i + 1}/{int(total/BATCH_SIZE) + 1}...")
        
        try:
            # STEP 1: Search for the Accessions to get UIDs (Integers)
            # We construct a query: "ID1[Accession] OR ID2[Accession]..."
            search_term = " OR ".join([f"{acc}[ACCN]" for acc in batch])
            
            search_handle = Entrez.esearch(db="nuccore", term=search_term, retmax=len(batch))
            search_results = Entrez.read(search_handle)
            uids = search_results['IdList']
            
            if not uids:
                print("    (No matching UIDs found for this batch)")
                continue

            # STEP 2: Map the UIDs back to the original Accessions
            # We need this because esearch doesn't tell us which UID belongs to which Accession
            summary_handle_nuc = Entrez.esummary(db="nucleotide", id=",".join(uids))
            nuc_summaries = Entrez.read(summary_handle_nuc)
            
            uid_to_acc = {}
            for doc in nuc_summaries:
                uid = str(doc['Id']) # Ensure ID is string for matching
                # AccessionVersion is preferred, fallback to Caption
                acc = doc.get('AccessionVersion', doc.get('Caption'))
                uid_to_acc[uid] = acc

            # STEP 3: Link Nucleotide UIDs to Assembly UIDs
            # cmd="neighbor" finds linked entries in other databases
            link_handle = Entrez.elink(
                dbfrom="nuccore", 
                db="assembly", 
                id=",".join(uids),
                linkname="nuccore_assembly",
                cmd="neighbor"
            )
            link_results = Entrez.read(link_handle)
            
            # Gather all Assembly UIDs we need to fetch details for
            assembly_uids_to_fetch = set()
            nuc_to_asm_uids = {}

            for linkset in link_results:
                input_nuc_uid = str(linkset['IdList'][0])
                linked_asm_uids = []
                
                if linkset.get("LinkSetDb"):
                    for link_db in linkset["LinkSetDb"]:
                        if link_db["DbTo"] == "assembly":
                            for link in link_db["Link"]:
                                asm_uid = str(link["Id"])
                                linked_asm_uids.append(asm_uid)
                                assembly_uids_to_fetch.add(asm_uid)
                
                nuc_to_asm_uids[input_nuc_uid] = linked_asm_uids

            # STEP 4: Get the actual Assembly Accession strings (GCA_...)
            uid_to_gca = {}
            if assembly_uids_to_fetch:
                asm_handle = Entrez.esummary(db="assembly", id=",".join(assembly_uids_to_fetch))
                asm_summaries = Entrez.read(asm_handle)
                
                for asm_doc in asm_summaries['DocumentSummarySet']['DocumentSummary']:
                    # Assembly uses 'uid' or 'Id' depending on API version
                    asm_uid = str(asm_doc.get('uid', asm_doc.get('Id')))
                    gca = asm_doc.get('AssemblyAccession', asm_doc.get('Synonym'))
                    uid_to_gca[asm_uid] = gca

            # STEP 5: Build the result
            for nuc_uid, asm_uids in nuc_to_asm_uids.items():
                original_acc = uid_to_acc.get(nuc_uid)
                if original_acc:
                    gca_list = [uid_to_gca.get(u) for u in asm_uids if u in uid_to_gca]
                    results[original_acc] = list(set(gca_list)) # Unique list

            time.sleep(0.5) # Respect NCBI rate limits

        except Exception as e:
            print(f"  ! Error in batch: {e}")

    return results

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Example list provided by you
    file_path = '/Users/jfdcd/projects/2025/disentangle_taxonomy/seq_ids.txt'
    with open(file_path) as f:
        input_ids = [line.strip() for line in f if line.strip()]
        input_ids = input_ids[:200]
    
    # If reading from a file, uncomment below:
    # with open("my_list.txt", "r") as f:
    #    input_ids = f.readlines()

    mapping = fetch_assembly_ids(input_ids)

    output_filename = "seq_to_assembly_results.csv"
    print(f"\nWriting to {output_filename}...")
    
    with open(output_filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Sequence_ID", "Assembly_ID"])
        
        for seq_id in input_ids:
            seq_id = seq_id.strip()
            if not seq_id: continue
            
            # Check if we found it; if not, mark Not Found
            found = mapping.get(seq_id)
            
            # Note: Some might match partial versions if esearch fuzzy matched, 
            # so we check mapping keys carefully.
            if not found:
                # Fallback: check if versionless key exists (e.g. CM000733)
                base_id = seq_id.split('.')[0]
                found = mapping.get(base_id)

            if found:
                writer.writerow([seq_id, ";".join(found)])
            else:
                writer.writerow([seq_id, "Not Found"])

    print("Done.")