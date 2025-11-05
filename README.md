# Fetch taxonomic annotations of Pfam and Rfam families

To get NCBI account + API key: https://support.nlm.nih.gov/kbArticle/?pn=KA-05317. This is useful to avoid NCBI request limits. Note that at CNRS we can get a NCBI account through Janus.

## Usage

For Rfam families:

```bash
python rfam_taxonomy.py RF00162 --email [YOUR_EMAIL] --api-key [YOUR_API_KEY] --out RF00162.tax.txt
```

For Pfam families:

```bash
python -u pfam_taxonomy.py PF00027 --page-size 150
```

(Note: `-u` means output will be unbuffered, so you can see progress in real time.)