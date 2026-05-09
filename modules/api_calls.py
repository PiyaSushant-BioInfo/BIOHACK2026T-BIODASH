"""
TBAnalytica API Module
Handles ALL live API calls with data quality scoring, caching,
and multi-source reconciliation. Every query — known or unknown —
makes live API calls. Local DB is a fallback, not a source of truth.
"""

import hashlib
import json
import os
import time
from datetime import datetime, timedelta
from enum import Enum
from pathlib import Path
from typing import Optional

import requests
from dotenv import load_dotenv
from pydantic import BaseModel, Field

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from schema import DataSource, DataQualityScore

load_dotenv(Path(__file__).resolve().parent.parent / ".env")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "")
UNIPROT_BASE = os.getenv("UNIPROT_BASE_URL", "https://rest.uniprot.org")
PDB_BASE = os.getenv("PDB_BASE_URL", "https://data.rcsb.org/rest/v1")
ALPHAFOLD_BASE = os.getenv("ALPHAFOLD_BASE_URL", "https://alphafold.ebi.ac.uk/api")

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BLAST_BASE = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

CACHE_DIR = Path(__file__).resolve().parent.parent / ".cache"
CACHE_DIR.mkdir(exist_ok=True)

TB_REVIEWED_UNIPROT = {
    "rpoB": "P9WGY9",
    "katG": "P9WIE5",
    "inhA": "P9WGR1",
    "gyrA": "P9WG51",
    "gyrB": "P9WG49",
    "pncA": "P9WIB7",
    "embB": "P9WGE7",
}

# ---------------------------------------------------------------------------
# Data Quality Scoring System
# ---------------------------------------------------------------------------


# DataSource and DataQualityScore are imported from schema.py


SOURCE_CREDIBILITY: dict[DataSource, float] = {
    DataSource.WHO_CATALOGUE: 30.0,
    DataSource.NCBI_REFSEQ: 25.0,
    DataSource.UNIPROT_REVIEWED: 25.0,
    DataSource.PDB: 20.0,
    DataSource.NCBI_GENBANK: 15.0,
    DataSource.UNIPROT_UNREVIEWED: 10.0,
    DataSource.ALPHAFOLD: 8.0,
    DataSource.LOCAL_DB: 5.0,
}


def _score_source(source: DataSource) -> float:
    return SOURCE_CREDIBILITY.get(source, 5.0)


def _score_citations(count: int) -> float:
    if count > 50:
        return 20.0
    if count >= 10:
        return 15.0
    if count >= 1:
        return 8.0
    return 0.0


def _score_completeness(
    has_protein_seq: bool = False,
    has_nucleotide_seq: bool = False,
    has_3d_experimental: bool = False,
    has_3d_predicted: bool = False,
    has_resistance_annotations: bool = False,
    has_binding_site: bool = False,
) -> float:
    s = 0.0
    if has_protein_seq:
        s += 10.0
    if has_nucleotide_seq:
        s += 10.0
    if has_3d_experimental:
        s += 10.0
    elif has_3d_predicted:
        s += 5.0
    if has_resistance_annotations:
        s += 5.0
    if has_binding_site:
        s += 5.0
    return min(s, 30.0)


def _score_recency(last_updated: Optional[datetime]) -> float:
    if last_updated is None:
        return 5.0
    age = datetime.now() - last_updated
    if age < timedelta(days=30):
        return 20.0
    if age < timedelta(days=180):
        return 15.0
    if age < timedelta(days=365):
        return 10.0
    if age < timedelta(days=730):
        return 5.0
    return -5.0


def build_quality_score(
    source: DataSource,
    citations: int = 0,
    has_protein_seq: bool = False,
    has_nucleotide_seq: bool = False,
    has_3d_experimental: bool = False,
    has_3d_predicted: bool = False,
    has_resistance_annotations: bool = False,
    has_binding_site: bool = False,
    last_updated: Optional[datetime] = None,
    review_status: str = "unknown",
    extra_details: Optional[dict] = None,
) -> DataQualityScore:
    score = DataQualityScore(
        source=source,
        source_score=_score_source(source),
        citation_score=_score_citations(citations),
        completeness_score=_score_completeness(
            has_protein_seq, has_nucleotide_seq,
            has_3d_experimental, has_3d_predicted,
            has_resistance_annotations, has_binding_site,
        ),
        recency_score=_score_recency(last_updated),
        review_status=review_status,
        details=extra_details or {},
    )
    return score.compute()


# ---------------------------------------------------------------------------
# 7. Cache System
# ---------------------------------------------------------------------------


def _cache_key(namespace: str, identifier: str) -> str:
    raw = f"{namespace}:{identifier}"
    return hashlib.sha256(raw.encode()).hexdigest()[:24]


def cache_with_timestamp(
    key: str,
    data: dict,
    quality_score: DataQualityScore,
    ttl_hours: int = 24,
) -> None:
    entry = {
        "data": data,
        "quality_score": quality_score.model_dump(mode="json"),
        "cached_at": datetime.now().isoformat(),
        "ttl_hours": ttl_hours,
    }
    path = CACHE_DIR / f"{key}.json"
    with open(path, "w") as f:
        json.dump(entry, f, indent=2, default=str)


def _cache_get(key: str) -> Optional[tuple[dict, DataQualityScore]]:
    path = CACHE_DIR / f"{key}.json"
    if not path.exists():
        return None
    with open(path, "r") as f:
        entry = json.load(f)
    cached_at = datetime.fromisoformat(entry["cached_at"])
    ttl = timedelta(hours=entry.get("ttl_hours", 24))
    if datetime.now() - cached_at > ttl:
        path.unlink(missing_ok=True)
        return None
    qs = DataQualityScore.model_validate(entry["quality_score"])
    return entry["data"], qs


def _cache_set(namespace: str, identifier: str, data: dict, qs: DataQualityScore, ttl_hours: int = 24) -> None:
    key = _cache_key(namespace, identifier)
    cache_with_timestamp(key, data, qs, ttl_hours)


def _cache_lookup(namespace: str, identifier: str) -> Optional[tuple[dict, DataQualityScore]]:
    key = _cache_key(namespace, identifier)
    return _cache_get(key)


# ---------------------------------------------------------------------------
# Helper: NCBI params
# ---------------------------------------------------------------------------


def _ncbi_params() -> dict:
    params = {}
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    if NCBI_EMAIL:
        params["email"] = NCBI_EMAIL
    return params


def _rate_limit():
    delay = 0.11 if NCBI_API_KEY else 0.34
    time.sleep(delay)


# ---------------------------------------------------------------------------
# 1. NCBI Entrez API
# ---------------------------------------------------------------------------


def fetch_ncbi_sequence(
    accession_id: str,
    db: str = "protein",
    score_quality: bool = True,
) -> tuple[dict, DataQualityScore]:
    """Fetch protein or nucleotide sequence from NCBI with quality scoring."""

    cached = _cache_lookup("ncbi_seq", f"{db}:{accession_id}")
    if cached:
        return cached

    params = {
        **_ncbi_params(),
        "db": db,
        "id": accession_id,
        "rettype": "fasta",
        "retmode": "text",
    }
    resp = requests.get(f"{NCBI_BASE}/efetch.fcgi", params=params, timeout=30)
    _rate_limit()

    if resp.status_code != 200:
        qs = build_quality_score(DataSource.LOCAL_DB, review_status="fetch_failed")
        return {}, qs

    lines = resp.text.strip().split("\n")
    header = lines[0] if lines else ""
    sequence = "".join(lines[1:]) if len(lines) > 1 else ""

    is_refseq = any(accession_id.startswith(p) for p in ["NP_", "NM_", "NC_", "XP_", "XM_"])
    source = DataSource.NCBI_REFSEQ if is_refseq else DataSource.NCBI_GENBANK

    summary = _fetch_ncbi_summary(accession_id, db)
    update_date = _parse_ncbi_date(summary.get("updatedate", ""))

    citations = 0
    if score_quality:
        citations = fetch_pubmed_citations(accession_id)

    qs = build_quality_score(
        source=source,
        citations=citations,
        has_protein_seq=(db == "protein" and bool(sequence)),
        has_nucleotide_seq=(db == "nucleotide" and bool(sequence)),
        last_updated=update_date,
        review_status="RefSeq curated" if is_refseq else "GenBank automated",
        extra_details={"accession": accession_id, "length": len(sequence)},
    )

    data = {
        "accession": accession_id,
        "header": header,
        "sequence": sequence,
        "length": len(sequence),
        "is_refseq": is_refseq,
        "db": db,
    }

    _cache_set("ncbi_seq", f"{db}:{accession_id}", data, qs, ttl_hours=24)
    return data, qs


def _fetch_ncbi_summary(accession_id: str, db: str = "protein") -> dict:
    params = {
        **_ncbi_params(),
        "db": db,
        "id": accession_id,
        "retmode": "json",
    }
    try:
        resp = requests.get(f"{NCBI_BASE}/esummary.fcgi", params=params, timeout=15)
        _rate_limit()
        if resp.status_code == 200:
            data = resp.json()
            result = data.get("result", {})
            uids = result.get("uids", [])
            if uids:
                return result.get(uids[0], {})
    except Exception:
        pass
    return {}


def _parse_ncbi_date(date_str: str) -> Optional[datetime]:
    if not date_str:
        return None
    for fmt in ("%Y/%m/%d", "%Y-%m-%d", "%Y/%m/%d %H:%M"):
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    return None


def fetch_ncbi_for_known_variant(
    variant_name: str,
    gene_name: str,
) -> tuple[dict, DataQualityScore]:
    """Search NCBI for the latest sequence data for a KNOWN variant gene."""

    cached = _cache_lookup("ncbi_known", f"{gene_name}:{variant_name}")
    if cached:
        return cached

    gene_id = search_ncbi_gene(gene_name)
    if not gene_id:
        qs = build_quality_score(DataSource.LOCAL_DB, review_status="gene_not_found")
        return {}, qs

    link_params = {
        **_ncbi_params(),
        "dbfrom": "gene",
        "db": "protein",
        "id": gene_id,
        "linkname": "gene_protein_refseq",
        "retmode": "json",
    }
    try:
        resp = requests.get(f"{NCBI_BASE}/elink.fcgi", params=link_params, timeout=20)
        _rate_limit()
        protein_ids = []
        if resp.status_code == 200:
            data = resp.json()
            linksets = data.get("linksets", [])
            for ls in linksets:
                for ldb in ls.get("linksetdbs", []):
                    protein_ids.extend(str(l.get("id", l)) if isinstance(l, dict) else str(l) for l in ldb.get("links", []))
    except Exception:
        protein_ids = []

    if protein_ids:
        best_data, best_qs = fetch_ncbi_sequence(protein_ids[0], db="protein")
        if best_data.get("sequence"):
            _cache_set("ncbi_known", f"{gene_name}:{variant_name}", best_data, best_qs, ttl_hours=24)
            return best_data, best_qs

    nuc_seq = fetch_ncbi_gene_sequence(gene_id)
    if nuc_seq:
        data = {"gene_id": gene_id, "sequence": nuc_seq, "length": len(nuc_seq), "db": "nucleotide"}
        qs = build_quality_score(
            source=DataSource.NCBI_GENBANK,
            has_nucleotide_seq=True,
            review_status="nucleotide_only",
        )
        _cache_set("ncbi_known", f"{gene_name}:{variant_name}", data, qs, ttl_hours=24)
        return data, qs

    qs = build_quality_score(DataSource.LOCAL_DB, review_status="no_ncbi_data")
    return {}, qs


def blast_sequence(
    sequence: str,
    database: str = "nr",
    program: str = "blastp",
    min_quality_score: float = 60.0,
    max_results: int = 10,
) -> list[tuple[dict, DataQualityScore]]:
    """Submit a BLAST search, wait for results, score each hit for quality.
    Returns hits sorted by quality score, filtered by min_quality_score."""

    put_params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": sequence,
        "FORMAT_TYPE": "JSON2",
        "HITLIST_SIZE": str(max_results * 2),
        "ENTREZ_QUERY": "Mycobacterium tuberculosis[ORGN]",
    }
    try:
        resp = requests.post(BLAST_BASE, data=put_params, timeout=30)
    except Exception:
        return []

    rid = None
    for line in resp.text.split("\n"):
        if line.strip().startswith("RID ="):
            rid = line.split("=")[1].strip()
            break
    if not rid:
        return []

    results_json = None
    for attempt in range(30):
        time.sleep(10)
        get_params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "JSON2",
        }
        try:
            resp = requests.get(BLAST_BASE, params=get_params, timeout=30)
            if resp.status_code == 200:
                text = resp.text
                if "Status=WAITING" in text:
                    continue
                if "Status=FAILED" in text or "Status=UNKNOWN" in text:
                    return []
                try:
                    results_json = resp.json()
                    break
                except Exception:
                    if '"BlastOutput2"' in text:
                        results_json = json.loads(text)
                        break
        except Exception:
            continue

    if not results_json:
        return []

    scored_hits = []
    try:
        reports = results_json.get("BlastOutput2", [])
        if not reports:
            return []
        search = reports[0].get("report", {}).get("results", {}).get("search", {})
        hits = search.get("hits", [])
    except (KeyError, IndexError):
        return []

    for hit in hits[:max_results * 2]:
        desc = hit.get("description", [{}])[0]
        accession = desc.get("accession", "")
        title = desc.get("title", "")

        hsps = hit.get("hsps", [{}])
        best_hsp = hsps[0] if hsps else {}
        identity = best_hsp.get("identity", 0)
        align_len = best_hsp.get("align_len", 1)
        pct_identity = (identity / align_len * 100) if align_len else 0

        is_refseq = any(accession.startswith(p) for p in ["NP_", "WP_", "XP_"])
        source = DataSource.NCBI_REFSEQ if is_refseq else DataSource.NCBI_GENBANK

        qs = build_quality_score(
            source=source,
            has_protein_seq=(program == "blastp"),
            has_nucleotide_seq=(program == "blastn"),
            review_status="RefSeq" if is_refseq else "GenBank",
            extra_details={
                "accession": accession,
                "percent_identity": round(pct_identity, 2),
                "e_value": best_hsp.get("evalue", 0),
                "bit_score": best_hsp.get("bit_score", 0),
            },
        )

        if qs.raw_score >= min_quality_score:
            hit_data = {
                "accession": accession,
                "title": title,
                "percent_identity": round(pct_identity, 2),
                "e_value": best_hsp.get("evalue", 0),
                "bit_score": best_hsp.get("bit_score", 0),
                "align_length": align_len,
                "query_coverage": best_hsp.get("query_len", 0),
            }
            scored_hits.append((hit_data, qs))

    scored_hits.sort(key=lambda x: x[1].raw_score, reverse=True)
    return scored_hits[:max_results]


def fetch_pubmed_citations(accession_id: str) -> int:
    """Search PubMed for papers citing a sequence accession."""

    cached = _cache_lookup("pubmed_cit", accession_id)
    if cached:
        return cached[0].get("count", 0)

    params = {
        **_ncbi_params(),
        "db": "pubmed",
        "term": f"{accession_id}[All Fields] AND tuberculosis[Title/Abstract]",
        "rettype": "count",
        "retmode": "json",
    }
    try:
        resp = requests.get(f"{NCBI_BASE}/esearch.fcgi", params=params, timeout=15)
        _rate_limit()
        if resp.status_code == 200:
            data = resp.json()
            count = int(data.get("esearchresult", {}).get("count", 0))
            qs = build_quality_score(DataSource.NCBI_REFSEQ)
            _cache_set("pubmed_cit", accession_id, {"count": count}, qs, ttl_hours=168)
            return count
    except Exception:
        pass
    return 0


# ---------------------------------------------------------------------------
# Backward-compatible wrappers (used by other modules)
# ---------------------------------------------------------------------------


def fetch_ncbi_gene_sequence(gene_id: str) -> Optional[str]:
    """Fetch a nucleotide sequence by gene ID. Legacy interface."""
    params = {
        **_ncbi_params(),
        "db": "nucleotide",
        "id": gene_id,
        "rettype": "fasta",
        "retmode": "text",
    }
    try:
        resp = requests.get(f"{NCBI_BASE}/efetch.fcgi", params=params, timeout=30)
        _rate_limit()
        if resp.status_code == 200:
            lines = resp.text.strip().split("\n")
            return "".join(lines[1:]) if len(lines) > 1 else None
    except Exception:
        pass
    return None


def search_ncbi_gene(gene_name: str, organism: str = "Mycobacterium tuberculosis") -> Optional[str]:
    """Search NCBI Gene for a gene ID. Legacy interface."""
    params = {
        **_ncbi_params(),
        "db": "gene",
        "term": f"{gene_name}[Gene Name] AND {organism}[Organism]",
        "retmax": 1,
        "retmode": "json",
    }
    try:
        resp = requests.get(f"{NCBI_BASE}/esearch.fcgi", params=params, timeout=30)
        _rate_limit()
        if resp.status_code == 200:
            data = resp.json()
            id_list = data.get("esearchresult", {}).get("idlist", [])
            return id_list[0] if id_list else None
    except Exception:
        pass
    return None


def fetch_ncbi_protein_sequence(protein_id: str) -> Optional[str]:
    """Legacy wrapper: fetch protein FASTA from NCBI."""
    data, _ = fetch_ncbi_sequence(protein_id, db="protein", score_quality=False)
    return data.get("sequence") or None


def fetch_ncbi_metadata(accession_id: str) -> dict:
    """Fetch rich metadata for an NCBI accession (protein or nucleotide).

    Tries the protein database first, then nucleotide.  Also parses the
    FASTA header line for organism/strain/gene hints that NCBI's esummary
    sometimes omits.

    Returns a dict with normalised keys:
        accession, title, organism, strain, gene, taxid, header, sequence,
        db  (which NCBI database hit)
    All values default to "" if unavailable.
    """

    result: dict = {
        "accession": accession_id,
        "title": "",
        "organism": "",
        "strain": "",
        "gene": "",
        "taxid": "",
        "header": "",
        "sequence": "",
        "db": "",
    }

    # --- try protein first, then nucleotide ---
    for db in ("protein", "nucleotide"):
        summary = _fetch_ncbi_summary(accession_id, db=db)
        if summary:
            result["title"] = summary.get("title", "")
            result["organism"] = summary.get("organism", "")
            result["taxid"] = str(summary.get("taxid", ""))
            result["db"] = db

            # esummary sometimes stores extra organism detail
            extra_org = summary.get("extra", "")
            if extra_org and not result["organism"]:
                result["organism"] = extra_org

            break  # got a hit

    # --- fetch FASTA for header parsing + sequence ---
    for db in ("protein", "nucleotide"):
        seq_data, _ = fetch_ncbi_sequence(accession_id, db=db, score_quality=False)
        if seq_data.get("sequence"):
            result["header"] = seq_data.get("header", "")
            result["sequence"] = seq_data.get("sequence", "")
            if not result["db"]:
                result["db"] = db
            break

    # --- parse FASTA header for strain / gene hints ---
    header = result["header"]
    if header:
        # Typical FASTA: >acc description [Organism strain]
        # e.g. >VCU48907.1 DNA-directed RNA polymerase subunit beta
        #       [Mycobacterium tuberculosis H37Rv]
        import re

        bracket = re.search(r"\[([^\]]+)\]", header)
        if bracket:
            org_str = bracket.group(1)
            if not result["organism"]:
                result["organism"] = org_str

            # Extract strain from organism string
            # e.g. "Mycobacterium tuberculosis H37Rv" → strain = "H37Rv"
            strain_match = re.search(
                r"Mycobacterium\s+tuberculosis\s+(.+)", org_str, re.IGNORECASE
            )
            if strain_match:
                result["strain"] = strain_match.group(1).strip()

        # Gene name from title portion (before bracket)
        title_part = header.split("[")[0] if "[" in header else header
        # Common TB gene product keywords → gene name mapping
        _PRODUCT_TO_GENE = {
            "rna polymerase": "rpoB",
            "rpob": "rpoB",
            "catalase-peroxidase": "katG",
            "katg": "katG",
            "inha": "inhA",
            "enoyl": "inhA",
            "gyrase": "gyrA",
            "gyra": "gyrA",
            "embB": "embB",
            "arabinosyltransferase": "embB",
            "pnca": "pncA",
            "pyrazinamidase": "pncA",
            "16s rrna": "rrs",
            "rrs": "rrs",
        }
        title_lower = title_part.lower()
        for keyword, gene in _PRODUCT_TO_GENE.items():
            if keyword in title_lower:
                result["gene"] = gene
                break

    return result


# ---------------------------------------------------------------------------
# 2. UniProt API
# ---------------------------------------------------------------------------


def fetch_uniprot_protein(
    uniprot_id: str,
    reviewed_only: bool = False,
) -> tuple[dict, DataQualityScore]:
    """Fetch protein data from UniProt with quality scoring.
    Checks reviewed (SwissProt) vs unreviewed (TrEMBL) status."""

    cached = _cache_lookup("uniprot", uniprot_id)
    if cached:
        return cached

    url = f"{UNIPROT_BASE}/uniprotkb/{uniprot_id}.json"
    try:
        resp = requests.get(url, timeout=30)
    except Exception:
        qs = build_quality_score(DataSource.LOCAL_DB, review_status="fetch_failed")
        return {}, qs

    if resp.status_code != 200:
        qs = build_quality_score(DataSource.LOCAL_DB, review_status="not_found")
        return {}, qs

    entry = resp.json()
    is_reviewed = entry.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"

    if reviewed_only and not is_reviewed:
        alt_data, alt_qs = _find_reviewed_alternative(entry)
        if alt_data:
            return alt_data, alt_qs

    source = DataSource.UNIPROT_REVIEWED if is_reviewed else DataSource.UNIPROT_UNREVIEWED

    sequence = entry.get("sequence", {}).get("value", "")
    features = entry.get("features", [])
    has_binding = any(f.get("type") in ("Binding site", "Active site") for f in features)
    has_resistance = any(
        "resistance" in (f.get("description", "") or "").lower()
        for f in features
    )

    update_str = entry.get("entryAudit", {}).get("lastSequenceUpdateDate", "")
    update_date = None
    if update_str:
        try:
            update_date = datetime.strptime(update_str, "%Y-%m-%d")
        except ValueError:
            pass

    cross_refs = entry.get("uniProtKBCrossReferences", [])
    pdb_ids = [x.get("id") for x in cross_refs if x.get("database") == "PDB"]
    has_structure = len(pdb_ids) > 0

    citations = len(entry.get("references", []))

    qs = build_quality_score(
        source=source,
        citations=citations,
        has_protein_seq=bool(sequence),
        has_3d_experimental=has_structure,
        has_resistance_annotations=has_resistance,
        has_binding_site=has_binding,
        last_updated=update_date,
        review_status="manually reviewed" if is_reviewed else "auto-annotated",
        extra_details={
            "uniprot_id": uniprot_id,
            "protein_name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
            "organism": entry.get("organism", {}).get("scientificName", ""),
            "pdb_ids": pdb_ids[:5],
            "feature_count": len(features),
        },
    )

    data = {
        "uniprot_id": uniprot_id,
        "entry_type": entry.get("entryType", ""),
        "is_reviewed": is_reviewed,
        "protein_name": qs.details.get("protein_name", ""),
        "organism": qs.details.get("organism", ""),
        "sequence": sequence,
        "length": len(sequence),
        "features": [
            {"type": f.get("type", ""), "description": f.get("description", ""), "location": f.get("location", {})}
            for f in features
        ],
        "pdb_ids": pdb_ids,
        "cross_references": {
            x.get("database"): x.get("id")
            for x in cross_refs[:20]
        },
    }

    ttl = 24 if is_reviewed else 12
    _cache_set("uniprot", uniprot_id, data, qs, ttl_hours=ttl)
    return data, qs


def _find_reviewed_alternative(unreviewed_entry: dict) -> tuple[Optional[dict], Optional[DataQualityScore]]:
    """If an entry is unreviewed, search for a reviewed alternative in the same organism."""
    gene_names = []
    for g in unreviewed_entry.get("genes", []):
        for name_obj in g.get("geneName", []) if isinstance(g.get("geneName"), list) else [g.get("geneName", {})]:
            if isinstance(name_obj, dict):
                gene_names.append(name_obj.get("value", ""))
            elif isinstance(name_obj, str):
                gene_names.append(name_obj)

    organism = unreviewed_entry.get("organism", {}).get("scientificName", "Mycobacterium tuberculosis")

    for gene_name in gene_names:
        if not gene_name:
            continue
        search_url = (
            f"{UNIPROT_BASE}/uniprotkb/search"
            f"?query=gene:{gene_name}+AND+organism_name:{organism}+AND+reviewed:true"
            f"&format=json&size=1"
        )
        try:
            resp = requests.get(search_url, timeout=15)
            if resp.status_code == 200:
                results = resp.json().get("results", [])
                if results:
                    reviewed_id = results[0].get("primaryAccession", "")
                    if reviewed_id:
                        return fetch_uniprot_protein(reviewed_id, reviewed_only=False)
        except Exception:
            continue

    return None, None


def fetch_uniprot_for_known_variant(
    gene_name: str,
    organism: str = "Mycobacterium tuberculosis",
) -> tuple[dict, DataQualityScore]:
    """Fetch protein data for a known TB resistance gene.
    Always prefers reviewed (SwissProt) entries using the curated ID map."""

    known_id = TB_REVIEWED_UNIPROT.get(gene_name)
    if known_id:
        return fetch_uniprot_protein(known_id, reviewed_only=False)

    search_url = (
        f"{UNIPROT_BASE}/uniprotkb/search"
        f"?query=gene:{gene_name}+AND+organism_name:{organism}+AND+reviewed:true"
        f"&format=json&size=1"
    )
    try:
        resp = requests.get(search_url, timeout=15)
        if resp.status_code == 200:
            results = resp.json().get("results", [])
            if results:
                uid = results[0].get("primaryAccession", "")
                if uid:
                    return fetch_uniprot_protein(uid, reviewed_only=False)
    except Exception:
        pass

    search_url_unrev = (
        f"{UNIPROT_BASE}/uniprotkb/search"
        f"?query=gene:{gene_name}+AND+organism_name:{organism}"
        f"&format=json&size=1"
    )
    try:
        resp = requests.get(search_url_unrev, timeout=15)
        if resp.status_code == 200:
            results = resp.json().get("results", [])
            if results:
                uid = results[0].get("primaryAccession", "")
                if uid:
                    return fetch_uniprot_protein(uid, reviewed_only=False)
    except Exception:
        pass

    qs = build_quality_score(DataSource.LOCAL_DB, review_status="no_uniprot_data")
    return {}, qs


# Backward-compatible wrappers

def fetch_uniprot_entry(uniprot_id: str) -> Optional[dict]:
    """Legacy wrapper."""
    data, _ = fetch_uniprot_protein(uniprot_id)
    return data if data else None


def fetch_uniprot_sequence(uniprot_id: str) -> Optional[str]:
    """Legacy wrapper: return just the sequence string."""
    data, _ = fetch_uniprot_protein(uniprot_id)
    return data.get("sequence") or None


def fetch_uniprot_features(uniprot_id: str) -> list[dict]:
    """Legacy wrapper: return features list."""
    data, _ = fetch_uniprot_protein(uniprot_id)
    return data.get("features", [])


# ---------------------------------------------------------------------------
# 3. PDB API
# ---------------------------------------------------------------------------


def fetch_pdb_structure(
    pdb_id: str,
) -> tuple[dict, DataQualityScore]:
    """Fetch experimental 3D structure from PDB with quality scoring.
    Checks resolution and experimental method."""

    cached = _cache_lookup("pdb", pdb_id)
    if cached:
        return cached

    url = f"{PDB_BASE}/core/entry/{pdb_id}"
    try:
        resp = requests.get(url, timeout=30)
    except Exception:
        qs = build_quality_score(DataSource.LOCAL_DB, review_status="fetch_failed")
        return {}, qs

    if resp.status_code != 200:
        qs = build_quality_score(DataSource.LOCAL_DB, review_status="not_found")
        return {}, qs

    entry = resp.json()

    methods = entry.get("exptl", [])
    method_name = methods[0].get("method", "UNKNOWN") if methods else "UNKNOWN"

    resolution = None
    refine = entry.get("refine", [])
    if refine:
        resolution = refine[0].get("ls_d_res_high")
        if resolution is not None:
            resolution = float(resolution)

    resolution_bonus = 0.0
    if resolution is not None:
        if resolution < 2.0:
            resolution_bonus = 10.0
        elif resolution < 2.5:
            resolution_bonus = 5.0
        elif resolution < 3.0:
            resolution_bonus = 2.0
        else:
            resolution_bonus = -3.0

    deposition_date = entry.get("rcsb_accession_info", {}).get("deposit_date")
    update_date = None
    if deposition_date:
        try:
            update_date = datetime.strptime(deposition_date[:10], "%Y-%m-%d")
        except ValueError:
            pass

    entity_url = f"{PDB_BASE}/core/polymer_entity/{pdb_id}/1"
    sequence = ""
    try:
        eresp = requests.get(entity_url, timeout=15)
        if eresp.status_code == 200:
            edata = eresp.json()
            sequence = edata.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
    except Exception:
        pass

    qs = build_quality_score(
        source=DataSource.PDB,
        has_protein_seq=bool(sequence),
        has_3d_experimental=True,
        last_updated=update_date,
        review_status=f"experimental ({method_name})",
        extra_details={
            "pdb_id": pdb_id,
            "method": method_name,
            "resolution": resolution,
            "resolution_bonus": resolution_bonus,
        },
    )
    qs.completeness_score = min(qs.completeness_score + resolution_bonus, 30.0)
    qs.compute()

    data = {
        "pdb_id": pdb_id,
        "method": method_name,
        "resolution": resolution,
        "sequence": sequence,
        "length": len(sequence),
        "title": entry.get("struct", {}).get("title", ""),
        "deposition_date": deposition_date,
    }

    _cache_set("pdb", pdb_id, data, qs, ttl_hours=168)
    return data, qs


# Backward-compatible wrappers

def fetch_pdb_info(pdb_id: str) -> Optional[dict]:
    """Legacy wrapper."""
    data, _ = fetch_pdb_structure(pdb_id)
    return data if data else None


def fetch_pdb_sequence(pdb_id: str, entity_id: int = 1) -> Optional[str]:
    """Legacy wrapper."""
    data, _ = fetch_pdb_structure(pdb_id)
    return data.get("sequence") or None


# ---------------------------------------------------------------------------
# 4. AlphaFold API
# ---------------------------------------------------------------------------


def fetch_alphafold_structure(
    uniprot_id: str,
) -> tuple[dict, DataQualityScore]:
    """Fetch AlphaFold predicted structure. Only recommended when no
    experimental PDB structure exists. Incorporates pLDDT confidence."""

    cached = _cache_lookup("alphafold", uniprot_id)
    if cached:
        return cached

    url = f"{ALPHAFOLD_BASE}/prediction/{uniprot_id}"
    try:
        resp = requests.get(url, timeout=30)
    except Exception:
        qs = build_quality_score(DataSource.LOCAL_DB, review_status="fetch_failed")
        return {}, qs

    if resp.status_code != 200:
        qs = build_quality_score(DataSource.LOCAL_DB, review_status="not_found")
        return {}, qs

    raw = resp.json()
    prediction = raw[0] if isinstance(raw, list) else raw

    pdb_url = prediction.get("pdbUrl", "")
    pae_url = prediction.get("paeImageUrl", "")
    entry_id = prediction.get("entryId", "")
    model_date = prediction.get("modelCreatedDate", "")

    global_plddt = prediction.get("globalMetricValue")
    if global_plddt is None:
        confidence_list = prediction.get("confidenceAvgLocalScore")
        if confidence_list is not None:
            global_plddt = confidence_list
        else:
            global_plddt = 70.0

    if isinstance(global_plddt, (int, float)):
        avg_plddt = float(global_plddt)
    else:
        avg_plddt = 70.0

    plddt_bonus = 0.0
    plddt_note = ""
    if avg_plddt > 90:
        plddt_bonus = 5.0
        plddt_note = "Very high confidence"
    elif avg_plddt >= 70:
        plddt_bonus = 2.0
        plddt_note = "Confident"
    elif avg_plddt >= 50:
        plddt_bonus = 0.0
        plddt_note = "Low confidence"
    else:
        plddt_bonus = -5.0
        plddt_note = "Very low confidence — do not use"

    update_date = None
    if model_date:
        try:
            update_date = datetime.strptime(model_date[:10], "%Y-%m-%d")
        except ValueError:
            pass

    qs = build_quality_score(
        source=DataSource.ALPHAFOLD,
        has_protein_seq=True,
        has_3d_predicted=True,
        last_updated=update_date,
        review_status=f"predicted (pLDDT={avg_plddt:.1f})",
        extra_details={
            "entry_id": entry_id,
            "avg_plddt": avg_plddt,
            "plddt_note": plddt_note,
            "pdb_url": pdb_url,
        },
    )
    qs.completeness_score = min(qs.completeness_score + plddt_bonus, 30.0)
    qs.compute()

    data = {
        "entry_id": entry_id,
        "uniprot_id": uniprot_id,
        "pdb_url": pdb_url,
        "pae_image_url": pae_url,
        "avg_plddt": avg_plddt,
        "plddt_note": plddt_note,
        "model_date": model_date,
    }

    _cache_set("alphafold", uniprot_id, data, qs, ttl_hours=168)
    return data, qs


# Backward-compatible wrapper

def fetch_alphafold_prediction(uniprot_id: str) -> Optional[dict]:
    """Legacy wrapper."""
    data, _ = fetch_alphafold_structure(uniprot_id)
    return data if data else None


def fetch_alphafold_pdb_url(uniprot_id: str) -> Optional[str]:
    """Legacy wrapper."""
    data, _ = fetch_alphafold_structure(uniprot_id)
    return data.get("pdb_url") or None


# ---------------------------------------------------------------------------
# 5. WHO TB Mutation Catalogue
# ---------------------------------------------------------------------------

WHO_GRADE_MAP = {
    1: ("Associated with resistance", 95.0),
    2: ("Associated with resistance — interim", 88.0),
    3: ("Uncertain significance", 60.0),
    4: ("Not associated with resistance — interim", 40.0),
    5: ("Not associated with resistance", 30.0),
}


def fetch_who_mutation_data(
    gene: str,
    position: int,
    mutation: str,
) -> tuple[dict, DataQualityScore]:
    """Look up a mutation in the WHO 2022 TB Mutation Catalogue.
    The catalogue is the highest-quality reference source.
    Currently queries the local resistance DB which mirrors the catalogue;
    in production this would hit the WHO genomics API."""

    cached = _cache_lookup("who", f"{gene}:{position}:{mutation}")
    if cached:
        return cached

    from modules.variant_db import load_resistance_db
    db = load_resistance_db()
    gene_data = db["genes"].get(gene, {})

    mut_code = mutation
    if not mut_code and position:
        for code in gene_data.get("mutations", {}):
            if len(code) > 2 and code[1:-1].lstrip("-").isdigit():
                if int(code[1:-1]) == position:
                    mut_code = code
                    break

    mut_info = gene_data.get("mutations", {}).get(mut_code, {})

    if not mut_info:
        silent_info = gene_data.get("silent_precursors", {}).get(mut_code, {})
        if silent_info:
            data = {
                "gene": gene,
                "position": position,
                "mutation": mut_code,
                "who_grade": 3,
                "classification": "Uncertain significance (silent precursor)",
                "confidence": "uncertain",
                "risk_note": silent_info.get("note", ""),
            }
            qs = build_quality_score(
                source=DataSource.WHO_CATALOGUE,
                has_resistance_annotations=True,
                review_status="WHO grade 3",
                extra_details=data,
            )
            qs.raw_score = 60.0
            qs.confidence = "MODERATE"
            qs.use_for_analysis = True
            _cache_set("who", f"{gene}:{position}:{mutation}", data, qs, ttl_hours=168)
            return data, qs

        qs = build_quality_score(DataSource.LOCAL_DB, review_status="not_in_who_catalogue")
        return {}, qs

    resistance = mut_info.get("resistance", "low")
    who_grade = {"high": 1, "moderate": 2, "low": 3}.get(resistance, 3)
    grade_label, grade_score = WHO_GRADE_MAP.get(who_grade, ("Unknown", 50.0))

    drugs = mut_info.get("drugs", [])
    mic = mut_info.get("mic_fold_change", "")

    data = {
        "gene": gene,
        "position": position,
        "mutation": mut_code,
        "who_grade": who_grade,
        "classification": grade_label,
        "resistance_level": resistance,
        "drugs_affected": drugs,
        "mic_fold_change": mic,
        "frequency": mut_info.get("frequency", 0.0),
    }

    qs = build_quality_score(
        source=DataSource.WHO_CATALOGUE,
        has_resistance_annotations=True,
        review_status=f"WHO grade {who_grade}",
        extra_details=data,
    )
    qs.raw_score = max(qs.raw_score, grade_score)
    qs.compute()

    _cache_set("who", f"{gene}:{position}:{mutation}", data, qs, ttl_hours=168)
    return data, qs


# ---------------------------------------------------------------------------
# 6. Data Reconciliation
# ---------------------------------------------------------------------------


def reconcile_data_sources(
    local_data: dict,
    ncbi_data: tuple[dict, DataQualityScore],
    uniprot_data: tuple[dict, DataQualityScore],
    who_data: tuple[dict, DataQualityScore],
) -> tuple[dict, DataQualityScore, list[str]]:
    """Reconcile data from multiple sources, picking the highest-quality
    version and flagging discrepancies."""

    warnings: list[str] = []

    sources = []
    if ncbi_data[0]:
        sources.append(("NCBI", ncbi_data[0], ncbi_data[1]))
    if uniprot_data[0]:
        sources.append(("UniProt", uniprot_data[0], uniprot_data[1]))
    if who_data[0]:
        sources.append(("WHO", who_data[0], who_data[1]))
    if local_data:
        local_qs = build_quality_score(DataSource.LOCAL_DB)
        sources.append(("Local", local_data, local_qs))

    if not sources:
        qs = build_quality_score(DataSource.LOCAL_DB, review_status="no_sources")
        return local_data or {}, qs, ["No live data sources returned results — using local DB only"]

    sources.sort(key=lambda x: x[2].raw_score, reverse=True)
    best_name, best_data, best_qs = sources[0]

    sequences = {}
    for name, data, _ in sources:
        seq = data.get("sequence", "")
        if seq:
            sequences[name] = seq

    if len(sequences) > 1:
        unique_seqs = set(sequences.values())
        if len(unique_seqs) > 1:
            warnings.append(
                f"Sequence DISCREPANCY detected across sources: "
                f"{', '.join(sequences.keys())} — verify manually"
            )
        else:
            pass

    if local_data:
        local_seq = local_data.get("sequence", "")
        if local_seq:
            for name, data, qs in sources:
                if name == "Local":
                    continue
                live_seq = data.get("sequence", "")
                if live_seq and live_seq != local_seq:
                    warnings.append(
                        f"Local database OUTDATED — sequence differs from {name} "
                        f"(live source quality: {qs.raw_score:.0f}/100)"
                    )
                    break

    if who_data[0]:
        who_grade = who_data[0].get("who_grade")
        if who_grade and who_grade >= 3:
            warnings.append(
                f"WHO classification: uncertain significance (grade {who_grade}) — "
                f"treat analysis with caution"
            )

    for name, data, qs in sources:
        if qs.confidence == "REJECT":
            warnings.append(
                f"{name} data REJECTED (quality {qs.raw_score:.0f}/100) — "
                f"not used in analysis"
            )

    return best_data, best_qs, warnings


# ---------------------------------------------------------------------------
# Data Quality Summary Generator
# ---------------------------------------------------------------------------


def generate_quality_summary(
    quality_scores: dict[str, DataQualityScore],
) -> str:
    """Generate a formatted data quality summary for report inclusion."""

    lines = [
        "DATA QUALITY SUMMARY",
        "=" * 50,
    ]

    ok_scores = []
    warn_scores = []

    for label, qs in quality_scores.items():
        line = qs.summary_line(label)
        if qs.confidence in ("HIGH", "MODERATE"):
            ok_scores.append(line)
        else:
            warn_scores.append((label, qs, line))

    for line in ok_scores:
        lines.append(line)

    for label, qs, line in warn_scores:
        lines.append("")
        lines.append(line)
        if qs.confidence == "LOW":
            lines.append(f"   -> Treat {label} analysis with caution")
        elif qs.confidence == "REJECT":
            lines.append(f"   -> {label} data REJECTED — not used")

    if quality_scores:
        avg = sum(qs.raw_score for qs in quality_scores.values()) / len(quality_scores)
        if avg >= 75:
            overall = "HIGH"
        elif avg >= 60:
            overall = "MODERATE"
        else:
            overall = "LOW"
        lines.append("")
        lines.append(f"Overall data confidence: {overall} (avg score: {avg:.0f}/100)")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# High-level: fetch all data for a gene with quality scoring
# ---------------------------------------------------------------------------


def fetch_all_gene_data(
    gene_name: str,
    score_quality: bool = True,
) -> tuple[dict, dict[str, DataQualityScore], list[str]]:
    """Fetch data from ALL sources for a gene, score quality, reconcile.
    Returns (best_data, quality_scores_by_source, warnings)."""

    from modules.variant_db import load_resistance_db

    quality_scores: dict[str, DataQualityScore] = {}
    warnings: list[str] = []

    db = load_resistance_db()
    gene_info = db["genes"].get(gene_name, {})
    local_data = {"gene": gene_name, **gene_info}

    ncbi_data, ncbi_qs = fetch_ncbi_for_known_variant("", gene_name)
    quality_scores[f"{gene_name} (NCBI)"] = ncbi_qs

    uniprot_data, uniprot_qs = fetch_uniprot_for_known_variant(gene_name)
    quality_scores[f"{gene_name} (UniProt)"] = uniprot_qs

    pdb_id = gene_info.get("pdb_id")
    pdb_data: dict = {}
    if pdb_id:
        pdb_data, pdb_qs = fetch_pdb_structure(pdb_id)
        quality_scores[f"{gene_name} (PDB {pdb_id})"] = pdb_qs
    else:
        uniprot_id = gene_info.get("uniprot_id") or TB_REVIEWED_UNIPROT.get(gene_name)
        if uniprot_id:
            af_data, af_qs = fetch_alphafold_structure(uniprot_id)
            quality_scores[f"{gene_name} (AlphaFold)"] = af_qs

    who_data_default = ({}, build_quality_score(DataSource.LOCAL_DB))
    first_mutation = list(gene_info.get("mutations", {}).keys())
    if first_mutation:
        mut_code = first_mutation[0]
        if len(mut_code) > 2 and mut_code[1:-1].lstrip("-").isdigit():
            pos = int(mut_code[1:-1])
            who_d, who_q = fetch_who_mutation_data(gene_name, pos, mut_code)
            quality_scores[f"{gene_name} (WHO)"] = who_q
            who_data_default = (who_d, who_q)

    best, best_qs, recon_warnings = reconcile_data_sources(
        local_data,
        (ncbi_data, ncbi_qs),
        (uniprot_data, uniprot_qs),
        who_data_default,
    )
    warnings.extend(recon_warnings)

    combined = {
        "gene": gene_name,
        "ncbi": ncbi_data,
        "uniprot": uniprot_data,
        "pdb": pdb_data,
        "who": who_data_default[0],
        "best_source": best,
    }

    return combined, quality_scores, warnings


# ---------------------------------------------------------------------------
# Batch helpers
# ---------------------------------------------------------------------------


def batch_fetch_sequences(gene_ids: list[str], db: str = "nucleotide") -> dict[str, str]:
    """Legacy batch fetcher."""
    results = {}
    for gid in gene_ids:
        if db == "nucleotide":
            seq = fetch_ncbi_gene_sequence(gid)
        else:
            seq = fetch_ncbi_protein_sequence(gid)
        if seq:
            results[gid] = seq
    return results


def fetch_ncbi_snp_info(gene_name: str) -> list[str]:
    """Search NCBI SNP database for a gene. Legacy interface."""
    params = {
        **_ncbi_params(),
        "db": "snp",
        "term": f"{gene_name}[Gene Name] AND Mycobacterium tuberculosis[Organism]",
        "retmax": 50,
        "retmode": "json",
    }
    try:
        resp = requests.get(f"{NCBI_BASE}/esearch.fcgi", params=params, timeout=30)
        _rate_limit()
        if resp.status_code == 200:
            data = resp.json()
            return data.get("esearchresult", {}).get("idlist", [])
    except Exception:
        pass
    return []
