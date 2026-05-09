"""
TBAnalytica Treatment Module
Generates treatment recommendations based on variant analysis,
similarity scores, and resistance profiles. Follows WHO 2022 guidelines
adapted for Nepal NTP (National Tuberculosis Programme).

First-line:  Isoniazid (H), Rifampicin (R), Pyrazinamide (Z), Ethambutol (E)
Second-line: Fluoroquinolones (Lfx, Mfx), Bedaquiline (Bdq), Linezolid (Lzd),
             Clofazimine (Cfz), Cycloserine (Cs), Delamanid (Dlm), Pretomanid (Pa)
"""

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from schema import TBVariant, ComparisonResult, DrugSensitivity

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


# ---------------------------------------------------------------------------
# Drug metadata
# ---------------------------------------------------------------------------

DRUG_ABBREVIATIONS: dict[str, str] = {
    "isoniazid": "H",
    "rifampicin": "R",
    "pyrazinamide": "Z",
    "ethambutol": "E",
    "levofloxacin": "Lfx",
    "moxifloxacin": "Mfx",
    "bedaquiline": "Bdq",
    "linezolid": "Lzd",
    "clofazimine": "Cfz",
    "cycloserine": "Cs",
    "delamanid": "Dlm",
    "pretomanid": "Pa",
    "amikacin": "Am",
    "streptomycin": "S",
    "ethionamide": "Eto",
    "para_aminosalicylic_acid": "PAS",
    "imipenem": "Ipm",
    "capreomycin": "Cm",
}

FIRST_LINE = {"isoniazid", "rifampicin", "pyrazinamide", "ethambutol"}

SECOND_LINE = {
    "levofloxacin", "moxifloxacin", "bedaquiline", "linezolid",
    "clofazimine", "cycloserine", "delamanid", "pretomanid",
}

DRUG_CLASSES: dict[str, str] = {
    "isoniazid": "first-line",
    "rifampicin": "first-line",
    "pyrazinamide": "first-line",
    "ethambutol": "first-line",
    "levofloxacin": "fluoroquinolone",
    "moxifloxacin": "fluoroquinolone",
    "bedaquiline": "group-A",
    "linezolid": "group-A",
    "clofazimine": "group-B",
    "cycloserine": "group-B",
    "delamanid": "group-C",
    "pretomanid": "nitroimidazole",
    "amikacin": "injectable",
    "streptomycin": "injectable",
    "ethionamide": "thioamide",
    "para_aminosalicylic_acid": "group-C",
    "imipenem": "carbapenem",
    "capreomycin": "injectable",
}

KNOWN_INTERACTIONS: list[dict] = [
    {
        "drugs": {"bedaquiline", "moxifloxacin"},
        "severity": "moderate",
        "warning": "Both prolong QT interval — ECG monitoring required weekly for first month then biweekly",
    },
    {
        "drugs": {"bedaquiline", "delamanid"},
        "severity": "moderate",
        "warning": "Additive QT prolongation — ECG monitoring required at least weekly; consider avoiding co-administration unless no alternatives",
    },
    {
        "drugs": {"linezolid", "isoniazid"},
        "severity": "low",
        "warning": "Both can cause peripheral neuropathy — monitor closely and supplement pyridoxine",
    },
    {
        "drugs": {"linezolid", "ethionamide"},
        "severity": "low",
        "warning": "Increased risk of serotonin syndrome (rare) — monitor for agitation, tremor, diarrhea",
    },
    {
        "drugs": {"rifampicin", "bedaquiline"},
        "severity": "high",
        "warning": "Rifampicin significantly reduces bedaquiline levels (CYP3A4 induction) — DO NOT co-administer",
    },
    {
        "drugs": {"rifampicin", "delamanid"},
        "severity": "moderate",
        "warning": "Rifampicin may reduce delamanid levels — monitor therapeutic response",
    },
    {
        "drugs": {"rifampicin", "linezolid"},
        "severity": "moderate",
        "warning": "Rifampicin reduces linezolid levels by ~30% — consider dose adjustment or TDM",
    },
    {
        "drugs": {"cycloserine", "ethionamide"},
        "severity": "moderate",
        "warning": "Both can cause CNS toxicity — monitor for seizures, psychosis, and depression",
    },
    {
        "drugs": {"amikacin", "capreomycin"},
        "severity": "high",
        "warning": "Additive nephrotoxicity and ototoxicity — never co-administer",
    },
    {
        "drugs": {"pyrazinamide", "levofloxacin"},
        "severity": "low",
        "warning": "Both can raise uric acid levels — monitor for gout symptoms",
    },
]


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_treatment_protocols() -> dict:
    path = DATA_DIR / "treatment_protocols.json"
    with open(path, "r") as f:
        return json.load(f)


def get_nepal_availability() -> dict:
    protocols = load_treatment_protocols()
    return protocols.get("nepal_drug_availability", {})


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _determine_regimen_key(resistance_class: str, resistant_drugs: list[str]) -> str:
    rc = resistance_class.upper().replace("-", "_")
    if rc == "XDR":
        return "xdr"
    if rc in ("PRE_XDR",):
        return "pre_xdr"
    if rc == "MDR":
        return "mdr"

    resistant_lower = {d.lower() for d in resistant_drugs}
    has_inh = "isoniazid" in resistant_lower
    has_rif = "rifampicin" in resistant_lower

    if has_rif and has_inh:
        return "mdr"
    if has_inh:
        return "isoniazid_resistant"
    if has_rif:
        return "mdr"
    return "susceptible"


def _extract_drugs_from_protocol(protocol: dict) -> tuple[list[str], list[str]]:
    """Extract primary and alternative drug lists from a protocol."""
    regimen = protocol.get("regimen", {})
    primary: list[str] = []
    alternative: list[str] = []

    if "intensive_phase" in regimen:
        intensive = regimen["intensive_phase"].get("drugs", [])
        continuation = regimen.get("continuation_phase", {}).get("drugs", [])
        primary = list(dict.fromkeys(intensive + continuation))
    elif "all_oral_shorter" in regimen:
        primary = regimen["all_oral_shorter"].get("drugs", [])
    elif "primary" in regimen:
        primary = regimen["primary"].get("drugs", [])

    if "alternative" in regimen:
        alternative = regimen["alternative"].get("drugs", [])
    elif "longer_regimen" in regimen:
        lr = regimen["longer_regimen"]
        alternative = lr.get("group_a", []) + lr.get("group_b", [])
    elif "salvage" in regimen:
        alternative = regimen["salvage"].get("drugs", [])

    return primary, alternative


def _get_nepal_avail_map(drugs: list[str], availability: dict) -> dict:
    result = {}
    for drug in drugs:
        key = drug.lower().replace(" ", "_")
        info = availability.get(key, availability.get(drug, {}))
        result[drug] = info if info else {
            "available": False, "source": "Unknown", "cost_category": "N/A",
        }
    return result


def _resistant_drug_set(variant: TBVariant) -> set[str]:
    """Collect all drugs the variant is resistant to."""
    drugs = set()
    all_known_drugs = set(DRUG_ABBREVIATIONS.keys())

    for drug, status in variant.drug_resistance.items():
        if status == DrugSensitivity.RESISTANT:
            dl = drug.lower()
            if dl in all_known_drugs:
                drugs.add(dl)

    from modules.variant_db import load_resistance_db
    db = load_resistance_db()
    for rg in variant.resistance_genes:
        gene_data = db["genes"].get(rg.gene_name, {})
        mutations_db = gene_data.get("mutations", {})
        for m in rg.mutations:
            code = f"{m.reference_amino_acid}{m.position}{m.mutant_amino_acid}"
            for d in mutations_db.get(code, {}).get("drugs", []):
                drugs.add(d.lower())
            if m.drug_affected:
                for d in m.drug_affected.split(","):
                    d_clean = d.strip().lower()
                    if d_clean in all_known_drugs:
                        drugs.add(d_clean)

    return drugs


def _precursor_drug_set(variant: TBVariant) -> set[str]:
    drugs = set()
    for m in variant.all_mutations():
        if m.one_step_away_risk and m.one_step_away_drug:
            for d in m.one_step_away_drug.split(","):
                drugs.add(d.strip().lower())
    return drugs


def _classify_resistance(variant: TBVariant) -> str:
    resistant = _resistant_drug_set(variant)
    has_rif = "rifampicin" in resistant
    has_inh = "isoniazid" in resistant
    has_fq = any(d in resistant for d in ["levofloxacin", "moxifloxacin"])
    has_inj = any(d in resistant for d in ["amikacin", "kanamycin", "capreomycin"])

    if has_rif and has_inh and has_fq and has_inj:
        return "XDR"
    if has_rif and has_inh and has_fq:
        return "Pre-XDR"
    if has_rif and has_inh:
        return "MDR"
    if resistant:
        return "MONO-RESISTANT"
    return "SUSCEPTIBLE"


def _compute_confidence(
    variant: TBVariant,
    nepal_avail: dict,
    recommended_drugs: list[str],
    similarity_score: float | None = None,
) -> float:
    score = 1.0

    all_mutations = variant.all_mutations()
    if not all_mutations:
        score *= 0.5

    well_known = {"katG", "rpoB", "inhA", "embB", "pncA", "gyrA", "gyrB", "rrs"}
    gene_names = {rg.gene_name for rg in variant.resistance_genes}
    if gene_names:
        characterized = sum(1 for g in gene_names if g in well_known)
        score *= 0.5 + 0.5 * (characterized / len(gene_names))

    if recommended_drugs:
        unavailable = sum(
            1 for d in recommended_drugs
            if not nepal_avail.get(d, {}).get("available", False)
        )
        score *= 1 - 0.2 * unavailable / len(recommended_drugs)

    if similarity_score is not None:
        if similarity_score < 60:
            score *= 0.6
        elif similarity_score < 75:
            score *= 0.8

    return round(min(max(score, 0.0), 1.0), 3)


# ---------------------------------------------------------------------------
# 1. get_treatment_for_known_variant
# ---------------------------------------------------------------------------

def get_treatment_for_known_variant(variant_id: str) -> dict:
    """Look up treatment protocol for a known variant.

    Returns regimen, duration, monitoring, and flags drugs
    the variant is resistant to.
    """
    from modules.variant_db import get_variant_by_id
    variant = get_variant_by_id(variant_id)
    if not variant:
        return {"error": f"Variant {variant_id} not found"}

    protocols = load_treatment_protocols()
    availability = get_nepal_availability()

    resistant = _resistant_drug_set(variant)
    resistant_list = sorted(resistant)
    res_class = _classify_resistance(variant)
    key = _determine_regimen_key(res_class, resistant_list)
    protocol = protocols["protocols"].get(key, protocols["protocols"]["susceptible"])

    primary, alternative = _extract_drugs_from_protocol(protocol)

    contraindicated_in_primary = [d for d in primary if d.lower() in resistant]
    safe_primary = [d for d in primary if d.lower() not in resistant]

    nepal_avail = _get_nepal_avail_map(
        list(dict.fromkeys(primary + alternative)), availability,
    )

    duration = protocol.get(
        "total_duration_weeks",
        protocol.get("total_duration_weeks_shorter", 26),
    )
    monitoring = protocol.get("monitoring", [])

    precursors = _precursor_drug_set(variant)
    drugs_to_monitor = sorted(precursors - resistant)

    interactions = check_drug_interactions(safe_primary)

    return {
        "variant_id": variant_id,
        "resistance_class": res_class,
        "who_category": protocol.get("who_category", ""),
        "recommended_regimen": safe_primary,
        "full_protocol_regimen": primary,
        "alternative_regimen": alternative,
        "contraindicated_drugs": resistant_list,
        "contraindicated_in_primary": contraindicated_in_primary,
        "drugs_to_monitor": drugs_to_monitor,
        "duration_weeks": duration,
        "monitoring_notes": monitoring,
        "drug_interactions": interactions,
        "nepal_availability": nepal_avail,
        "nepal_ntp_compliant": _check_ntp_compliance(safe_primary, nepal_avail),
        "regimen_string": format_regimen_string(safe_primary),
    }


# ---------------------------------------------------------------------------
# 2. infer_treatment_from_similarity
# ---------------------------------------------------------------------------

def infer_treatment_from_similarity(
    comparison: ComparisonResult,
    new_variant: TBVariant,
) -> dict:
    """Infer treatment by starting from the matched variant's protocol
    and adjusting for the new variant's resistance mutations.

    Steps:
      1. Get matched variant's treatment protocol
      2. Check new variant's resistance profile
      3. Remove drugs new variant is resistant to
      4. Substitute alternatives
      5. Flag confidence level
    """
    from modules.variant_db import get_variant_by_id

    matched = get_variant_by_id(comparison.matched_variant_id)
    if matched:
        base = get_treatment_for_known_variant(comparison.matched_variant_id)
    else:
        base = _build_default_protocol(new_variant)

    new_resistant = _resistant_drug_set(new_variant)
    base_resistant = set(base.get("contraindicated_drugs", []))
    additional_resistant = new_resistant - base_resistant

    recommended = base.get("recommended_regimen", [])
    alternative = base.get("alternative_regimen", [])

    newly_removed = [d for d in recommended if d.lower() in additional_resistant]
    safe = [d for d in recommended if d.lower() not in new_resistant]

    substitutions = []
    if newly_removed:
        available_alt = [
            d for d in alternative
            if d.lower() not in new_resistant
        ]
        for removed in newly_removed:
            removed_class = DRUG_CLASSES.get(removed.lower(), "")
            sub = _find_substitute(removed_class, available_alt, safe)
            if sub:
                safe.append(sub)
                available_alt.remove(sub)
                substitutions.append({"removed": removed, "substituted": sub, "reason": f"Resistance to {removed}"})

    warnings: list[str] = []
    if additional_resistant:
        warnings.append(
            f"New variant has additional resistance to: {', '.join(sorted(additional_resistant))}"
        )
    if newly_removed and not substitutions:
        warnings.append("Could not find substitutes for all removed drugs — specialist review required")

    score = comparison.weighted_final_score
    if score > 90:
        confidence = "HIGH"
    elif score >= 75:
        confidence = "MODERATE"
    elif score >= 60:
        confidence = "LOW"
    else:
        confidence = "VERY_LOW"
        warnings.append("Low similarity — treatment inference is unreliable, phenotypic AST required")

    precursors = _precursor_drug_set(new_variant)
    drugs_to_monitor = sorted(precursors - new_resistant)

    interactions = check_drug_interactions(safe)

    availability = get_nepal_availability()
    nepal_avail = _get_nepal_avail_map(
        list(dict.fromkeys(safe + alternative)), availability,
    )

    return {
        "variant_id": new_variant.variant_id,
        "based_on_variant": comparison.matched_variant_id,
        "similarity_score": score,
        "resistance_class": _classify_resistance(new_variant),
        "confidence": confidence,
        "recommended_regimen": safe,
        "alternative_regimen": [d for d in alternative if d.lower() not in new_resistant],
        "contraindicated_drugs": sorted(new_resistant),
        "drugs_to_monitor": drugs_to_monitor,
        "substitutions": substitutions,
        "warnings": warnings,
        "drug_interactions": interactions,
        "nepal_availability": nepal_avail,
        "nepal_ntp_compliant": _check_ntp_compliance(safe, nepal_avail),
        "regimen_string": format_regimen_string(safe),
    }


def _find_substitute(drug_class: str, available: list[str], already_in: list[str]) -> str | None:
    in_lower = {d.lower() for d in already_in}
    same_class = [
        d for d in available
        if DRUG_CLASSES.get(d.lower(), "") == drug_class and d.lower() not in in_lower
    ]
    if same_class:
        return same_class[0]
    if available:
        other = [d for d in available if d.lower() not in in_lower]
        return other[0] if other else None
    return None


def _build_default_protocol(variant: TBVariant) -> dict:
    """Build a protocol from scratch when no matched variant exists."""
    res_class = _classify_resistance(variant)
    resistant = sorted(_resistant_drug_set(variant))
    protocols = load_treatment_protocols()
    key = _determine_regimen_key(res_class, resistant)
    protocol = protocols["protocols"].get(key, protocols["protocols"]["susceptible"])
    primary, alternative = _extract_drugs_from_protocol(protocol)
    return {
        "recommended_regimen": primary,
        "alternative_regimen": alternative,
        "contraindicated_drugs": resistant,
        "monitoring_notes": protocol.get("monitoring", []),
    }


# ---------------------------------------------------------------------------
# 3. generate_treatment_recommendation
# ---------------------------------------------------------------------------

def generate_treatment_recommendation(
    variant: TBVariant,
    resistance_class: str | ComparisonResult | None = None,
    similarity_score: float = 0.0,
) -> dict:
    """Generate a complete treatment recommendation.

    Accepts either a resistance_class string (legacy) or a ComparisonResult.
    Returns a dict with all fields needed by downstream consumers
    (charts, reports, format_treatment_summary).
    """
    if isinstance(resistance_class, ComparisonResult):
        comparison = resistance_class
        res_class = _classify_resistance(variant)
        similarity_score = comparison.weighted_final_score
    else:
        comparison = None
        res_class = resistance_class or _classify_resistance(variant)

    protocols = load_treatment_protocols()
    availability = get_nepal_availability()

    resistant = _resistant_drug_set(variant)
    resistant_list = sorted(resistant)
    key = _determine_regimen_key(res_class, resistant_list)
    protocol = protocols["protocols"].get(key, protocols["protocols"]["susceptible"])

    primary, alternative = _extract_drugs_from_protocol(protocol)

    safe_primary = [d for d in primary if d.lower() not in resistant]
    safe_alternative = [d for d in alternative if d.lower() not in resistant]

    if comparison and comparison.matched_variant_id and comparison.matched_variant_id != "none":
        inferred = infer_treatment_from_similarity(comparison, variant)
        safe_primary = inferred["recommended_regimen"]
        safe_alternative = inferred.get("alternative_regimen", safe_alternative)

    nepal_avail = _get_nepal_avail_map(
        list(dict.fromkeys(primary + alternative + safe_primary + safe_alternative)),
        availability,
    )

    duration = protocol.get(
        "total_duration_weeks",
        protocol.get("total_duration_weeks_shorter", 26),
    )
    monitoring = protocol.get("monitoring", [])

    precursors = _precursor_drug_set(variant)
    drugs_to_monitor = sorted(precursors - resistant)

    warnings: list[str] = []
    unavailable_primary = [
        d for d in safe_primary
        if not nepal_avail.get(d, {}).get("available", False)
    ]
    if unavailable_primary:
        warnings.append(
            f"Drugs not available in Nepal: {', '.join(unavailable_primary)}. "
            f"Contact NTC for procurement or consider alternatives."
        )

    if drugs_to_monitor:
        warnings.append(
            f"Precursor mutations detected for: {', '.join(drugs_to_monitor)}. "
            f"Monitor for resistance emergence during treatment."
        )

    if variant.has_xdr_profile():
        warnings.append("XDR profile — refer to NTC or GENETUP for specialist management")
    elif variant.has_mdr_profile():
        warnings.append("MDR profile — treatment at designated MDR centre required")

    interactions = check_drug_interactions(safe_primary)
    if interactions:
        warnings.extend(interactions)

    confidence = _compute_confidence(variant, nepal_avail, safe_primary, similarity_score)

    next_steps = _generate_next_steps(res_class, confidence, drugs_to_monitor)

    ntp_compliant = _check_ntp_compliance(safe_primary, nepal_avail)

    return {
        "variant_id": variant.variant_id,
        "resistance_class": res_class,
        "recommended_regimen": safe_primary,
        "alternative_regimen": safe_alternative,
        "contraindicated_drugs": resistant_list,
        "drugs_to_avoid": resistant_list,
        "drugs_to_monitor": drugs_to_monitor,
        "duration_weeks": duration,
        "duration_months": round(duration / 4.33),
        "monitoring_notes": monitoring,
        "nepal_availability": nepal_avail,
        "confidence_score": confidence,
        "confidence": _confidence_label(confidence),
        "warnings": warnings,
        "next_steps": next_steps,
        "nepal_ntp_compliant": ntp_compliant,
        "regimen_string": format_regimen_string(safe_primary),
        "who_category": protocol.get("who_category", ""),
    }


def _confidence_label(score: float) -> str:
    if score >= 0.85:
        return "HIGH"
    if score >= 0.65:
        return "MODERATE"
    if score >= 0.40:
        return "LOW"
    return "VERY_LOW"


def _generate_next_steps(
    res_class: str,
    confidence: float,
    drugs_to_monitor: list[str],
) -> list[str]:
    steps = []

    if confidence < 0.65:
        steps.append("Obtain phenotypic antimicrobial susceptibility testing (AST) to confirm resistance profile")

    rc = res_class.upper().replace("-", "_")
    if rc in ("XDR", "PRE_XDR"):
        steps.append("Refer to National Tuberculosis Centre (NTC) for specialist management")
        steps.append("Consider compassionate-use access for unavailable drugs via WHO/GDF")
    elif rc == "MDR":
        steps.append("Register patient in national eDRTB system")
        steps.append("Initiate treatment at designated MDR treatment centre")
    else:
        steps.append("Initiate treatment and schedule month-2 sputum follow-up")

    if drugs_to_monitor:
        steps.append(
            f"Monitor for resistance emergence to: {', '.join(drugs_to_monitor)}. "
            f"Repeat DST if treatment response is poor."
        )

    steps.append("Ensure DOT (Directly Observed Therapy) throughout treatment")
    steps.append("Report case to district TB register and NTP HMIS")

    return steps


def _check_ntp_compliance(drugs: list[str], nepal_avail: dict) -> bool:
    """Check if all recommended drugs are available through NTP channels."""
    for drug in drugs:
        info = nepal_avail.get(drug, {})
        if not info.get("available", False):
            return False
    return True


# ---------------------------------------------------------------------------
# 4. check_drug_interactions
# ---------------------------------------------------------------------------

def check_drug_interactions(drugs: list[str]) -> list[str]:
    """Check for known drug-drug interactions in the regimen.

    Returns a list of warning strings for any detected interactions.
    """
    drug_set = {d.lower() for d in drugs}
    warnings: list[str] = []

    for interaction in KNOWN_INTERACTIONS:
        pair = {d.lower() for d in interaction["drugs"]}
        if pair.issubset(drug_set):
            severity = interaction["severity"].upper()
            warnings.append(f"[{severity}] {interaction['warning']}")

    return warnings


# ---------------------------------------------------------------------------
# 5. format_regimen_string
# ---------------------------------------------------------------------------

def format_regimen_string(drugs: list[str]) -> str:
    """Format a drug list using standard TB shorthand notation.

    For first-line regimens: 2HRZE/4HR format
    For second-line regimens: lists abbreviations
    """
    if not drugs:
        return ""

    drugs_lower = [d.lower() for d in drugs]

    first_line_intensive = ["isoniazid", "rifampicin", "pyrazinamide", "ethambutol"]
    first_line_continuation = ["isoniazid", "rifampicin"]

    if set(drugs_lower) == set(first_line_intensive):
        return "2HRZE/4HR"

    if set(first_line_intensive).issubset(set(drugs_lower)):
        extra = [d for d in drugs_lower if d not in first_line_intensive]
        extra_abbr = "".join(DRUG_ABBREVIATIONS.get(d, d[:3].title()) for d in extra)
        if extra_abbr:
            return f"2HRZE+{extra_abbr}/4HR+{extra_abbr}"
        return "2HRZE/4HR"

    has_intensive = "pyrazinamide" in drugs_lower or "ethambutol" in drugs_lower
    has_rif = "rifampicin" in drugs_lower
    has_inh = "isoniazid" in drugs_lower

    if has_rif and has_inh and has_intensive:
        intensive_drugs = []
        continuation_drugs = []
        for d in drugs_lower:
            abbr = DRUG_ABBREVIATIONS.get(d, d[:3].title())
            intensive_drugs.append(abbr)
            if d not in ("pyrazinamide", "ethambutol"):
                continuation_drugs.append(abbr)
        return f"2{''.join(intensive_drugs)}/4{''.join(continuation_drugs)}"

    abbrs = []
    for d in drugs_lower:
        abbrs.append(DRUG_ABBREVIATIONS.get(d, d[:3].title()))
    return " + ".join(abbrs)


# ---------------------------------------------------------------------------
# Treatment centres
# ---------------------------------------------------------------------------

def get_treatment_centers(resistance_class: str) -> list[dict]:
    """Return Nepal treatment centres appropriate for the resistance class."""
    protocols = load_treatment_protocols()
    centers = protocols.get("nepal_treatment_centers", [])

    rc = resistance_class.upper().replace("-", "_")
    if rc in ("XDR", "PRE_XDR"):
        return [c for c in centers if "Pre-XDR/XDR treatment" in c.get("services", [])]
    if rc == "MDR":
        return [c for c in centers if "MDR treatment" in c.get("services", [])]
    return centers


# ---------------------------------------------------------------------------
# Legacy public API
# ---------------------------------------------------------------------------

def determine_regimen_key(resistance_class: str, resistant_drugs: list[str]) -> str:
    """Legacy alias."""
    return _determine_regimen_key(resistance_class, resistant_drugs)


def format_treatment_summary(rec: dict) -> str:
    """Format a treatment recommendation dict as readable text."""
    lines = [
        f"=== Treatment Recommendation for {rec['variant_id']} ===",
        f"Resistance Class: {rec['resistance_class']}",
    ]

    if "confidence_score" in rec:
        lines.append(f"Confidence: {rec['confidence_score'] * 100:.1f}%")
    elif "confidence" in rec:
        lines.append(f"Confidence: {rec['confidence']}")

    if rec.get("regimen_string"):
        lines.append(f"Regimen: {rec['regimen_string']}")

    lines.append("")
    lines.append("Recommended Regimen:")
    for drug in rec.get("recommended_regimen", []):
        avail = rec.get("nepal_availability", {}).get(drug, {})
        status = "Available" if avail.get("available") else "NOT AVAILABLE"
        lines.append(f"  - {drug} [{status}]")

    alt = rec.get("alternative_regimen", [])
    if alt:
        lines.append("")
        lines.append("Alternative Regimen:")
        for drug in alt:
            avail = rec.get("nepal_availability", {}).get(drug, {})
            status = "Available" if avail.get("available") else "NOT AVAILABLE"
            lines.append(f"  - {drug} [{status}]")

    lines.append("")
    lines.append(f"Duration: {rec.get('duration_weeks', '?')} weeks")

    contra = rec.get("contraindicated_drugs", rec.get("drugs_to_avoid", []))
    if contra:
        lines.append("")
        lines.append("Contraindicated (resistant):")
        for drug in contra:
            lines.append(f"  - {drug}")

    monitor_drugs = rec.get("drugs_to_monitor", [])
    if monitor_drugs:
        lines.append("")
        lines.append("Monitor for resistance emergence:")
        for drug in monitor_drugs:
            lines.append(f"  - {drug}")

    monitoring = rec.get("monitoring_notes", [])
    if monitoring:
        lines.append("")
        lines.append("Monitoring:")
        for note in monitoring:
            lines.append(f"  * {note}")

    warnings = rec.get("warnings", [])
    if warnings:
        lines.append("")
        lines.append("Warnings:")
        for w in warnings:
            lines.append(f"  ! {w}")

    next_steps = rec.get("next_steps", [])
    if next_steps:
        lines.append("")
        lines.append("Next Steps:")
        for i, step in enumerate(next_steps, 1):
            lines.append(f"  {i}. {step}")

    if rec.get("nepal_ntp_compliant") is not None:
        lines.append("")
        ntp = "Yes" if rec["nepal_ntp_compliant"] else "No — some drugs unavailable through NTP"
        lines.append(f"Nepal NTP Compliant: {ntp}")

    return "\n".join(lines)
