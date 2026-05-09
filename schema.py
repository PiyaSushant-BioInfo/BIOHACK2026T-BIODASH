"""
TBAnalytica Shared Schema
All data structures used across modules. Pydantic models with
validation, type hints, and JSON serialization.
"""

from __future__ import annotations

import json
from datetime import datetime
from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field, field_validator, model_validator


# ---------------------------------------------------------------------------
# Enums
# ---------------------------------------------------------------------------

class DrugSensitivity(str, Enum):
    RESISTANT = "resistant"
    SENSITIVE = "sensitive"
    INTERMEDIATE = "intermediate"


class ResistanceLevel(str, Enum):
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    CRITICAL = "critical"


class ConfidenceLevel(str, Enum):
    HIGH = "HIGH"
    MODERATE = "MODERATE"
    LOW = "LOW"


class RiskLevel(str, Enum):
    LOW = "LOW"
    MODERATE = "MODERATE"
    HIGH = "HIGH"
    CRITICAL = "CRITICAL"


class RiskColor(str, Enum):
    GREEN = "green"
    YELLOW = "yellow"
    RED = "red"
    BLACK = "black"


class DataSource(str, Enum):
    WHO_CATALOGUE = "WHO Mutation Catalogue 2022"
    NCBI_REFSEQ = "NCBI RefSeq"
    NCBI_GENBANK = "NCBI GenBank"
    UNIPROT_REVIEWED = "UniProt/SwissProt (manually reviewed)"
    UNIPROT_UNREVIEWED = "UniProt/TrEMBL (auto-annotated)"
    PDB = "Protein Data Bank"
    ALPHAFOLD = "AlphaFold (predicted)"
    LOCAL_DB = "Local Database"


# ---------------------------------------------------------------------------
# DataQualityScore
# ---------------------------------------------------------------------------

class DataQualityScore(BaseModel):
    """Quality assessment for a single piece of fetched data."""

    source: DataSource
    raw_score: float = Field(default=0.0, ge=0.0, le=100.0)
    source_score: float = Field(default=0.0, ge=0.0, le=30.0)
    citation_score: float = Field(default=0.0, ge=0.0, le=20.0)
    completeness_score: float = Field(default=0.0, ge=0.0, le=40.0)
    recency_score: float = Field(default=-5.0, le=20.0)
    review_status: str = Field(default="unknown")
    confidence: str = Field(default="LOW")
    use_for_analysis: bool = Field(default=False)
    details: dict = Field(default_factory=dict)

    def compute(self) -> DataQualityScore:
        self.raw_score = round(min(
            self.source_score + self.citation_score
            + self.completeness_score + self.recency_score,
            100.0,
        ), 2)
        if self.raw_score < 0:
            self.raw_score = 0.0
        if self.raw_score >= 75:
            self.confidence = "HIGH"
        elif self.raw_score >= 60:
            self.confidence = "MODERATE"
        elif self.raw_score >= 40:
            self.confidence = "LOW"
        else:
            self.confidence = "REJECT"
        self.use_for_analysis = self.raw_score >= 40
        return self

    def summary_line(self, label: str = "") -> str:
        icon = {
            "HIGH": "[OK]",
            "MODERATE": "[!!]",
            "LOW": "[??]",
            "REJECT": "[XX]",
        }.get(self.confidence, "[??]")
        return (
            f"{icon} {label}: {self.source.value} | "
            f"Quality: {self.raw_score:.0f}/100 | {self.confidence}"
        )

    def to_json(self) -> str:
        return self.model_dump_json(indent=2)

    @classmethod
    def from_json(cls, raw: str | dict) -> DataQualityScore:
        if isinstance(raw, str):
            raw = json.loads(raw)
        return cls.model_validate(raw)


QUALITY_GATE_THRESHOLD = 60


# ---------------------------------------------------------------------------
# 3. Mutation
# ---------------------------------------------------------------------------

class Mutation(BaseModel):
    """Single nucleotide / codon-level mutation in a resistance gene."""

    position: int = Field(..., description="1-based codon or nucleotide position (negative for promoter regions)")
    reference_codon: str = Field(default="", description="Wild-type codon (e.g. AGC)")
    mutant_codon: str = Field(default="", description="Mutant codon (e.g. ACC)")
    reference_amino_acid: str = Field(default="", max_length=3, description="Wild-type amino acid one- or three-letter")
    mutant_amino_acid: str = Field(default="", max_length=3, description="Mutant amino acid one- or three-letter")
    is_synonymous: bool = Field(default=False, description="True when the nucleotide change does not alter the amino acid")
    is_resistance_conferring: bool = Field(default=False, description="True when this mutation is known to confer drug resistance")
    drug_affected: str = Field(default="", description="Primary drug whose efficacy is reduced by this mutation")
    resistance_level: ResistanceLevel = Field(default=ResistanceLevel.LOW, description="Severity of resistance conferred")
    one_step_away_risk: bool = Field(default=False, description="True when a single additional SNP could confer resistance")
    one_step_away_drug: Optional[str] = Field(default=None, description="Drug at risk if the one-step mutation occurs")

    @field_validator("reference_codon", "mutant_codon")
    @classmethod
    def codon_must_be_valid_chars(cls, v: str) -> str:
        if v and not all(c in "ACGTUacgtu" for c in v):
            raise ValueError(f"Codon contains invalid nucleotide characters: {v}")
        return v.upper()

    @field_validator("reference_amino_acid", "mutant_amino_acid")
    @classmethod
    def amino_acid_uppercase(cls, v: str) -> str:
        return v.upper() if v else v

    @model_validator(mode="after")
    def set_synonymous_flag(self) -> Mutation:
        if self.reference_amino_acid and self.mutant_amino_acid:
            if self.reference_amino_acid == self.mutant_amino_acid and not self.is_synonymous:
                self.is_synonymous = True
        return self

    def short_code(self) -> str:
        if self.reference_amino_acid and self.mutant_amino_acid:
            return f"{self.reference_amino_acid}{self.position}{self.mutant_amino_acid}"
        return f"{self.reference_codon}{self.position}{self.mutant_codon}"

    def to_json(self) -> str:
        return self.model_dump_json(indent=2)

    @classmethod
    def from_json(cls, raw: str | dict) -> Mutation:
        if isinstance(raw, str):
            raw = json.loads(raw)
        return cls.model_validate(raw)


# ---------------------------------------------------------------------------
# 2. ResistanceGene
# ---------------------------------------------------------------------------

class ResistanceGene(BaseModel):
    """A gene known to harbour TB drug-resistance mutations."""

    gene_name: str = Field(..., description="Gene symbol (e.g. rpoB, katG, inhA, gyrA, gyrB, pncA)")
    mutations: list[Mutation] = Field(default_factory=list, description="Known mutations in this gene")
    drug_target: str = Field(default="", description="Drug or drug class this gene product is the target of")

    @field_validator("gene_name")
    @classmethod
    def gene_name_not_empty(cls, v: str) -> str:
        if not v.strip():
            raise ValueError("gene_name must not be empty")
        return v.strip()

    def resistance_mutations(self) -> list[Mutation]:
        return [m for m in self.mutations if m.is_resistance_conferring]

    def silent_mutations(self) -> list[Mutation]:
        return [m for m in self.mutations if m.is_synonymous]

    def one_step_risks(self) -> list[Mutation]:
        return [m for m in self.mutations if m.one_step_away_risk]

    def to_json(self) -> str:
        return self.model_dump_json(indent=2)

    @classmethod
    def from_json(cls, raw: str | dict) -> ResistanceGene:
        if isinstance(raw, str):
            raw = json.loads(raw)
        return cls.model_validate(raw)


# ---------------------------------------------------------------------------
# 1. TBVariant
# ---------------------------------------------------------------------------

class TBVariant(BaseModel):
    """Complete representation of a Mycobacterium tuberculosis variant."""

    variant_id: str = Field(..., description="Unique identifier for this variant")
    lineage: str = Field(default="", description="Phylogenetic lineage (1, 2, 3, 4, Beijing, CAS, etc.)")
    name: str = Field(default="", description="Human-readable strain or variant name")
    drug_resistance: dict[str, DrugSensitivity] = Field(
        default_factory=dict,
        description="Mapping of drug name to resistance status",
    )
    resistance_genes: list[ResistanceGene] = Field(
        default_factory=list,
        description="Genes harboring resistance-associated mutations",
    )
    protein_sequences: dict[str, str] = Field(
        default_factory=dict,
        description="Mapping of protein name to amino-acid sequence",
    )
    nucleotide_sequences: dict[str, str] = Field(
        default_factory=dict,
        description="Mapping of gene name to nucleotide sequence",
    )
    treatment_protocol: str = Field(default="", description="Reference to a treatment protocol ID")
    source: str = Field(default="", description="NCBI accession, UniProt ID, or other data source")
    last_updated: datetime = Field(default_factory=datetime.now, description="Timestamp of last record update")

    @field_validator("variant_id")
    @classmethod
    def variant_id_not_empty(cls, v: str) -> str:
        if not v.strip():
            raise ValueError("variant_id must not be empty")
        return v.strip()

    def all_mutations(self) -> list[Mutation]:
        return [m for gene in self.resistance_genes for m in gene.mutations]

    def resistant_drugs(self) -> list[str]:
        return [drug for drug, status in self.drug_resistance.items() if status == DrugSensitivity.RESISTANT]

    def susceptible_drugs(self) -> list[str]:
        return [drug for drug, status in self.drug_resistance.items() if status == DrugSensitivity.SENSITIVE]

    def has_mdr_profile(self) -> bool:
        resistant = {d.lower() for d in self.resistant_drugs()}
        return "isoniazid" in resistant and "rifampicin" in resistant

    def has_xdr_profile(self) -> bool:
        resistant = {d.lower() for d in self.resistant_drugs()}
        mdr = "isoniazid" in resistant and "rifampicin" in resistant
        fq = any(d in resistant for d in ["levofloxacin", "moxifloxacin"])
        injectable = any(d in resistant for d in ["amikacin", "kanamycin", "capreomycin"])
        return mdr and fq and injectable

    def to_json(self) -> str:
        return self.model_dump_json(indent=2)

    def to_dict(self) -> dict:
        return self.model_dump(mode="json")

    @classmethod
    def from_json(cls, raw: str | dict) -> TBVariant:
        if isinstance(raw, str):
            raw = json.loads(raw)
        return cls.model_validate(raw)


# ---------------------------------------------------------------------------
# 4. ComparisonResult
# ---------------------------------------------------------------------------

class ComparisonResult(BaseModel):
    """Result of comparing a new/unknown variant against a known reference."""

    new_variant_id: str = Field(..., description="ID of the query (new/unknown) variant")
    matched_variant_id: str = Field(default="", description="ID of the closest matching known variant")
    protein_similarity_score: float = Field(default=0.0, ge=0.0, le=100.0, description="Protein-level sequence identity (%)")
    gene_similarity_score: float = Field(default=0.0, ge=0.0, le=100.0, description="Gene-level sequence identity (%)")
    weighted_final_score: float = Field(default=0.0, ge=0.0, le=100.0, description="Combined weighted similarity score")
    confidence_level: ConfidenceLevel = Field(default=ConfidenceLevel.LOW, description="Confidence in the match")
    gene_changes: list[Mutation] = Field(default_factory=list, description="All detected gene-level changes")
    silent_mutations: list[Mutation] = Field(default_factory=list, description="Synonymous mutations detected")
    resistance_mutations: list[Mutation] = Field(default_factory=list, description="Known resistance-conferring mutations")
    novel_mutations: list[Mutation] = Field(default_factory=list, description="Mutations not present in any known variant")
    treatment_recommendation: str = Field(default="", description="Treatment protocol ID or summary recommendation")

    @model_validator(mode="after")
    def compute_weighted_score(self) -> ComparisonResult:
        if self.weighted_final_score == 0.0 and (self.protein_similarity_score or self.gene_similarity_score):
            self.weighted_final_score = round(
                0.6 * self.protein_similarity_score + 0.4 * self.gene_similarity_score, 2
            )
        return self

    @model_validator(mode="after")
    def derive_confidence(self) -> ComparisonResult:
        if self.confidence_level == ConfidenceLevel.LOW and self.weighted_final_score > 0:
            if self.weighted_final_score >= 90:
                self.confidence_level = ConfidenceLevel.HIGH
            elif self.weighted_final_score >= 70:
                self.confidence_level = ConfidenceLevel.MODERATE
        return self

    def to_json(self) -> str:
        return self.model_dump_json(indent=2)

    def to_dict(self) -> dict:
        return self.model_dump(mode="json")

    @classmethod
    def from_json(cls, raw: str | dict) -> ComparisonResult:
        if isinstance(raw, str):
            raw = json.loads(raw)
        return cls.model_validate(raw)


# ---------------------------------------------------------------------------
# 5. RiskScore
# ---------------------------------------------------------------------------

class RiskScore(BaseModel):
    """Quantified resistance-risk assessment for a variant."""

    variant_id: str = Field(..., description="Variant this score belongs to")
    score: float = Field(default=0.0, ge=0.0, le=100.0, description="Numeric risk score 0-100")
    level: RiskLevel = Field(default=RiskLevel.LOW, description="Categorical risk level")
    color: RiskColor = Field(default=RiskColor.GREEN, description="Color code for UI display")
    factors: list[str] = Field(default_factory=list, description="Individual factors contributing to the score")

    @model_validator(mode="after")
    def sync_level_and_color(self) -> RiskScore:
        if self.score < 25:
            self.level = RiskLevel.LOW
            self.color = RiskColor.GREEN
        elif self.score < 50:
            self.level = RiskLevel.MODERATE
            self.color = RiskColor.YELLOW
        elif self.score < 75:
            self.level = RiskLevel.HIGH
            self.color = RiskColor.RED
        else:
            self.level = RiskLevel.CRITICAL
            self.color = RiskColor.BLACK
        return self

    def to_json(self) -> str:
        return self.model_dump_json(indent=2)

    def to_dict(self) -> dict:
        return self.model_dump(mode="json")

    @classmethod
    def from_json(cls, raw: str | dict) -> RiskScore:
        if isinstance(raw, str):
            raw = json.loads(raw)
        return cls.model_validate(raw)


# ---------------------------------------------------------------------------
# 6. ClinicalReport
# ---------------------------------------------------------------------------

class ClinicalReport(BaseModel):
    """Full clinical report combining variant analysis, risk assessment,
    and treatment guidance.  Serves both the doctor-facing and
    patient-facing report views."""

    patient_id: str = Field(..., description="Anonymised patient identifier")
    variant: TBVariant = Field(..., description="Identified TB variant")
    comparison_result: Optional[ComparisonResult] = Field(default=None, description="Comparison with closest known variant")
    risk_score: RiskScore = Field(..., description="Computed risk assessment")
    treatment: str = Field(default="", description="Treatment protocol ID or regimen summary")
    doctor_report: str = Field(default="", description="Narrative clinical report for the physician")
    patient_report: str = Field(default="", description="Simplified report for the patient")
    generated_at: datetime = Field(default_factory=datetime.now, description="Report generation timestamp")
    data_quality: dict[str, DataQualityScore] = Field(
        default_factory=dict,
        description="Quality scores per data source (label -> DataQualityScore)",
    )
    quality_warnings: list[str] = Field(
        default_factory=list,
        description="Data quality warnings — sources below threshold, gating decisions",
    )
    treatment_gated: bool = Field(
        default=False,
        description="True when treatment was restricted due to low data quality",
    )

    @field_validator("patient_id")
    @classmethod
    def patient_id_not_empty(cls, v: str) -> str:
        if not v.strip():
            raise ValueError("patient_id must not be empty")
        return v.strip()

    def to_json(self) -> str:
        return self.model_dump_json(indent=2)

    def to_dict(self) -> dict:
        return self.model_dump(mode="json")

    @classmethod
    def from_json(cls, raw: str | dict) -> ClinicalReport:
        if isinstance(raw, str):
            raw = json.loads(raw)
        return cls.model_validate(raw)
