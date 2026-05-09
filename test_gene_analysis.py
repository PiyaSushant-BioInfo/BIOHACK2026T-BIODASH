"""Integration test for gene_analysis module."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from modules.gene_analysis import (
    translate_sequence, find_all_mutations, classify_mutation,
    check_resistance_database, check_one_step_away,
    generate_gene_report, full_gene_analysis,
    find_silent_mutations, assess_silent_mutation_risk,
    perform_gene_analysis, compare_gene_sequences,
    translate_codon, identify_mutation_type,
    compute_gc_content, compute_codon_usage,
)
from schema import Mutation, ResistanceLevel, TBVariant, ResistanceGene
from modules.variant_db import load_known_variants


def test_translate_sequence():
    print("--- 1. translate_sequence ---")
    seq = "ATGAAATTT"
    aa = translate_sequence(seq)
    assert aa == "MKF", f"Expected MKF, got {aa}"
    print(f"  ATG AAA TTT -> {aa}")

    seq2 = "ATGTAAGGG"
    aa2 = translate_sequence(seq2)
    assert aa2 == "M*G", f"Expected M*G, got {aa2}"
    print(f"  ATG TAA GGG -> {aa2} (stop codon at pos 2)")

    all_codons = "TTTTTCTTATTGCTTCTCCTACTGATTATCATAATGGTCGTCGTAGTGTCTTCCTCATCGCCTCCCCCACCGACTACCACAACGGCTGCCGCAGCGTATTAGTACTAATAGTAATAGCAACAGCATCATCACAAACAGGCTGCCGCAGCGAATAGAAGAGGAGGATTGATGACGAATGTGTCTGGTGGCGTCGCCGACGGAGTAGCAGAAGGGGCGGTGGAGGG"
    aa_all = translate_sequence(all_codons)
    assert len(aa_all) > 50
    assert "*" in aa_all
    print(f"  All 64 codons translated: {len(aa_all)} AAs, contains stop: {'*' in aa_all}")

    assert translate_sequence("") == ""
    assert translate_sequence("AT") == ""
    print("  Edge cases OK")


def test_find_all_mutations():
    print("--- 2. find_all_mutations ---")
    ref = "TCG" * 530 + "TCG" + "GGG" * 50
    new = "TCG" * 530 + "TTG" + "GGG" * 50
    mutations = find_all_mutations(new, ref, "rpoB")
    assert len(mutations) == 1
    m = mutations[0]
    assert m.position == 531
    assert m.reference_codon == "TCG"
    assert m.mutant_codon == "TTG"
    assert m.reference_amino_acid == "S"
    assert m.mutant_amino_acid == "L"
    assert not m.is_synonymous
    assert m.is_resistance_conferring
    assert "rifampicin" in m.drug_affected
    print(f"  Found {m.short_code()} at pos {m.position}: resistance={m.is_resistance_conferring}, drug={m.drug_affected}")

    ref_silent = "TCG" * 513 + "TTC" + "AAA" * 50
    new_silent = "TCG" * 513 + "TTT" + "AAA" * 50
    silent_muts = find_all_mutations(new_silent, ref_silent, "rpoB")
    assert len(silent_muts) == 1
    sm = silent_muts[0]
    assert sm.position == 514
    assert sm.is_synonymous
    assert sm.reference_amino_acid == sm.mutant_amino_acid == "F"
    print(f"  Silent: {sm.short_code()} at pos {sm.position}: synonymous={sm.is_synonymous}, one_step={sm.one_step_away_risk}")

    same = find_all_mutations("ATGATGATG", "ATGATGATG", "rpoB")
    assert len(same) == 0
    print("  Identical sequences: 0 mutations")


def test_classify_mutation():
    print("--- 3. classify_mutation ---")
    missense = Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L")
    classified = classify_mutation(missense)
    assert not classified.is_synonymous
    print(f"  S531L -> synonymous={classified.is_synonymous}")

    silent = Mutation(position=514, reference_amino_acid="F", mutant_amino_acid="F",
                      reference_codon="TTC", mutant_codon="TTT")
    classified_s = classify_mutation(silent)
    assert classified_s.is_synonymous
    print(f"  F514F -> synonymous={classified_s.is_synonymous}")

    nonsense = Mutation(position=100, reference_amino_acid="W", mutant_amino_acid="*")
    classified_n = classify_mutation(nonsense)
    assert classified_n.is_resistance_conferring
    print(f"  W100* -> resistance={classified_n.is_resistance_conferring} (nonsense)")


def test_check_resistance_database():
    print("--- 4. check_resistance_database ---")
    result = check_resistance_database("rpoB", 531, "S531L")
    assert result["is_known_resistance"] is True
    assert "rifampicin" in result["drug_affected"]
    assert result["resistance_level"] == "high"
    print(f"  rpoB S531L: known={result['is_known_resistance']}, drugs={result['drug_affected']}, level={result['resistance_level']}")

    result2 = check_resistance_database("katG", 315, "S315T")
    assert result2["is_known_resistance"] is True
    assert "isoniazid" in result2["drug_affected"]
    print(f"  katG S315T: known={result2['is_known_resistance']}, drugs={result2['drug_affected']}")

    result3 = check_resistance_database("rpoB", 999, "X999Y")
    assert result3["is_known_resistance"] is False
    print(f"  rpoB X999Y: known={result3['is_known_resistance']} (novel)")

    result4 = check_resistance_database("gyrA", 94, "D94G")
    assert result4["is_known_resistance"] is True
    assert result4["in_critical_region"] is True
    print(f"  gyrA D94G: critical_region={result4['in_critical_region']}")


def test_check_one_step_away():
    print("--- 5. check_one_step_away ---")
    result = check_one_step_away("rpoB", 531, "TCA")
    print(f"  rpoB pos 531, codon TCA: precursor={result['is_precursor']}, risk={result['risk_level']}")
    if result["target_mutations"]:
        for t in result["target_mutations"]:
            print(f"    -> {t['mutation_code']} ({t['target_codon']}): drugs={t['drugs']}, nt_change={t['nucleotide_change']}")

    result2 = check_one_step_away("rpoB", 514, "TTT")
    print(f"  rpoB pos 514, codon TTT: precursor={result2['is_precursor']}, risk={result2['risk_level']}")

    result3 = check_one_step_away("rpoB", 100, "GGG")
    assert result3["is_precursor"] is False
    print(f"  rpoB pos 100 (not near resistance): precursor={result3['is_precursor']}")


def test_generate_gene_report():
    print("--- 6. generate_gene_report ---")
    mutations = [
        Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                 reference_codon="TCG", mutant_codon="TTG",
                 is_resistance_conferring=True, drug_affected="rifampicin",
                 resistance_level=ResistanceLevel.HIGH),
        Mutation(position=514, reference_amino_acid="F", mutant_amino_acid="F",
                 reference_codon="TTC", mutant_codon="TTT",
                 is_synonymous=True, one_step_away_risk=True,
                 one_step_away_drug="rifampicin"),
        Mutation(position=400, reference_amino_acid="A", mutant_amino_acid="A",
                 reference_codon="GCT", mutant_codon="GCC",
                 is_synonymous=True),
        Mutation(position=200, reference_amino_acid="G", mutant_amino_acid="D",
                 reference_codon="GGC", mutant_codon="GAC"),
    ]
    report = generate_gene_report(mutations, "rpoB")
    assert report["total_mutations"] == 4
    assert report["resistance_count"] == 1
    assert report["silent_precursor_count"] == 1
    assert report["silent_benign_count"] == 1
    assert report["novel_count"] == 1
    assert report["overall_risk"] == "high"
    print(f"  rpoB: {report['total_mutations']} mutations")
    print(f"    resistance: {report['resistance_count']}")
    print(f"    silent precursor: {report['silent_precursor_count']}")
    print(f"    silent benign: {report['silent_benign_count']}")
    print(f"    novel: {report['novel_count']}")
    print(f"    risk: {report['overall_risk']}")


def test_full_gene_analysis():
    print("--- 7. full_gene_analysis ---")
    known = load_known_variants()
    new_var = TBVariant(
        variant_id="TEST_NEW",
        resistance_genes=[
            ResistanceGene(gene_name="rpoB", mutations=[
                Mutation(position=531, reference_amino_acid="S", mutant_amino_acid="L",
                         is_resistance_conferring=True, drug_affected="rifampicin",
                         resistance_level=ResistanceLevel.HIGH),
            ], drug_target="rifampicin"),
            ResistanceGene(gene_name="katG", mutations=[
                Mutation(position=315, reference_amino_acid="S", mutant_amino_acid="T",
                         is_resistance_conferring=True, drug_affected="isoniazid",
                         resistance_level=ResistanceLevel.HIGH),
            ], drug_target="isoniazid"),
        ],
    )
    ref_var = known[0]
    result = full_gene_analysis(new_var, ref_var)
    assert result["new_variant_id"] == "TEST_NEW"
    assert result["genes_analyzed"] > 0
    print(f"  Analyzed {result['genes_analyzed']} genes")
    s = result["summary"]
    print(f"    Total mutations: {s['total_mutations']}")
    print(f"    Resistance: {s['resistance_mutations']}")
    print(f"    Silent: {s['silent_mutations']}")
    print(f"    Precursor: {s['precursor_mutations']}")
    print(f"    Novel: {s['novel_mutations']}")
    print(f"    Drugs affected: {s['drugs_affected']}")
    print(f"    Overall risk: {s['overall_risk']}")
    if result["warnings"]:
        for w in result["warnings"]:
            print(f"    WARNING: {w}")


def test_legacy_functions():
    print("--- 8. Legacy functions ---")
    silent = find_silent_mutations("TCGTCGTTT", "TCGTCGTTC")
    assert len(silent) == 1
    assert silent[0]["amino_acid"] == "F"
    print(f"  find_silent_mutations: {len(silent)} found")

    assessed = assess_silent_mutation_risk("rpoB", [
        {"codon_position": 514, "ref_codon": "TTC", "alt_codon": "TTT",
         "amino_acid": "F", "nucleotide_position": 1540}
    ])
    assert len(assessed) == 1
    assert assessed[0]["risk_level"] != ""
    print(f"  assess_silent_mutation_risk: risk={assessed[0]['risk_level']}, note={assessed[0]['note']}")

    assert translate_codon("ATG") == "M"
    assert translate_codon("TAA") == "*"
    assert identify_mutation_type("TCG", "TTG") == "missense"
    assert identify_mutation_type("TTC", "TTT") == "silent"
    assert identify_mutation_type("TGG", "TGA") == "nonsense"
    print("  translate_codon, identify_mutation_type OK")

    gc = compute_gc_content("GCGCATAT")
    assert 40 < gc < 60
    print(f"  GC content of GCGCATAT: {gc}%")

    usage = compute_codon_usage("ATGATGATG")
    assert "ATG" in usage
    print(f"  Codon usage: {usage}")

    comp = compare_gene_sequences("ATGATGATG", "ATGATTATG")
    assert comp["total_snps"] > 0
    print(f"  compare_gene_sequences: {comp['total_snps']} SNPs, identity={comp['identity_percent']}%")


if __name__ == "__main__":
    test_translate_sequence()
    test_find_all_mutations()
    test_classify_mutation()
    test_check_resistance_database()
    test_check_one_step_away()
    test_generate_gene_report()
    test_full_gene_analysis()
    test_legacy_functions()
    print()
    print("ALL TESTS PASSED")
