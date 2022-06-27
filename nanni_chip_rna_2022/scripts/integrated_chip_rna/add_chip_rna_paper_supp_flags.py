#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Add flags for ChIP-RNA supplemental files"
    )

    # Input data
    parser.add_argument(
        "-m",
        "--mel",
        dest="inM",
        required=True,
        help=(
            "Input CSV of D. melanogaster gene-level file with combo flag "
            "gene_sex_bias_ttest_foldchange and gene list flags"
        )
    )
    parser.add_argument(
        "-s",
        "--sim",
        dest="inS",
        required=True,
        help=(
            "Input CSV of D. simulans gene-level file with combo flag "
            "gene_sex_bias_ttest_foldchange and gene list flags"
        )
    )
    parser.add_argument(
        "-or",
        "--ortho",
        dest="inO",
        required=True,
        help=(
            "Input CSV of D. melangoster and D. simulans orthologous "
            "gene-level file with ttest+foldchange values and gene list flags"
        )
    )

    # Output data
    parser.add_argument(
        "-d",
        "--output-dir",
        dest="outDir",
        required=True,
        help="Output directory for files with new flags"
    )

    args = parser.parse_args()
    return args

def main():
    # Get inputs
    melDF = pd.read_csv(args.inM, low_memory=False)
    simDF = pd.read_csv(args.inS, low_memory=False)
    orthoDF = pd.read_csv(args.inO, low_memory=False)

    # Create CHiP flags:
    # flag_any_k4 = detected k4 in any sex (gene_k4 != "none")
    # flag_has_male_k4 = detected k4 in males (gene_k4 != "none" or "fem")
    # flag_has_female_k4 = detected k4 in males (gene_k4 != "none" or "male")
    # flag_male_limited_k4 = detected k4 in only male (gene_k4 == "male")
    # flag_female_limited_k4 = detected k4 in only female (gene_k4 == "fem")
    # flag_sex_limited_k4 = detected k4 in only male or only female (gene_k4 == "male" or "fem")
    # flag_any_k27 = detected k27 in any sex (gene_k27 != "none")
    # flag_has_male_k27 = detected k27 in males (gene_k27 != "none" or "fem")
    # flag_has_female_k27 = detected k27 in males (gene_k27 != "none" or "male")
    # flag_male_limited_k27 = detected k27 in only male (gene_k27 == "male")
    # flag_female_limited_k27 = detected k27 in only female (gene_k27 == "fem")
    # flag_sex_limited_k27 = detected k27 in only male or only female (gene_k27 == "male" or "fem")
    
    # Create RNA flags:
    #
    # flag_expressed = detected (APN>0) in at least one sex (gene_rna != "none")
    # flag_expressed_M = detected (APN>0) in males (gene_rna != "none" or "fem")
    # flag_expressed_F = detected (APN>0) in females (gene_rna != "none" or "male")
    # flag_sex_limited = detected (APN>0) in only one sex (gene_rna = "fem" or "male")
    # flag_male_limited = detected (APN>0) in only males (gene_rna = "male")
    # flag_female_limited = detected (APN>0) in only females (gene_rna = "fem")
    # flag_low_expressed = APN < 5
    # flag_sex_biased = detected in both sexes with statistical evidence of differential expression (gene_trend_ttest = male, female, or male_and_female)
    # flag_U = detected in both sexes with no statistical evidence of differential expression
    # flag_F = detected in both sexes with statistical significance with a trend towards female
    # flag_FP = detected in both sexes with flag_F + at least 2-fold change towards female
    # flag_M = detected in both sexes with statistical significance with a trend towards male
    # flag_MP = detected in both sexes with flag_M + at least 2-fold change towards male
    # flag_MF = detected in both sexes with at least one exonic region statistically significant in males and at least one in females
    # flag_MFP = detected in both sexes with flag_MF + at least 2-fold change in each of the sex-biased exonic regions
    # flag_conserved_sex_bias = detected in both sexes with statistical evidence of differential expression in both sexes in the same direction
    # flag_one2one_ortholog = mel_geneID only ortholog to one sim_geneID and vice versa (excluding paralogs)

### CHiP flags
    # flag_any_k4, flag_any_k4_mel, flag_any_k4_sim
    melDF['flag_any_k4'] = np.where(melDF['gene_k4']!="none", 1, 0)
    simDF['flag_any_k4'] = np.where(simDF['gene_k4']!="none", 1, 0)
    orthoDF['flag_any_k4_mel'] = np.where(orthoDF['gene_k4']!="none", 1, 0)
    orthoDF['flag_any_k4_sim'] = np.where(orthoDF['sim_gene_k4']!="none", 1, 0)
    orthoDF['flag_any_k4_both_species'] = np.where(
            (orthoDF['flag_any_k4_mel'] == 1)
            & (orthoDF['flag_any_k4_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_any_k4_either_species'] = np.where(
            (orthoDF['flag_any_k4_mel'] == 1)
            | (orthoDF['flag_any_k4_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_any_k4_mel_only'] = np.where(
            (orthoDF['flag_any_k4_mel'] == 1)
            & (orthoDF['flag_any_k4_sim'] == 0),
            1,
            0
    )
    orthoDF['flag_any_k4_sim_only'] = np.where(
            (orthoDF['flag_any_k4_mel'] == 0)
            & (orthoDF['flag_any_k4_sim'] == 1),
            1,
            0
    )

    # flag_has_male_k4, flag_has_male_k4_mel, flag_has_male_k4_sim
    melDF['flag_has_male_k4'] = np.where(~melDF['gene_k4'].isin(["none", "fem"]), 1, 0)
    simDF['flag_has_male_k4'] = np.where(~simDF['gene_k4'].isin(["none", "fem"]), 1, 0)
    orthoDF['flag_has_male_k4_mel'] = np.where(~orthoDF['gene_k4'].isin(["none", "fem"]), 1, 0)
    orthoDF['flag_has_male_k4_sim'] = np.where(~orthoDF['sim_gene_k4'].isin(["none", "fem"]), 1, 0)
    orthoDF['flag_has_male_k4_both_species'] = np.where(
            (orthoDF['flag_has_male_k4_mel'] == 1)
            & (orthoDF['flag_has_male_k4_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_has_male_k4_mel_only'] = np.where(
            (orthoDF['flag_has_male_k4_mel'] == 1)
            & (orthoDF['flag_has_male_k4_sim'] == 0),
            1,
            0
    )    
    orthoDF['flag_has_male_k4_sim_only'] = np.where(
            (orthoDF['flag_has_male_k4_mel'] == 0)
            & (orthoDF['flag_has_male_k4_sim'] == 1),
            1,
            0
    )

    # flag_has_female_k4, flag_has_female_k4_mel, flag_has_female_k4_sim
    melDF['flag_has_female_k4'] = np.where(~melDF['gene_k4'].isin(["none", "male"]), 1, 0)
    simDF['flag_has_female_k4'] = np.where(~simDF['gene_k4'].isin(["none", "male"]), 1, 0)
    orthoDF['flag_has_female_k4_mel'] = np.where(~orthoDF['gene_k4'].isin(["none", "male"]), 1, 0)
    orthoDF['flag_has_female_k4_sim'] = np.where(~orthoDF['sim_gene_k4'].isin(["none", "male"]), 1, 0)
    orthoDF['flag_has_female_k4_both_species'] = np.where(
            (orthoDF['flag_has_female_k4_mel'] == 1)
            & (orthoDF['flag_has_female_k4_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_has_female_k4_mel_only'] = np.where(
            (orthoDF['flag_has_female_k4_mel'] == 1)
            & (orthoDF['flag_has_female_k4_sim'] == 0),
            1,
            0
    )    
    orthoDF['flag_has_female_k4_sim_only'] = np.where(
            (orthoDF['flag_has_female_k4_mel'] == 0)
            & (orthoDF['flag_has_female_k4_sim'] == 1),
            1,
            0
    )

    # flag_male_limited_k4, flag_male_limited_k4_mel, flag_male_limited_k4_sim
    melDF['flag_male_limited_k4'] = np.where(melDF['gene_k4']=="male", 1, 0)
    simDF['flag_male_limited_k4'] = np.where(simDF['gene_k4']=="male", 1, 0)
    orthoDF['flag_male_limited_k4_mel'] = np.where(orthoDF['gene_k4']=="male", 1, 0)
    orthoDF['flag_male_limited_k4_sim'] = np.where(orthoDF['sim_gene_k4']=="male", 1, 0)
    orthoDF['flag_male_limited_k4_both_species'] = np.where(
            (orthoDF['flag_male_limited_k4_mel'] == 1)
            & (orthoDF['flag_male_limited_k4_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_male_limited_k4_mel_only'] = np.where(
            (orthoDF['flag_male_limited_k4_mel'] == 1)
            & (orthoDF['flag_male_limited_k4_sim'] == 0),
            1,
            0
    )    
    orthoDF['flag_male_limited_k4_sim_only'] = np.where(
            (orthoDF['flag_male_limited_k4_mel'] == 0)
            & (orthoDF['flag_male_limited_k4_sim'] == 1),
            1,
            0
    )

    # flag_female_limited_k4, flag_female_limited_k4_mel, flag_female_limited_k4_sim
    melDF['flag_female_limited_k4'] = np.where(melDF['gene_k4']=="fem", 1, 0)
    simDF['flag_female_limited_k4'] = np.where(simDF['gene_k4']=="fem", 1, 0)
    orthoDF['flag_female_limited_k4_mel'] = np.where(orthoDF['gene_k4']=="fem", 1, 0)
    orthoDF['flag_female_limited_k4_sim'] = np.where(orthoDF['sim_gene_k4']=="fem", 1, 0)
    orthoDF['flag_female_limited_k4_both_species'] = np.where(
            (orthoDF['flag_female_limited_k4_mel'] == 1)
            & (orthoDF['flag_female_limited_k4_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_female_limited_k4_mel_only'] = np.where(
            (orthoDF['flag_female_limited_k4_mel'] == 1)
            & (orthoDF['flag_female_limited_k4_sim'] == 0),
            1,
            0
    )    
    orthoDF['flag_female_limited_k4_sim_only'] = np.where(
            (orthoDF['flag_female_limited_k4_mel'] == 0)
            & (orthoDF['flag_female_limited_k4_sim'] == 1),
            1,
            0
    )

    # flag_sex_limited_k4, flag_sex_limited_k4_mel, flag_sex_limited_k4_sim
    melDF['flag_sex_limited_k4'] = np.where(
        (melDF['gene_k4']=="fem")
        | (melDF['gene_k4']=="male"),
        1,
        0
    )
    simDF['flag_sex_limited_k4'] = np.where(
        (simDF['gene_k4']=="fem")
        | (simDF['gene_k4']=="male"),
        1,
        0
    )
    orthoDF['flag_sex_limited_k4_mel'] = np.where(
        (orthoDF['gene_k4']=="fem")
        | (orthoDF['gene_k4']=="male"),
        1,
        0
    )
    orthoDF['flag_sex_limited_k4_sim'] = np.where(
        (orthoDF['sim_gene_k4']=="fem")
        | (orthoDF['sim_gene_k4']=="male"),
        1,
        0
    )
    orthoDF['flag_sex_limited_k4_both_species'] = np.where(
            (orthoDF['flag_sex_limited_k4_mel'] == 1)
            & (orthoDF['flag_sex_limited_k4_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_sex_limited_k4_mel_only'] = np.where(
            (orthoDF['flag_sex_limited_k4_mel'] == 1)
            & (orthoDF['flag_sex_limited_k4_sim'] == 0),
            1,
            0
    )    
    orthoDF['flag_sex_limited_k4_sim_only'] = np.where(
            (orthoDF['flag_sex_limited_k4_mel'] == 0)
            & (orthoDF['flag_sex_limited_k4_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_conserved_sex_limited_k4'] = np.where(
            (orthoDF['flag_female_limited_k4_both_species'] == 1)
            | (orthoDF['flag_male_limited_k4_both_species'] == 1),
            1,
            0
    )

    # flag_any_k27, flag_any_k27_mel, flag_any_k27_sim
    melDF['flag_any_k27'] = np.where(melDF['gene_k27']!="none", 1, 0)
    simDF['flag_any_k27'] = np.where(simDF['gene_k27']!="none", 1, 0)
    orthoDF['flag_any_k27_mel'] = np.where(orthoDF['gene_k27']!="none", 1, 0)
    orthoDF['flag_any_k27_sim'] = np.where(orthoDF['sim_gene_k27']!="none", 1, 0)
    orthoDF['flag_any_k27_both_species'] = np.where(
            (orthoDF['flag_any_k27_mel'] == 1)
            & (orthoDF['flag_any_k27_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_any_k27_mel_only'] = np.where(
            (orthoDF['flag_any_k27_mel'] == 1)
            & (orthoDF['flag_any_k27_sim'] == 0),
            1,
            0
    )
    orthoDF['flag_any_k27_sim_only'] = np.where(
            (orthoDF['flag_any_k27_mel'] == 0)
            & (orthoDF['flag_any_k27_sim'] == 1),
            1,
            0
    )

    # flag_has_male_k27, flag_has_male_k27_mel, flag_has_male_k27_sim
    melDF['flag_has_male_k27'] = np.where(~melDF['gene_k27'].isin(["none", "fem"]), 1, 0)
    simDF['flag_has_male_k27'] = np.where(~simDF['gene_k27'].isin(["none", "fem"]), 1, 0)
    orthoDF['flag_has_male_k27_mel'] = np.where(~orthoDF['gene_k27'].isin(["none", "fem"]), 1, 0)
    orthoDF['flag_has_male_k27_sim'] = np.where(~orthoDF['sim_gene_k27'].isin(["none", "fem"]), 1, 0)
    orthoDF['flag_has_male_k27_both_species'] = np.where(
            (orthoDF['flag_has_male_k27_mel'] == 1)
            & (orthoDF['flag_has_male_k27_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_has_male_k27_mel_only'] = np.where(
            (orthoDF['flag_has_male_k27_mel'] == 1)
            & (orthoDF['flag_has_male_k27_sim'] == 0),
            1,
            0
    )
    orthoDF['flag_has_male_k27_sim_only'] = np.where(
            (orthoDF['flag_has_male_k27_mel'] == 0)
            & (orthoDF['flag_has_male_k27_sim'] == 1),
            1,
            0
    )

    # flag_has_female_k27, flag_has_female_k27_mel, flag_has_female_k27_sim
    melDF['flag_has_female_k27'] = np.where(~melDF['gene_k27'].isin(["none", "male"]), 1, 0)
    simDF['flag_has_female_k27'] = np.where(~simDF['gene_k27'].isin(["none", "male"]), 1, 0)
    orthoDF['flag_has_female_k27_mel'] = np.where(~orthoDF['gene_k27'].isin(["none", "male"]), 1, 0)
    orthoDF['flag_has_female_k27_sim'] = np.where(~orthoDF['sim_gene_k27'].isin(["none", "male"]), 1, 0)
    orthoDF['flag_has_female_k27_both_species'] = np.where(
            (orthoDF['flag_has_female_k27_mel'] == 1)
            & (orthoDF['flag_has_female_k27_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_has_female_k27_mel_only'] = np.where(
            (orthoDF['flag_has_female_k27_mel'] == 1)
            & (orthoDF['flag_has_female_k27_sim'] == 0),
            1,
            0
    )
    orthoDF['flag_has_female_k27_sim_only'] = np.where(
            (orthoDF['flag_has_female_k27_mel'] == 0)
            & (orthoDF['flag_has_female_k27_sim'] == 1),
            1,
            0
    )

    # flag_male_limited_k27, flag_male_limited_k27_mel, flag_male_limited_k27_sim
    melDF['flag_male_limited_k27'] = np.where(melDF['gene_k27']=="male", 1, 0)
    simDF['flag_male_limited_k27'] = np.where(simDF['gene_k27']=="male", 1, 0)
    orthoDF['flag_male_limited_k27_mel'] = np.where(orthoDF['gene_k27']=="male", 1, 0)
    orthoDF['flag_male_limited_k27_sim'] = np.where(orthoDF['sim_gene_k27']=="male", 1, 0)
    orthoDF['flag_male_limited_k27_both_species'] = np.where(
            (orthoDF['flag_male_limited_k27_mel'] == 1)
            & (orthoDF['flag_male_limited_k27_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_male_limited_k27_mel_only'] = np.where(
            (orthoDF['flag_male_limited_k27_mel'] == 1)
            & (orthoDF['flag_male_limited_k27_sim'] == 0),
            1,
            0
    )
    orthoDF['flag_male_limited_k27_sim_only'] = np.where(
            (orthoDF['flag_male_limited_k27_mel'] == 0)
            & (orthoDF['flag_male_limited_k27_sim'] == 1),
            1,
            0
    )

    # flag_male_limited_k27, flag_male_limited_k27_mel, flag_male_limited_k27_sim
    melDF['flag_female_limited_k27'] = np.where(melDF['gene_k27']=="fem", 1, 0)
    simDF['flag_female_limited_k27'] = np.where(simDF['gene_k27']=="fem", 1, 0)
    orthoDF['flag_female_limited_k27_mel'] = np.where(orthoDF['gene_k27']=="fem", 1, 0)
    orthoDF['flag_female_limited_k27_sim'] = np.where(orthoDF['sim_gene_k27']=="fem", 1, 0)
    orthoDF['flag_female_limited_k27_both_species'] = np.where(
            (orthoDF['flag_female_limited_k27_mel'] == 1)
            & (orthoDF['flag_female_limited_k27_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_female_limited_k27_mel_only'] = np.where(
            (orthoDF['flag_female_limited_k27_mel'] == 1)
            & (orthoDF['flag_female_limited_k27_sim'] == 0),
            1,
            0
    )
    orthoDF['flag_female_limited_k27_sim_only'] = np.where(
            (orthoDF['flag_female_limited_k27_mel'] == 0)
            & (orthoDF['flag_female_limited_k27_sim'] == 1),
            1,
            0
    )

    # flag_sex_limited_k27, flag_sex_limited_k27_mel, flag_sex_limited_k27_sim
    melDF['flag_sex_limited_k27'] = np.where(
        (melDF['gene_k27']=="fem")
        | (melDF['gene_k27']=="male"),
        1,
        0
    )
    simDF['flag_sex_limited_k27'] = np.where(
        (simDF['gene_k27']=="fem")
        | (simDF['gene_k27']=="male"),
        1,
        0
    )
    orthoDF['flag_sex_limited_k27_mel'] = np.where(
        (orthoDF['gene_k27']=="fem")
        | (orthoDF['gene_k27']=="male"),
        1,
        0
    )
    orthoDF['flag_sex_limited_k27_sim'] = np.where(
        (orthoDF['sim_gene_k27']=="fem")
        | (orthoDF['sim_gene_k27']=="male"),
        1,
        0
    )
    orthoDF['flag_sex_limited_k27_both_species'] = np.where(
            (orthoDF['flag_sex_limited_k27_mel'] == 1)
            & (orthoDF['flag_sex_limited_k27_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_sex_limited_k27_mel_only'] = np.where(
            (orthoDF['flag_sex_limited_k27_mel'] == 1)
            & (orthoDF['flag_sex_limited_k27_sim'] == 0),
            1,
            0
    )    
    orthoDF['flag_sex_limited_k27_sim_only'] = np.where(
            (orthoDF['flag_sex_limited_k27_mel'] == 0)
            & (orthoDF['flag_sex_limited_k27_sim'] == 1),
            1,
            0
    )
    orthoDF['flag_conserved_sex_limited_k27'] = np.where(
            (orthoDF['flag_female_limited_k27_both_species'] == 1)
            | (orthoDF['flag_male_limited_k27_both_species'] == 1),
            1,
            0
    )

### RNA flags

    # flag_expressed, flag_expressed_mel, flag_expressed_sim
    melDF['flag_expressed'] = np.where(melDF['gene_rna']!="none", 1, 0)
    simDF['flag_expressed'] = np.where(simDF['gene_rna']!="none", 1, 0)
    orthoDF['flag_expressed_mel'] = np.where(orthoDF['gene_rna']!="none", 1, 0)
    orthoDF['flag_expressed_sim'] = np.where(orthoDF['sim_gene_rna']!="none", 1, 0)

    # flag_expressed_M, flag_expressed_M_mel, flag_expressed_M_sim
    melDF['flag_expressed_M'] = np.where(~melDF['gene_rna'].isin(["none","fem"]), 1, 0)
    simDF['flag_expressed_M'] = np.where(~simDF['gene_rna'].isin(["none","fem"]), 1, 0)
    orthoDF['flag_expressed_M_mel'] = np.where(~orthoDF['gene_rna'].isin(["none","fem"]), 1, 0)
    orthoDF['flag_expressed_M_sim'] = np.where(~orthoDF['sim_gene_rna'].isin(["none","fem"]), 1, 0)

    # flag_expressed_F, flag_expressed_F_mel, flag_expressed_F_sim
    melDF['flag_expressed_F'] = np.where(~melDF['gene_rna'].isin(["none","male"]), 1, 0)
    simDF['flag_expressed_F'] = np.where(~simDF['gene_rna'].isin(["none","male"]), 1, 0)
    orthoDF['flag_expressed_F_mel'] = np.where(~orthoDF['gene_rna'].isin(["none","male"]), 1, 0)
    orthoDF['flag_expressed_F_sim'] = np.where(~orthoDF['sim_gene_rna'].isin(["none","male"]), 1, 0)

    # flag_sex_limited, flag_sex_limited_mel, flag_sex_limited_sim
    melDF['flag_sex_limited'] = np.where((melDF['gene_rna']=="fem")|(melDF['gene_rna']=="male"), 1, 0)
    simDF['flag_sex_limited'] = np.where((simDF['gene_rna']=="fem")|(simDF['gene_rna']=="male"), 1, 0)
    orthoDF['flag_sex_limited_mel'] = np.where((orthoDF['gene_rna']=="fem")|(orthoDF['gene_rna']=="male"), 1, 0)
    orthoDF['flag_sex_limited_sim'] = np.where((orthoDF['sim_gene_rna']=="fem")|(orthoDF['sim_gene_rna']=="male"), 1, 0)

    # flag_male_limited, flag_male_limited_mel, flag_male_limited_sim
    melDF['flag_male_limited'] = np.where((melDF['gene_rna']=="male"), 1, 0)
    simDF['flag_male_limited'] = np.where((simDF['gene_rna']=="male"), 1, 0)
    orthoDF['flag_male_limited_mel'] = np.where((orthoDF['gene_rna']=="male"), 1, 0)
    orthoDF['flag_male_limited_sim'] = np.where((orthoDF['sim_gene_rna']=="male"), 1, 0)

    # flag_female_limited, flag_female_limited_mel, flag_female_limited_sim
    melDF['flag_female_limited'] = np.where((melDF['gene_rna']=="fem"), 1, 0)
    simDF['flag_female_limited'] = np.where((simDF['gene_rna']=="fem"), 1, 0)
    orthoDF['flag_female_limited_mel'] = np.where((orthoDF['gene_rna']=="fem"), 1, 0)
    orthoDF['flag_female_limited_sim'] = np.where((orthoDF['sim_gene_rna']=="fem"), 1, 0)

    # flag_low_express, flag_low_express_mel, flag_low_express_sim
    melDF['flag_low_express'] = np.where(melDF['flag_rna_detected05']==0, 1, 0)
    simDF['flag_low_express'] = np.where(simDF['flag_rna_detected05']==0, 1, 0)
    orthoDF['flag_low_express_mel'] = np.where(orthoDF['flag_rna_detected05']==0, 1, 0)
    orthoDF['flag_low_express_sim'] = np.where(orthoDF['sim_flag_rna_detected05']==0, 1, 0)

    # flag_sex_biased, flag_sex_biased_mel, flag_sex_biased_sim
    melDF['flag_sex_biased'] = np.where((melDF['flag_sex_limited']==0)&(melDF['gene_trend_ttest'].isin(["male", "female", "male_and_female"])), 1, 0)
    simDF['flag_sex_biased'] = np.where((simDF['flag_sex_limited']==0)&(simDF['gene_trend_ttest'].isin(["male", "female", "male_and_female"])), 1, 0)
    orthoDF['flag_sex_biased_mel'] = np.where((orthoDF['flag_sex_limited_mel']==0)&(orthoDF['gene_trend_ttest'].isin(["male", "female", "male_and_female"])), 1, 0)
    orthoDF['flag_sex_biased_sim'] = np.where((orthoDF['flag_sex_limited_sim']==0)&(orthoDF['sim_gene_trend_ttest'].isin(["male", "female", "male_and_female"])), 1, 0)

    # flag_U, flag_U_mel, flag_U_sim
    melDF['flag_U'] = np.where((melDF['flag_sex_limited']==0)&(melDF['gene_trend_ttest']=="unbiased"), 1, 0)
    simDF['flag_U'] = np.where((simDF['flag_sex_limited']==0)&(simDF['gene_trend_ttest']=="unbiased"), 1, 0)
    orthoDF['flag_U_mel'] = np.where((orthoDF['flag_sex_limited_mel']==0)&(orthoDF['gene_trend_ttest']=="unbiased"), 1, 0)
    orthoDF['flag_U_sim'] = np.where((orthoDF['flag_sex_limited_sim']==0)&(orthoDF['sim_gene_trend_ttest']=="unbiased"), 1, 0)

    # flag_M, flag_M_mel, flag_M_sim
    melDF['flag_M'] = np.where((melDF['flag_sex_limited']==0)&(melDF['gene_trend_ttest']=="male"), 1, 0)
    simDF['flag_M'] = np.where((simDF['flag_sex_limited']==0)&(simDF['gene_trend_ttest']=="male"), 1, 0)
    orthoDF['flag_M_mel'] = np.where((orthoDF['flag_sex_limited_mel']==0)&(orthoDF['gene_trend_ttest']=="male"), 1, 0)
    orthoDF['flag_M_sim'] = np.where((orthoDF['flag_sex_limited_sim']==0)&(orthoDF['sim_gene_trend_ttest']=="male"), 1, 0)
    orthoDF['flag_M_mel_only'] = np.where(
            (orthoDF['flag_M_mel']==1)
            & (orthoDF['flag_M_sim']==0),
            1,
            0
    )
    orthoDF['flag_M_sim_only'] = np.where(
            (orthoDF['flag_M_mel']==0)
            & (orthoDF['flag_M_sim']==1),
            1,
            0
    )
    orthoDF['flag_M_mel_U_sim'] = np.where(
            (orthoDF['flag_M_mel']==1)
            & (orthoDF['flag_U_sim']==1),
            1,
            0
    )
    orthoDF['flag_M_sim_U_mel'] = np.where(
            (orthoDF['flag_U_mel']==1)
            & (orthoDF['flag_M_sim']==1),
            1,
            0
    )

    # flag_MP, flag_MP_mel, flag_MP_sim
    melDF['flag_MP'] = np.where((melDF['flag_sex_limited']==0)&(melDF['gene_trend_ttest']=="male")&(melDF['gene_ratio2_ttest']=="male"), 1, 0)
    simDF['flag_MP'] = np.where((simDF['flag_sex_limited']==0)&(simDF['gene_trend_ttest']=="male")&(simDF['gene_ratio2_ttest']=="male"), 1, 0)
    orthoDF['flag_MP_mel'] = np.where((orthoDF['flag_sex_limited_mel']==0)&(orthoDF['gene_trend_ttest']=="male")&(orthoDF['gene_ratio2_ttest']=="male"), 1, 0)
    orthoDF['flag_MP_sim'] = np.where((orthoDF['flag_sex_limited_sim']==0)&(orthoDF['sim_gene_trend_ttest']=="male")&(orthoDF['sim_gene_ratio2_ttest']=="male"), 1, 0)

    # flag_F, flag_F_mel, flag_F_sim
    melDF['flag_F'] = np.where((melDF['flag_sex_limited']==0)&(melDF['gene_trend_ttest']=="female"), 1, 0)
    simDF['flag_F'] = np.where((simDF['flag_sex_limited']==0)&(simDF['gene_trend_ttest']=="female"), 1, 0)
    orthoDF['flag_F_mel'] = np.where((orthoDF['flag_sex_limited_mel']==0)&(orthoDF['gene_trend_ttest']=="female"), 1, 0)
    orthoDF['flag_F_sim'] = np.where((orthoDF['flag_sex_limited_sim']==0)&(orthoDF['sim_gene_trend_ttest']=="female"), 1, 0)
    orthoDF['flag_F_mel_only'] = np.where(
            (orthoDF['flag_F_mel']==1)
            & (orthoDF['flag_F_sim']==0),
            1,
            0
    )
    orthoDF['flag_F_sim_only'] = np.where(
            (orthoDF['flag_F_mel']==0)
            & (orthoDF['flag_F_sim']==1),
            1,
            0
    )
    orthoDF['flag_F_mel_U_sim'] = np.where(
            (orthoDF['flag_F_mel']==1)
            & (orthoDF['flag_U_sim']==1),
            1,
            0
    )
    orthoDF['flag_F_sim_U_mel'] = np.where(
            (orthoDF['flag_U_mel']==1)
            & (orthoDF['flag_F_sim']==1),
            1,
            0
    )
    orthoDF['flag_M_mel_F_sim'] = np.where(
            (orthoDF['flag_M_mel']==1)
            & (orthoDF['flag_F_sim']==1),
            1,
            0
    )
    orthoDF['flag_M_sim_F_mel'] = np.where(
            (orthoDF['flag_F_mel']==1)
            & (orthoDF['flag_M_sim']==1),
            1,
            0
    )

    # flag_FP, flag_FP_mel, flag_FP_sim
    melDF['flag_FP'] = np.where((melDF['flag_sex_limited']==0)&(melDF['gene_trend_ttest']=="female")&(melDF['gene_ratio2_ttest']=="female"), 1, 0)
    simDF['flag_FP'] = np.where((simDF['flag_sex_limited']==0)&(simDF['gene_trend_ttest']=="female")&(simDF['gene_ratio2_ttest']=="female"), 1, 0)
    orthoDF['flag_FP_mel'] = np.where((orthoDF['flag_sex_limited_mel']==0)&(orthoDF['gene_trend_ttest']=="female")&(orthoDF['gene_ratio2_ttest']=="female"), 1, 0)
    orthoDF['flag_FP_sim'] = np.where((orthoDF['flag_sex_limited_sim']==0)&(orthoDF['sim_gene_trend_ttest']=="female")&(orthoDF['sim_gene_ratio2_ttest']=="female"), 1, 0)

    # flag_MF, flag_MF_mel, flag_MF_sim
    melDF['flag_MF'] = np.where((melDF['flag_sex_limited']==0)&(melDF['gene_trend_ttest']=="male_and_female"), 1, 0)
    simDF['flag_MF'] = np.where((simDF['flag_sex_limited']==0)&(simDF['gene_trend_ttest']=="male_and_female"), 1, 0)
    orthoDF['flag_MF_mel'] = np.where((orthoDF['flag_sex_limited_mel']==0)&(orthoDF['gene_trend_ttest']=="male_and_female"), 1, 0)
    orthoDF['flag_MF_sim'] = np.where((orthoDF['flag_sex_limited_sim']==0)&(orthoDF['sim_gene_trend_ttest']=="male_and_female"), 1, 0)
    orthoDF['flag_MF_mel_only'] = np.where(
            (orthoDF['flag_MF_mel']==1)
            & (orthoDF['flag_MF_sim']==0),
            1,
            0
    )
    orthoDF['flag_MF_sim_only'] = np.where(
            (orthoDF['flag_MF_mel']==0)
            & (orthoDF['flag_MF_sim']==1),
            1,
            0
    )
    orthoDF['flag_MF_mel_U_sim'] = np.where(
            (orthoDF['flag_MF_mel']==1)
            & (orthoDF['flag_U_sim']==1),
            1,
            0
    )
    orthoDF['flag_MF_sim_U_mel'] = np.where(
            (orthoDF['flag_U_mel']==1)
            & (orthoDF['flag_MF_sim']==1),
            1,
            0
    )

    # flag_MFP, flag_MFP_mel, flag_MFP_sim
    melDF['flag_MFP'] = np.where((melDF['flag_sex_limited']==0)&(melDF['gene_trend_ttest']=="male_and_female")&(melDF['gene_ratio2_ttest']=="male_and_female"), 1, 0)
    simDF['flag_MFP'] = np.where((simDF['flag_sex_limited']==0)&(simDF['gene_trend_ttest']=="male_and_female")&(simDF['gene_ratio2_ttest']=="male_and_female"), 1, 0)
    orthoDF['flag_MFP_mel'] = np.where((orthoDF['flag_sex_limited_mel']==0)&(orthoDF['gene_trend_ttest']=="male_and_female")&(orthoDF['gene_ratio2_ttest']=="male_and_female"), 1, 0)
    orthoDF['flag_MFP_sim'] = np.where((orthoDF['flag_sex_limited_sim']==0)&(orthoDF['sim_gene_trend_ttest']=="male_and_female")&(orthoDF['sim_gene_ratio2_ttest']=="male_and_female"), 1, 0)

    # flag_sex_bias_both_species
    orthoDF["flag_sex_bias_both_species"] = np.where(
        (orthoDF["xsome"]==orthoDF["sim_xsome"])
        & (orthoDF["flag_sex_biased_mel"] == 1)
        & (orthoDF["flag_sex_biased_sim"] == 1),
        1,
        0
    )
    # flag_conserved_sex_bias, flag_conserved_male_bias, flag_conserved_female_bias
    orthoDF["flag_conserved_sex_bias"] = np.where(
        (orthoDF["xsome"]==orthoDF["sim_xsome"])
        & (orthoDF["flag_sex_limited_mel"] + orthoDF["flag_sex_limited_sim"] == 0)
        & (orthoDF["gene_trend_ttest"]==orthoDF["sim_gene_trend_ttest"])
        & (orthoDF["flag_sex_bias_both_species"] == 1),
        1,
        0
    )
    orthoDF["flag_conserved_male_bias"] = np.where(
        (orthoDF["flag_conserved_sex_bias"]==1)
        & (orthoDF["gene_trend_ttest"]=="male"),
        1,
        0
    )
    orthoDF["flag_conserved_female_bias"] = np.where(
        (orthoDF["flag_conserved_sex_bias"]==1)
        & (orthoDF["gene_trend_ttest"]=="female"),
        1,
        0
    )
    orthoDF["flag_conserved_male_and_female"] = np.where(
        (orthoDF["flag_MF_mel"]==1)
        & (orthoDF["flag_MF_sim"]==1),
        1,
        0
    )    
    orthoDF["flag_conserved_unbiased"] = np.where(
        (orthoDF["flag_U_mel"]==1)
        & (orthoDF["flag_U_sim"]==1),
        1,
        0
    )
    # flag_sex_bias_mel_only, flag_sex_bias_sim_only
    orthoDF["flag_sex_bias_mel_only"] = np.where(
        (orthoDF["xsome"]==orthoDF["sim_xsome"])
        & (orthoDF["flag_sex_biased_mel"] == 1)
        & (orthoDF["flag_sex_biased_sim"] == 0),
        1,
        0
    )
    orthoDF["flag_sex_bias_sim_only"] = np.where(
        (orthoDF["xsome"]==orthoDF["sim_xsome"])
        & (orthoDF["flag_sex_biased_mel"] == 0)
        & (orthoDF["flag_sex_biased_sim"] == 1),
        1,
        0
    )

    # Make all divergence flags
    # flag_switch_M_U, flag_switch_F_U, flag_switch_M_F, flag_diverged_sex_bias
    orthoDF["flag_switch_M_U"] = np.where(
        (orthoDF["xsome"]==orthoDF["sim_xsome"])
        & (orthoDF["flag_M_mel_U_sim"] + orthoDF["flag_M_sim_U_mel"] > 0),
        1,
        0
    )
    orthoDF["flag_switch_F_U"] = np.where(
        (orthoDF["xsome"]==orthoDF["sim_xsome"])
        & (orthoDF["flag_F_mel_U_sim"] + orthoDF["flag_F_sim_U_mel"] > 0),
        1,
        0
    )
    orthoDF["flag_switch_M_F"] = np.where(
        (orthoDF["xsome"]==orthoDF["sim_xsome"])
        & (orthoDF["flag_M_mel_F_sim"] + orthoDF["flag_M_sim_F_mel"] > 0),
        1,
        0
    )
    orthoDF["flag_diverged_sex_bias"] = np.where(
        (orthoDF["xsome"]==orthoDF["sim_xsome"])
        & (orthoDF["flag_sex_biased_mel"] + orthoDF["flag_sex_biased_sim"] > 0)
        & (orthoDF["flag_conserved_sex_bias"] == 0),
        1,
        0
    )

    # flag_one2one_ortholog
    orthoDF["num_mel_gene_ID_2_sim_geneID"] = orthoDF.groupby("sim_geneID")["mel_geneID"].transform("nunique")
    orthoDF["num_sim_gene_ID_2_mel_geneID"] = orthoDF.groupby("mel_geneID")["sim_geneID"].transform("nunique")
    orthoDF["flag_one2one_ortholog"] = np.where(
            (orthoDF["num_mel_gene_ID_2_sim_geneID"]==1)
            & (orthoDF["num_sim_gene_ID_2_mel_geneID"]==1),
            1,
            0
    )

    # Drop unnecessary flags
    melDF = melDF.drop(columns={
        "gene_k4_sex", "gene_k4_bias", "gene_k27_sex", "gene_k27_bias",
        "sum_fragments_rna_detected05", "unb", "num_ratio2_male", "num_ratio2_fem",
        "sum_fragments_ratio2", "gene_ratio2", "num_ratio_male", "num_ratio_fem",
        "sum_fragments_ratio_trend", "gene_ratio_trend", "unbias", "male_exp",
        "fem_exp", "male_full", "male_swch", "fem_full", "fem_swch",
        "sum_male_express", "sum_female_express", "num_fragments_express",
        "gene_sex_bias", "flag_fem_full", "flag_fem_exp", "flag_male_full",
        "flag_male_exp", "flag_male_swch", "flag_fem_swch", "flag_unb_k4e",
        "organism_abbreviation", "secondary_FBgn_s_", "annotation_ID",
        "secondary_annotation_ID_s_", "PERCENT", "XCHR", "gene_id_x", "gene_id_y",
        "flag_graze2014_sim_m_bias", "flag_graze2014_sim_f_bias"
    })
    simDF = simDF.drop(columns={
        "gene_k4_sex", "gene_k4_bias", "gene_k27_sex", "gene_k27_bias",
        "sum_fragments_rna_detected05", "unb", "num_ratio2_male", "num_ratio2_fem",
        "sum_fragments_ratio2", "gene_ratio2", "num_ratio_male", "num_ratio_fem",
        "sum_fragments_ratio_trend", "gene_ratio_trend", "male_exp",
        "fem_exp", "male_full", "male_swch", "fem_full", "fem_swch",
        "sum_male_express", "sum_female_express", "num_fragments_express",
        "gene_sex_bias", "flag_fem_full", "flag_fem_exp", "flag_male_full",
        "flag_male_exp", "flag_male_swch", "flag_fem_swch", "flag_unb_k4e", 
        "organism_abbreviation", "secondary_FBgn_s_", "annotation_ID",
        "secondary_annotation_ID_s_", "PERCENT", "XCHR"
    })
    orthoDF = orthoDF.drop(columns=[
        "gene_id_x", "gene_id_y", "gene_k4_sex", "gene_k4_bias", "gene_k27_sex",
        "gene_k27_bias", "sim_gene_k4_sex", "sim_gene_k4_bias", "sim_gene_k27_sex",
        "gene_k27_bias", "XCHR", "PERCENT", "sim_fbgn", "symbol", "sim_symbol",
        "sum_fragments_rna_detected05", "sim_sum_fragments_rna_detected05", "unb",
        "num_ratio2_male", "num_ratio2_fem", "sum_fragments_ratio2", "gene_ratio2",
        "sim_num_ratio2_male", "sim_num_ratio2_fem", "sim_sum_fragments_ratio2",
        "sim_gene_ratio2", "num_ratio_male", "num_ratio_fem", "sum_fragments_ratio_trend",
        "gene_ratio_trend", "sim_num_ratio_male", "sim_num_ratio_fem", "sim_sum_fragments_ratio_trend",
        "sim_gene_ratio_trend", "unbias", "male_exp", "fem_exp", "male_full",
        "male_swch", "fem_full", "fem_swch", "sum_male_express", "sum_female_express",
        "num_fragments_express", "gene_sex_bias", "flag_fem_full", "flag_fem_exp",
        "flag_male_full", "flag_male_exp", "flag_male_swch", "flag_fem_swch",
        "flag_unb_k4e", "sim_male_exp", "sim_fem_exp", "sim_male_full", "sim_male_swch",
        "sim_fem_full", "sim_fem_swch", "sim_sum_male_express", "sim_sum_female_express",
        "sim_num_fragments_express", "sim_gene_sex_bias", "sim_flag_fem_full",
        "sim_flag_fem_exp", "sim_flag_male_full", "sim_flag_male_exp", "sim_flag_male_swch",
        "sim_flag_fem_swch", "sim_flag_unb_k4e", "chrom", "start", "end", "strand",
        "organism_abbreviation", "secondary_FBgn_s_", "annotation_ID", "secondary_annotation_ID_s_",
        "species_ratio2_ttest", "species_ratio2_ttest_sex", "species_trend_ttest",
        "species_trend_ttest_sex", "flag_graze2014_sim_m_bias", "flag_graze2014_sim_f_bias",
        
    ])

    # Output files to designated output directory
    melDF.to_csv(args.outDir+"/dmel_chip_rna_flags.csv", index=False)
    simDF.to_csv(args.outDir+"/dsim_chip_rna_flags.csv", index=False)
    orthoDF.to_csv(args.outDir+"/dmel_dsim_ortholog_chip_rna_flags.csv", index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
