#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
            description="Evaluate transcript lengths and BLAST results."
    )

    # Input data
    parser.add_argument(
        "-l",
        "--length",
        dest="inLen",
        required=True,
        help="Input sequence length file (query of BLAST)."
    )
    parser.add_argument(
        "-b",
        "--blast",
        dest="inBlast",
        required=True,
        help="BLAST CSV output file."
    )
    parser.add_argument(
        "-f",
        "--flag",
        dest="inFlag",
        required=False,
        type=str,
        action="append",
        help=(
            "String to search for, flag, and count in subject names of BLAST "
            "results. Use a new argument for each string to flag."
        )
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output-prefix",
        dest="outPrefix",
        required=True,
        help="Output prefix."
    )

    args = parser.parse_args()
    return args

def main():
    # Get intput files
    lengthDF = pd.read_csv(args.inLen, low_memory=False)
    blastDF = pd.read_csv(args.inBlast, sep="\t", low_memory=False)

    # Check how many query sequences in the blast output have multiple hits to the same subject
#    blastQSub = blastDF.groupby(
#            ["qseqid", "stitle"]
#        )["qstart"].count().reset_index().rename(
#                columns={"qstart": "num_q_hit_to_s"}
#    )
#    blastQSub["num_q_hit_to_s"].value_counts()

    # Flag blast hits that are full matches (where length == qlen)
    blastDF["flag_full_match"] = np.where(
        blastDF["length"] == blastDF["qlen"],
        1,
        0
    )

    # Flag blast hits that are at least 95%, 90%, 85%, 80%, and 75% the query length
    blastDF["perc_qlen"] = blastDF["length"] / blastDF["qlen"]
    blastDF["flag_95_perc_qlen"] = np.where(
        blastDF["perc_qlen"] >= 0.95,
        1,
        0
    )
    blastDF["flag_90_perc_qlen"] = np.where(
        blastDF["perc_qlen"] >= 0.90,
        1,
        0
    )
    blastDF["flag_85_perc_qlen"] = np.where(
        blastDF["perc_qlen"] >= 0.85,
        1,
        0
    )
    blastDF["flag_80_perc_qlen"] = np.where(
        blastDF["perc_qlen"] >= 0.80,
        1,
        0
    )
    blastDF["flag_75_perc_qlen"] = np.where(
        blastDF["perc_qlen"] >= 0.75,
        1,
        0
    )

    # Open log file
    logFile = open(args.outPrefix + ".log", "w")

    # Add flags for strings in blast results if requested
    if args.inFlag is not None:
        if type(args.inFlag) == str:
            flagList = list()
            flagList.append(args.inFlag)
        else:
            flagList = args.inFlag
        logFile.write("Counts for flags in full BLAST output file:\n\n")
        for flag in flagList:
            blastDF["flag_hit_to_" + flag.replace(" ", "_")] = np.where(
                blastDF["stitle"].str.contains(flag),
                1,
                0
            )
            logFile.write("{}\n{}\n\n".format(
                    "flag_hit_to_" + flag.replace(" ", "_"),
                    blastDF.groupby("qseqid")[
                            "flag_hit_to_" + flag.replace(" ", "_")
                        ].max().reset_index()[
                                "flag_hit_to_" + flag.replace(" ", "_")
                            ].value_counts(sort=False).to_string()
            ))


    # Output number of transcripts < 500, >= 500nt and <1000nt, and >=1000nt
    ranges = [0, 500, 1000, lengthDF["length"].max()+1]
    counts = lengthDF.groupby(pd.cut(lengthDF["length"], ranges, right=False))["length"].count()
    logFile.write("Read Length Counts:\n\n{}\n\n".format(counts.to_string()))

    # Begin output CSV with seq length counts
    outputDF = pd.DataFrame(
        [[
            len(lengthDF),
            len(lengthDF[lengthDF["length"]<500]),
            len(lengthDF[(lengthDF["length"]>=500) & (lengthDF["length"]<1000)]),
            len(lengthDF[lengthDF["length"]>=1000])
        ]],
        columns=[
            "num_seq",
            "num_seq_lt_500nt",
            "num_seq_ge_500_lt_1000",
            "num_seq_ge_1000",
    ])

    # Of the transcripts >=1000nt, describe the blast results
    selectLenDF = lengthDF[lengthDF["length"]>=1000].copy()

    # Subset blast results only for queries in the selected sequences
    selectBlast = blastDF[blastDF["qseqid"].isin(selectLenDF["name"])].copy()
    outputDF["num_seq_ge_1000_w_blast_hit"] = selectBlast["qseqid"].nunique()

    # Count flags for strings in selected blast results if requested
    if args.inFlag is not None:
        logFile.write("Counts for flags in selected BLAST output file of reads >=1000nt:\n\n")
        for flag in flagList:
            tempFlag = selectBlast.groupby("qseqid")[
                    "flag_hit_to_" + flag.replace(" ", "_")
                ].max().reset_index()
            try:
                outputDF["num_seq_ge_1000_hit_to_" + flag.replace(" ", "_")] = tempFlag[
                        "flag_hit_to_" + flag.replace(" ", "_")].value_counts()[1]
            except:
                outputDF["num_seq_ge_1000_hit_to_" + flag.replace(" ", "_")] = 0
            logFile.write("{}\n{}\n\n".format(
                    "flag_hit_to_" + flag.replace(" ", "_"),
                    tempFlag["flag_hit_to_" + flag.replace(" ", "_")].value_counts(sort=False).to_string()
            ))

    # Sort the selected blast results and select the top hit for each query sequence
    selectBlast = selectBlast.sort_values(["bitscore", "evalue"], ascending=[False, True])
    selectBestHit = selectBlast.groupby("qseqid").first().reset_index()
    # Sort best hits by bitscore and evalue
    selectBestHit = selectBestHit.sort_values(["bitscore", "evalue"], ascending=[False, True])

    # Add to output the counts of best hits with perc_qlen flags
    outputDF["num_seq_ge_1000_best_hit_95_perc_qlen"] = selectBestHit["flag_95_perc_qlen"].value_counts()[1]
    outputDF["num_seq_ge_1000_best_hit_90_perc_qlen"] = selectBestHit["flag_90_perc_qlen"].value_counts()[1]
    outputDF["num_seq_ge_1000_best_hit_85_perc_qlen"] = selectBestHit["flag_85_perc_qlen"].value_counts()[1]
    outputDF["num_seq_ge_1000_best_hit_80_perc_qlen"] = selectBestHit["flag_80_perc_qlen"].value_counts()[1]
    outputDF["num_seq_ge_1000_best_hit_75_perc_qlen"] = selectBestHit["flag_75_perc_qlen"].value_counts()[1]
    outputDF["mean_perc_qlen_seq_ge_1000_best_hit"] = selectBestHit["perc_qlen"].mean()
    outputDF["median_perc_qlen_seq_ge_1000_best_hit"] = selectBestHit["perc_qlen"].median()

    # Count flags for strings in selected best hit blast results if requested
    if args.inFlag is not None:
        logFile.write("Counts for flags in selected BLAST output file of reads >=1000nt:\n\n")
        for flag in flagList:
            selectBestHit["flag_seq_ge_1000_best_hit_to_" + flag.replace(" ", "_")] = np.where(
                selectBestHit["stitle"].str.contains(flag),
                1,
                0
            )
            try:
                outputDF["num_seq_ge_1000_best_hit_to_" + flag.replace(" ", "_")] = selectBestHit[
                        "flag_seq_ge_1000_best_hit_to_" + flag.replace(" ", "_")].value_counts()[1]
            except:
                outputDF["num_seq_ge_1000_best_hit_to_" + flag.replace(" ", "_")] = 0
            logFile.write("{}\n{}\n\n".format(
                    "flag_seq_ge_1000_best_hit_to_" + flag.replace(" ", "_"),
                    selectBestHit["flag_seq_ge_1000_best_hit_to_" + flag.replace(" ", "_")].value_counts(sort=False).to_string()
            ))

    # Output counts and subset of blast results to selected sequences
    outputDF.to_csv(args.outPrefix + "_counts.csv", index = False)
    selectBlast.to_csv(args.outPrefix + "_seq_ge_1000_blast.csv", index=False)
    selectBestHit.to_csv(args.outPrefix + "_seq_ge_1000_best_hit.csv", index=False)

    logFile.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
