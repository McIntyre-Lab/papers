#!/usr/bin/env python3

import argparse
import sys
import os
import subprocess

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in values to create a venn diagram")

    # Input data
    parser.add_argument("-v", "--values", dest="values", required=True, help="Comma separated list (no spaces) of values for venn diagram in the order of (10, 01, 11) for 2 categories (100, 010, 110, 001, 101, 011, 111) for 3 categories and (1000, 0100, 1100, 0010, 1010, 0110, 1110, 0001, 1001, 0101, 1101, 0011, 1011, 0111, 1111) for 4 categories")
    parser.add_argument("-l", "--labels", dest="labels", default="A,B,C", help="Comma separated list (no spaces) of labels for the categories (default A,B,C)")
    parser.add_argument("-p", "--scriptPath", dest="path", default=".", help="Path to scripts directory that contains venn_diagram.R script")

    # Output data
    parser.add_argument("-o", "--output-prefix", dest="outFile", required=True, help="Output file prefix (including path)")

    args = parser.parse_args()
    return args

def main():
    # Split value and label arguments
    nums = args.values.split(",")
    names = args.labels.split(",")

    # Check that label and value numbers match
    if(len(nums)!=pow(2,len(names))-1):
        print("ERROR : Incorrect number of values/labels")
        sys.exit()

    # Check for R script
    if(not os.path.exists(args.path+"/venn_diagram.R")):
        print("ERROR : Invalid path to venn_diagram.R")
        sys.exit()
    else:
        script = args.path+"/venn_diagram.R"

    # Create venn diagram using R script
    subprocess.call(['Rscript',script,",".join(names),",".join(nums),args.outFile])


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()