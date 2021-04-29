#!/usr/bin/env python3

import matplotlib

## Fix for error on HPG with _tkinter
## Solution found on https://help.rc.ufl.edu/doc/Python3
matplotlib.use('Agg')

from matplotlib_venn import venn3, venn2
import argparse
import sys

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in values to create a venn diagram")

    # Input data
    parser.add_argument("-n", "--num-categories", dest="num", default=2, help="Number of categories in venn diagram (can be 2 or 3, default is 2)")
    parser.add_argument("-v", "--values", dest="values", required=True, help="Comma separated list (no spaces) of values for venn diagram in the order of Ab,aB,AB for 2 categories or Abc,aBc,ABc,abC,AbC,aBC,ABC for 3 categories")
    parser.add_argument("-l", "--labels", dest="labels", default="A,B,C", help="Comma separated list (no spaces) of labels for the categories (default A,B,C)")

    # Output data
    parser.add_argument("-o", "--output-prefix", dest="outFile", required=True, help="Output file name (should have .png ending)")

    args = parser.parse_args()
    return args
def main():
    # Split value and label arguments
    nums = args.values.split(",")
    names = args.labels.split(",")

    # Create venn diagram
    if int(args.num) == 2 :
        if len(nums) != 3 or len(names) < 2 : # need < 2 because the default of labels will always be 3
            print("ERROR : Incorrect number of values/labels")
            sys.exit()
        venn2(subsets=(int(nums[0]),int(nums[1]),int(nums[2])), set_labels=(names[0],names[1]))
    elif int(args.num) == 3 :
        if len(nums) != 7 or len(names) != 3 :
            print("ERROR : Incorrect number of values/labels")
            sys.exit()
        venn3((int(nums[0]),int(nums[1]),int(nums[2]),int(nums[3]),int(nums[4]),int(nums[5]),int(nums[6])), set_labels = (names[0],names[1],names[2]))
    else :
        print("ERROR : Number of categories must be 2 or 3")
        sys.exit()

    # Print venn diagram to output png file
    matplotlib.pyplot.title(args.outFile.split("/")[len(args.outFile.split("/"))-1])
    matplotlib.pyplot.savefig(args.outFile)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
