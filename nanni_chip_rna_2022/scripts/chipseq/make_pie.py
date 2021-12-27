#!/usr/bin/env python3

import matplotlib

## Fix for error on HPG with _tkinter
## Solution found on https://help.rc.ufl.edu/doc/Python3
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import argparse
import sys

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Make pie chart")

    # Input data
    parser.add_argument("-l", "--labels", dest="inLabels", required=True, help="Comma separated list (no spaces) of labels for pie chart, must be in the same order as the values")
    parser.add_argument("-v", "--values", dest="inValues", required=True, help="Comma separated list (no spaces) of values to be used in the pie chart, must be in the same order as the labels")
    parser.add_argument("-c", "--colors", dest="inColors", required=False, help="Comma separated list (no spaces) of colors to use in the pie chart, be in the same order as the labels and values. List of valid color names on https://matplotlib.org/examples/color/named_colors.html", default=None)

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output file prefix for pie chart")
    parser.add_argument("-f", "--format", dest="outFormat", required=False, help="Output format for output file", default="png", choices=["png","tiff","eps"])

    args = parser.parse_args()
    return args
def main():
    # Split arguments
    names = args.inLabels.split(",")
    nums = list(map(int,args.inValues.split(",")))
    if args.inColors != None :
        colors = args.inColors.split(",")
    else :
        colors = None

    # Check number of arguments all match
    if len(names) != len(nums) :
        print("ERROR : Number of labels does not match number of values provided.")
        sys.exit()
    if colors != None and len(colors) != len(names) :
        print("ERROR : Number of colors does not match number of labels/values provided.")
        sys.exit()
        
    # Create pie chart
    plt.pie(nums, labels=names, colors=colors, autopct=lambda p: '{:.2f}% ({:.0f})'.format(p,(p/100)*sum(nums)) )
    plt.axis("equal")
#    plt.title(args.outFile.split("/")[len(args.outFile.split("/"))-1])
    plt.savefig(args.outFile+"."+args.outFormat, format=args.outFormat)
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
