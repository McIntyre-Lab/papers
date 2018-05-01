#!/usr/bin/env python

#	DESCRIPTION: This program builds the database file required by GFFutils for GFF manipulation

# Built-in packages
import argparse
import os

# Add-on packages
import gffutils

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Provide the path to the GFF file for creating the database")
    parser.add_argument("--gff", dest="gffInput", action='store', required=True, help="Input GFF file")
    
    args = parser.parse_args()
    return args

def main():
    gff_fn=args.gffInput
    db_fn=gff_fn + '.db'

    gffutils.create_db(gff_fn, db_fn,merge_strategy='create_unique')
    #gffutils.create_db(gff_fn, db_fn)

                        
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

