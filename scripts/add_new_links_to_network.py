#!/usr/bin/env python
# Built-in packages
import argparse
from argparse import RawDescriptionHelpFormatter
import logging
import sys

# Add-on packages
import numpy as np

# McLab Packages
import mclib_Python as mclib
import semnet

def getOptions():
    """ Function to pull in arguments """
    description = """"This script starts with a begining pathway in the format
    of a 'path' file. It then builds SAS code for SEMs after adding new links
    in between genes already in the network. New links will include:

        (a) Exogenous to endogenous genes
        (b) Exogenous to exogenous genes
        (c) Endogenous to endogenous genes
        (d) Endogenous to exogenous genes

    Current links will be kept.
    """

    parser = argparse.ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-p", dest="pname", action='store', required=True, help="Name of the 'path' file [Required]")
    parser.add_argument("-l", dest="lname", action='store', required=True, help="Path to sas library that has the SAS dataset that will be analyzed [Required]")
    parser.add_argument("-m", dest="mname", action='store', required=True, help="Name of SAS dataset that will be analyzed [Required]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name with PATH of the output file ending with '.sas'; NOTE: model number will be appended to the filename [Required]")
    parser.add_argument("-t", dest="template", action='store', required=False, help="Name of the PROC CALIS template file [Optional]")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Name of the log file [Optional]; NOTE: if no file is provided logging information will be output to STDOUT") 
    args = parser.parse_args()
    return(args)

if __name__ == '__main__':
    args = getOptions()

    # Turn on logging
    logger = logging.getLogger()
    mclib.logger.setLogger(logger, args.log)

    # Initialize base variable list
    path = semnet.createPath(args.pname)

    # Run baseline model
    args.gname=None
    semnet.models_baseline.baseline(path, args)

    ## Add new links for gene existing in the network
    semnet.models_newlink.add_newlinks_beta_to_beta(path, args)
    semnet.models_newlink.add_newlinks_gamma_to_beta(path, args)
    semnet.models_newlink.add_newlinks_beta_to_gamma(path, args)
    semnet.models_newlink.add_newlinks_gamma_to_gamma(path, args)
