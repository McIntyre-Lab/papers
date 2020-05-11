#!/usr/bin/env python2.7
# Daniel B. Stribling
# 2017-07
# Searches a VCF file for the header line, and renames the first sample
#   column.

import os
import sys
import argparse
import vcf

SCRIPT_DESCRIPTION = 'Renames the sample column of a VCF file containing a single sample.'

def rename_sample(vcf_file, output, new_sample_name, verbose=False, **kwargs):
    if verbose:
        print 'Reading info from VCF file: %s' % vcf_file
        print 'Outputting to file: %s' % output
        
    with open(vcf_file, 'r') as vcf_file_obj, open(output, 'w') as output_obj:
        for line in vcf_file_obj:
            #If this is the header line, operate on it.
            if line.startswith('#') and not line.startswith('##'):
                split_line = line.split()
                if len(split_line) != 10:
                    message = ('Cannot Correctly Parse Header Line:\n' + line.rstrip()
                               + '\nExiting...\n')
                    raise NotImplementedError(message)           
 
                if verbose:
                    print 'Renaming Sample: %s to %s' % (split_line[9], new_sample_name)
                split_line[9] = new_sample_name
                new_line = '\t'.join(split_line)
                output_obj.write(new_line + '\n')     
            #Otherwise, write the line as is.
            else:
                output_obj.write(line)

        if verbose:
            print 'VCF Writing Complete.\n'


# Arg Parsing Methods
def dir_exists(dir_name, dir_type=''):
    if '~' in dir_name:
        dir_name = os.path.expanduser(dir_name)
    dir_name = os.path.abspath(dir_name)
    if not os.path.isdir(dir_name):
        message = ('Provided %sDirectory: %s Does not exist.\n' % (dir_type.title()+' ', dir_name)
                   + 'Exiting...\n')
        raise argparse.ArgumentTypeError(message)
    return dir_name

def input_file_exists(file_name):
    dir_exists(os.path.dirname(file_name), 'input')
    if not os.path.isfile(file_name):
        message = 'Provided Input File: %s Does not exist.\nExiting...\n' % file_name
        raise argparse.ArgumentTypeError(message)
    return os.path.abspath(file_name)

def output_file_exists(file_name):
    YES_ANSWERS = ['y', 'yes', 'go', 'continue','absolutely']
    YES_ANSWERS += ([answer.title() for answer in YES_ANSWERS] 
                    + [answer.upper() for answer in YES_ANSWERS])
    dir_exists(os.path.dirname(file_name), 'output')
    if os.path.isfile(file_name):
        response = raw_input('Output File: %s Already Exists, Overwirte? : ' % file_name).strip()
        if response not in YES_ANSWERS:
            message = 'Overwrite Not Confirmed.\nExiting...\n'
            raise argparse.ArgumentTypeError(message)
    return os.path.abspath(file_name)


def parse_args_optional(print_help=False):
    parser = argparse.ArgumentParser(description=SCRIPT_DESCRIPTION)
    parser.add_argument('vcf_file', type=input_file_exists,
                        help='Input VCF File')
    parser.add_argument('output',
                        help='Output File')
    parser.add_argument('new_sample_name',
                        help='New Name for Sample Column')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print Verbose Output')

    if print_help:
        parser.print_help()
        return {}

    namespace = parser.parse_args()
    args_info = vars(namespace)
    return args_info


# Call Main Function
if __name__ == '__main__':
    args_info = parse_args_optional()
    if not args_info:
        sys.exit()
    rename_sample(**args_info)




