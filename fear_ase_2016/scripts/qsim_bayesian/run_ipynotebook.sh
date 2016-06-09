#!/bin/bash
set -o nounset
set -e 

# Run the IPython Notebook
cd $(dirname $0)
ipython notebook
