#!/bin/bash

PROJ=/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/sex_specific_splicing


# Check if sufficient arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <column_name>"
    exit 1
fi

# Assign arguments to variables
input_file="$1"
column_name="$2"

# Find the column number of the specified column
col_num=$(head -1 "$input_file" | tr ',' '\n' | nl -v 0 | grep -w "$column_name" | awk '{print $1 + 1}')

# Check if the column was found
if [ -z "$col_num" ]; then
    echo "Column '$column_name' not found in file '$input_file'."
    exit 1
fi

# Find the widest row in the specified column
widest_row=$(awk -v col="$col_num" -F',' '
NR > 1 { if (length($col) > max_width) { max_width = length($col); row = $0 } }
END { print row, max_width }
' "$input_file")

# Output the results
if [ -n "$widest_row" ]; then
    echo "Input File: $input_file"
    echo "Widest Row Width in '$column_name' column: $(echo "$widest_row" | awk '{print $NF}')"
else
    echo "No data found in '$column_name' column in file '$input_file'."
fi


