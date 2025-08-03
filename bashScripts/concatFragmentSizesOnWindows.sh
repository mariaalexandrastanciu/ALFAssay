#!/bin/bash
header_file=$1
output_file=$2
echo "$1"
echo "$2"

FIRST=

for FILE in $1
do
        exec 5<"$FILE" # Open file
        read LINE <&5 # Read first line
        [ -z "$FIRST" ] && echo "$LINE" # Print it only from first file
        FIRST="no"

        cat <&5 # Print the rest directly to standard output
        exec 5<&- # Close file
        # Redirect stdout for this section into file.out
done > $2

echo "Done bash script"