#!/bin/bash

# Script to find files with the same name in two different directories
# Usage: ./check_dp.sh /path/to/folder1 /path/to/folder2

# Check if two arguments are provided
if [ $# -ne 2 ]; then
    echo "Error: Two directory paths are required."
    echo "Usage: $0 /path/to/folder1 /path/to/folder2"
    exit 1
fi

# Check if both arguments are valid directories
if [ ! -d "$1" ]; then
    echo "Error: '$1' is not a valid directory."
    exit 1
fi

if [ ! -d "$2" ]; then
    echo "Error: '$2' is not a valid directory."
    exit 1
fi

# Get the absolute paths
FOLDER1=$(realpath "$1")
FOLDER2=$(realpath "$2")

echo "Checking for duplicate filenames between:"
echo "  $FOLDER1"
echo "  $FOLDER2"
echo

# Get list of files in both directories (just the filenames, not the full paths)
FILES1=$(find "$FOLDER1" -type f -printf "%f\n" | sort)
FILES2=$(find "$FOLDER2" -type f -printf "%f\n" | sort)

# Use comm to find common files
DUPLICATES=$(comm -12 <(echo "$FILES1") <(echo "$FILES2"))

# Count duplicates
DUP_COUNT=$(echo "$DUPLICATES" | grep -v '^$' | wc -l)

echo "Found $DUP_COUNT duplicate filename(s):"
echo

# Print the duplicates with their full paths
if [ $DUP_COUNT -gt 0 ]; then
    echo "$DUPLICATES" | while read -r filename; do
        if [ -n "$filename" ]; then
            echo "File: $filename"
            echo "  Found in: $FOLDER1/$(find "$FOLDER1" -name "$filename" | sed "s|$FOLDER1/||")"
            echo "  Found in: $FOLDER2/$(find "$FOLDER2" -name "$filename" | sed "s|$FOLDER2/||")"
            echo
        fi
    done
else
    echo "No duplicate filenames found."
fi

exit 0