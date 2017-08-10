#!/bin/bash

# Convert output of extractTIABs.py --json --metadata --no-title --no-abstract
# into YEAR<TAB>PUBDATE format.

for f in "$@"; do
    egrep '^[[:space:]]*"(id|pubdate)":' $f \
	| perl -pe 's/"//g' \
	| perl -pe 's/^\s*id:\s*(\d+),?\s*/$1\t/; s/^\s*pubdate:\s*//'
done
