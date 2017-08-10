#!/bin/bash

# Convert output of extractTIABs.py --json --metadata --no-title --no-abstract
# into PMID<TAB>YEAR format, where YEAR is the publication year.

for f in "$@"; do
    egrep '^[[:space:]]*"(id|pubdate)":' $f \
	| perl -pe 's/"//g' \
	| perl -pe 's/^\s*id:\s*(\d+),?\s*/$1\t/; s/^\s*pubdate:\s*//' \
	| perl -pe 's/\t(?:spring|summer|autumn|fall|winter) (\d{4})/\t$1/i' \
	| perl -pe 's/(\t[0-9]{4})\b.*/$1/'
done
