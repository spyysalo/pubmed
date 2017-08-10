#!/bin/bash

# Convert output of extractTIABs.py --json --metadata --no-title --no-abstract
# into YEAR<TAB>PUBDATE format.

cat medline16n0*.meta \
    | egrep '^[[:space:]]*"(id|pubdate)":' \
    | perl -pe 's/"//g' \
    | perl -pe 's/^\s*id:\s*(\d+),?\s*/$1\t/; s/^\s*pubdate:\s*//'
