#!/bin/bash

# Check PubMed distribution-style md5 checksums.

set -e
set -u

if [[ $# < 1 ]]; then
    echo "Usage: $0 DIR [DIR [...]]"
    exit 1
fi

for d in "$@"; do
    for f in "$d/"*.md5; do
	# assume that stripping the suffix ".md5" gives the file name that
	# the md5sum is for
	new=`md5sum < ${f%.md5} | perl -pe 's/ +- *$//'`
	old=`cat $f | perl -pe 's/^MD5 \(.*?\) = //'`
	if [[ "$new" != "$old" ]]; then
	    echo "Mismatch: $f: $new vs $old"
	fi
    done
done
