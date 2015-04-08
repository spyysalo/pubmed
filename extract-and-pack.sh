#!/bin/bash

# Extract texts from PubMed XML and package them as tar.gz.

set -e
set -u

DIR=`mktemp -d extract-XXX`

function cleanup {
    rm -rf "$DIR"
}
trap cleanup EXIT

if [[ $# < 1 ]]; then
    echo "Usage: $0 FILE [FILE [...]]"
    exit 1
fi

for f in "$@"; do
    echo -n "Extracting $f ... " >&2
    python extractTIABs.py -o "$DIR" "$f"
    echo "done." >&2

    t=`basename $f .xml.gz`
    t="${t%.xml}.tar.gz"
    if [[ -e "$t" ]]; then
	echo "Error: $t exists already, won't clobber."
	exit 1
    fi

    echo -n "Packing to $t ... " >&2
    tar czf "$t" -C "$DIR" .
    echo "done." >&2

    echo -n "Cleaning up ... " >&2
    rm -rf "$DIR"
    mkdir "$DIR"
    echo "done." >&2
done
