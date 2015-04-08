#!/bin/bash

# Extract texts from PubMed XML and package them as tar.gz.

set -e
set -u

# http://stackoverflow.com/a/246128
SCRIPTDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

TMPDIR=`mktemp -d extract-XXX`

function cleanup {
    rm -rf "$TMPDIR"
}
trap cleanup EXIT

if [[ $# < 1 ]]; then
    echo "Usage: $0 FILE [FILE [...]]"
    exit 1
fi

for f in "$@"; do
    echo -n "Extracting $f ... " >&2
    python "$SCRIPTDIR/extractTIABs.py" -o "$TMPDIR" "$f"
    echo "done." >&2

    t=`basename $f .xml.gz`
    t="${t%.xml}.tar.gz"
    if [[ -e "$t" ]]; then
	echo "Error: $t exists already, won't clobber."
	exit 1
    fi

    echo -n "Packing to $t ... " >&2
    tar czf "$t" -C "$TMPDIR" .
    echo "done." >&2

    echo -n "Cleaning up ... " >&2
    rm -rf "$TMPDIR"
    mkdir "$TMPDIR"
    echo "done." >&2
done
