#!/usr/bin/env python

# Script to extract the title and abstract text from the given PubMed
# XML file.

import sys
import os
from pubmedtools import getPmidTitleAbstract

if len(sys.argv) != 2:
    print >> sys.stderr, "Usage:", sys.argv[0], "FILE"
    sys.exit(1)
pmxmlfn = sys.argv[1]

(PMID, title, abstract) = getPmidTitleAbstract(pmxmlfn)

sys.stdout.write(title.encode("utf-8")+"\n")

if abstract == "":
    print >> sys.stderr, "Warning:", sys.argv[0], "title only for PMID %s" % PMID
else:
    sys.stdout.write(abstract.encode("utf-8")+"\n")

