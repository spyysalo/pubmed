#!/usr/bin/env python

import sys
import re
import os

import gzip

import psyco

try:
    import xml.etree.ElementTree as ET
except ImportError: 
    import cElementTree as ET   

# Name of root directory under which extracted titles and abstracts
# should be stored
outputRoot = "texts/"

if len(sys.argv) < 2:
    print >> sys.stderr, "Usage:", sys.argv[0], "[-gt PMID] FILES"
    sys.exit(1)

# Quick & dirty support for limiting extraction to a range of
# PMIDs greater than the supplied lower limit.
lower_limit_PMID = None
if len(sys.argv) > 3 and sys.argv[1] == "-gt":
    lower_limit_PMID = int(sys.argv[2])
    sys.argv = [sys.argv[0]]+sys.argv[3:]

# maybe this will give us a bit of something
psyco.full()

output_count, skipped_count = 0,0

for fn in sys.argv[1:]:
    # create a directory for this package; we don't want to have all
    # the files in a single directory.
    filename = os.path.basename(fn)
    m = re.match(r'(medline\d\dn\d+)\.xml(?:\.gz)?$', filename)
    assert m, "ERROR: unexpected filename '%s'" % filename
    filenamebase = m.group(1)
    outputDir = outputRoot+filenamebase+"/"
    
    # TODO: fail gracefully on problems
    os.mkdir(outputDir)

    # if the extension suggests a zip, wrap
    input = fn
    gzipped = False
    if re.search('\.gz$', fn):
        input = gzip.GzipFile(fn)
        gzipped = True

    for event, element in ET.iterparse(input):
        # we're only interested in tracking citation end events
        if event != "end" or element.tag != "MedlineCitation":
            continue

        citation = element

        # the citation element should have exactly one PMID child
        PMIDs = citation.findall("PMID")
        assert len(PMIDs) == 1, "ERROR: expected 1 PMID, got %d" % len(PMIDs)
        PMID = PMIDs[0]

        # if a "greater than" range has been specified, check that
        # we're in the range
        if lower_limit_PMID is not None and int(PMID.text) <= lower_limit_PMID:
            print >> sys.stderr, "Note: skipping %s (lower limit %d)" % (PMID.text, lower_limit_PMID)
            skipped_count += 1
            # clear out the element; we're not going to use it.
            citation.clear()
            continue

        # likewise, there should be exactly one Article child
        articles = citation.findall("Article")
        assert len(articles) == 1, "ERROR: %d articles for PMID %s" % (len(articles), PMID.text)
        article = articles[0]

        # further, Article should have a single ArticleTitle
        articleTitles = article.findall("ArticleTitle")
        assert len(articleTitles) == 1, "ERROR: %d titles for PMID %s" % (len(articleTitles, PMID.text))
        articleTitle = articleTitles[0]

        # also, Article typically (but not always) contains an Abstract
        abstracts = article.findall("Abstract")
        assert len(abstracts) in (0,1), "ERROR: %d abstracts for PMID %s" % (len(abstracts, PMID.text))
        abstract = None
        if abstracts != []:
            abstract = abstracts[0]

        # if there's no Abstract, try to look for <OtherAbstract> in
        # the citation (not the article) element, which seems to be
        # used in some cases
        if abstract is None:
            otherAbstracts = citation.findall("OtherAbstract")
            # This happens a few times.
            if len(otherAbstracts) > 1:
                print >> sys.stderr, "NOTE: %d 'other' abstracts for PMID %s. Only printing first." % (len(otherAbstracts), PMID.text)
            if otherAbstracts != []:
                abstract = otherAbstracts[0]

        abstractText = None
        textAbstract = ""
        if abstract is not None:
            # if there's an Abstract, it should contain an AbstractText
            abstractTexts = abstract.findall("AbstractText")
            assert len(abstractTexts) != 0, "ERROR: %d abstract texts for PMID %s" % (len(abstractTexts), PMID.text)

            if len(abstractTexts) == 1:
                abstractText = abstractTexts[0]
                textAbstract = abstractText.text
            else:
                # recent versions of PubMed data may contain multiple
                # AbstractText elements for structured abstracts. In these
                # cases, "label" attributes give structured abstract
                # section headers and should be combined into the text.
                assert len(abstractTexts) > 1, "INTERNAL ERROR"
                print >> sys.stderr, "NOTE: multiple <AbstractText>s for %s" % PMID.text
                sectionTexts = []
                for at in abstractTexts:
                    t = ""
                    if "Label" in at.attrib:
                        t = at.attrib["Label"] + ":"
                    else:
                        print >> sys.stderr, "Warning: missing 'Label' for multiple <AbstractText>s in %s" % PMID.text
                    if at.text is not None:
                        t += " " + at.text
                    else:
                        print >> sys.stderr, "Warning: missing text for one of multiple <AbstractText>s in %s" % PMID.text
                    sectionTexts.append(t)
                textAbstract = " ".join(sectionTexts)

        # OK, we've got all we need. Now we just need the texts
        textPMID     = PMID.text
        textTitle    = articleTitle.text        

        # bit of sanity checking
        assert re.match(r'^\d+$', textPMID), "ERROR: unexpected characters in PMID: '%s'" % textPMID

        # output title and abstract into a file
        outputFile = outputDir + textPMID + ".txt"
        out = open(outputFile, "w")

        print >> out, textTitle.encode("UTF-8")
        if textAbstract:
            print >> out, textAbstract.encode("UTF-8")
        else:
            print >> sys.stderr, "No abstract for %s" % textPMID
        out.close()

        output_count += 1
        
        # finally, clear out the used data; we don't need it anymore.
        citation.clear()

    # if we were wrapping a .gz, close the GzipFile
    if gzipped:
        input.close()

print >> sys.stderr, "Done. Output texts for %d PMIDs, skipped %d." % (output_count, skipped_count)
