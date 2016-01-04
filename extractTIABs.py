#!/usr/bin/env python

from __future__ import with_statement

import sys
import re
import os
import gzip

try:
    import xml.etree.ElementTree as ET
except ImportError: 
    import cElementTree as ET   

options = None
output_count, skipped_count = 0,0

def argparser():
    import argparse

    ap=argparse.ArgumentParser(description="Extract per-document title and abstract texts from PubMed distribution XML file.")
    ap.add_argument("-o", "--output-dir", metavar="DIR", default="texts", help="Output directory (default \"texts/\", \"-\" for stdout)")
    ap.add_argument("-gt", "--PMID-greater-than", metavar="PMID", default=None, help="Only process citations with PMIDs greater than the given value.")
    ap.add_argument("-lt", "--PMID-lower-than", metavar="PMID", default=None, help="Only process citations with PMIDs lower than the given value.")
    ap.add_argument("-sa", "--single-line-abstract", default=False, action="store_true", help="Always output abstract on single line.")
    ap.add_argument("-na", "--no-abstract", default=False, action="store_true", help="Don't output abstracts (titles only).")
    ap.add_argument("-nt", "--no-title", default=False, action="store_true", help="Don't output titles (abstracts only).")
    ap.add_argument("-nc", "--no-colon", default=False, action="store_true", help="Don't add a colon to structured abstract headings.")
    ap.add_argument("-v", "--verbose", default=False, action="store_true", help="Verbose output.")
    ap.add_argument("files", metavar="FILE", nargs="+", help="Input PubMed distribution XML file(s).")

    return ap

def find_only(element, match):
    """Return the only matching child of the given element.

    Fail on assert if there are no or multiple matches.
    """
    found = element.findall(match)
    assert len(found) == 1, 'Error: expected 1 %s, got %d' % (match, len(found))
    return found[0]

def find_abstract(citation, PMID, options):
    """Return the Abstract element for given Article, or None if none."""
    # basic case: exactly one <Abstract> in <Article>.
    article = find_only(citation, 'Article')
    abstracts = article.findall('Abstract')
    assert len(abstracts) in (0,1), 'ERROR: %d abstracts for PMID %s' % \
        (len(abstracts, PMID.text))
    abstract = None if not abstracts else abstracts[0]

    # if there's no <Abstract>, look for <OtherAbstract> in the
    # citation (*not* the article).
    if abstract is None:
        otherAbstracts = citation.findall('OtherAbstract')
        # This happens a few times.
        if len(otherAbstracts) > 1 and options.verbose:
            print >> sys.stderr, "NOTE: %d 'other' abstracts for PMID %s. Only printing first." % (len(otherAbstracts), PMID.text)
        if otherAbstracts != []:
            abstract = otherAbstracts[0]

    return abstract

def get_abstract_text(abstract, PMID, options):
    """Return the text of the given Abstract element."""
    if abstract is None:
        return ''

    # Abstract should contain an AbstractText
    abstractTexts = abstract.findall('AbstractText')
    assert len(abstractTexts) != 0, "ERROR: %d abstract texts for PMID %s" % \
                                 (len(abstractTexts), PMID.text)

    # Basic case: exactly one <AbstractText>; just return the content.
    if len(abstractTexts) == 1:
        return abstractTexts[0].text

    # Recent versions of PubMed data may contain multiple AbstractText
    # elements for structured abstracts. In these cases, "label"
    # attributes give structured abstract section headers and should
    # be combined into the text.
    assert len(abstractTexts) > 1, "INTERNAL ERROR"
    if options.verbose:
        print >> sys.stderr, "NOTE: multiple <AbstractText>s for %s" % PMID.text

    sectionTexts = []
    for at in abstractTexts:
        # there may be instances of empty <AbstractText> elements in
        # the data (see e.g. PMID 20619000 in the PubMed 2012
        # baseline). Skip those with the special "empty" label
        # "UNLABELLED" entirely; print the label only for the rest.
        if ((at.text is None or at.text.strip() == "") and
            at.attrib.get("Label","") == "UNLABELLED"):
            if options.verbose:
                print >> sys.stderr, "NOTE: skipping empty <AbstractText>s with label \"UNLABELLED\" in %s" % PMID.text
            continue

        t = ""
        if "Label" not in at.attrib:
            print >> sys.stderr, "Warning: missing 'Label' for multiple <AbstractText>s in %s" % PMID.text
        elif at.attrib["Label"] == "UNLABELLED":
            if options.verbose:
                print >> sys.stderr, "NOTE: skipping <AbstractText> Label \"UNLABELLED\" in %s" % PMID.text
        else:
            t = at.attrib["Label"]
            if not options.no_colon:
                t += ":"

        if at.text is None or at.text.strip() == "":
            print >> sys.stderr, "NOTE: empty text for one of multiple <AbstractText>s in %s" % PMID.text
        else:
            if not t:
                t = at.text
            elif options.single_line_abstract:
                t += " " + at.text
            else:
                t += "\n"+ at.text
        sectionTexts.append(t)

    if options.single_line_abstract:
        return " ".join(sectionTexts)
    else:
        return "\n".join(sectionTexts)

def skip_pmid(PMID, options, out=None):
    """Return True if PMID should be skipped by options, False otherwise."""
    if out is None:
        out = sys.stderr
    if ((options.PMID_greater_than is not None and
         PMID <= options.PMID_greater_than) or
        (options.PMID_lower_than is not None and
         PMID >= options.PMID_lower_than)):
        if options.verbose:
            print >> out, "Note: skipping %d" % PMID,
            if options.PMID_greater_than is not None:
                print >> out, "(lower limit %d)" % options.PMID_greater_than,
            if options.PMID_lower_than is not None:
                print >> out, "(upper limit %d)" % options.PMID_lower_than,
            print >> out
        return True
    else:
        return False

def process(fn):
    global options, output_count, skipped_count

    # create a directory for this package; we don't want to have all
    # the files in a single directory.
    filename = os.path.basename(fn)
    m = re.match(r'([A-Za-z0-9_-]+)\.xml(?:\.gz)?$', filename)
    assert m, "ERROR: unexpected filename '%s'" % filename
    filenamebase = m.group(1)

    if options.output_dir != '-':
        outputDir = os.path.join(options.output_dir, filenamebase)
        # TODO: fail gracefully on problems
        os.mkdir(outputDir)
    else:
        outputDir = None    # use STDOUT

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
        PMID = find_only(citation, 'PMID')

        # if a PMID range has been specified, check that we're in it
        if skip_pmid(int(PMID.text), options):
            skipped_count += 1
            citation.clear()    # Won't need this
            continue

        # likewise, there should be exactly one Article child and
        # Article should have a single ArticleTitle. (Abstract is a
        # bit trickier.)
        article = find_only(citation, 'Article')
        articleTitle = find_only(article, 'ArticleTitle')
        abstract = find_abstract(citation, PMID, options)

        # We've got all the elements we need. Now we just need the texts
        textPMID = PMID.text
        textTitle = articleTitle.text
        textAbstract = get_abstract_text(abstract, PMID, options)

        # sanity
        assert re.match(r'^\d+$', textPMID), \
            "ERROR: unexpected PMID format: '%s'" % textPMID

        if outputDir is not None:
            # output title and abstract into a file
            outputFile = os.path.join(outputDir, textPMID + ".txt")
            out = open(outputFile, "w")
        else:
            out = sys.stdout

        if not options.no_title:
            print >> out, textTitle.encode("UTF-8")
        if not options.no_abstract:
            if textAbstract:
                print >> out, textAbstract.encode("UTF-8")
            elif options.verbose:
                print >> sys.stderr, "No abstract for %s" % textPMID

        if outputDir is not None:
            out.close()

        output_count += 1
        
        # finally, clear out the used data; we don't need it anymore.
        citation.clear()

    # if we were wrapping a .gz, close the GzipFile
    if gzipped:
        input.close()

def main(argv):
    global options, output_count, skipped_count
    options = arg = argparser().parse_args(argv[1:])

    if options.no_title and options.no_abstract:
        print >> sys.stderr, "Error: both --no-abstract and --no-title given"
        return 1

    if options.PMID_greater_than is not None:
        options.PMID_greater_than = int(options.PMID_greater_than)
    if options.PMID_lower_than is not None:
        options.PMID_lower_than = int(options.PMID_lower_than)

    for fn in arg.files:
        process(fn)

    print >> sys.stderr, "Done. Output texts for %d PMIDs, skipped %d." % (output_count, skipped_count)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
