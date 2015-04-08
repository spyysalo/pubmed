# Intended to keep generally useful functions for dealing with PubMed XML data.

import sys

try:
    import xml.etree.ElementTree as ET
except ImportError: 
    import cElementTree as ET 

# Helper function for grabbing a text representation of the PMID,
# title and abstract of a given file. NOTE: this isn't really
# written as a library function; it may fail on assert and/or
# print stuff on stderr. Use at your own risk.

def getPmidTitleAbstract(fn):
    try:
        tree = ET.parse(fn)
    except IOError:
        print >> sys.stderr, "ERROR: failed to open PubMed XML file %s" % fn
        raise

    root = tree.getroot()

    # This largely follows code from the parsePubmed.py script in the
    # full-PubMed event extraction project (data/PubMed2009).

    citations = tree.findall(".//MedlineCitation")
    assert len(citations) == 1, "ERROR: expected 1 MedlineCitation, got %d" % len(citations)
    citation = citations[0]    

    # the citation element should have exactly one PMID child
    PMIDs = citation.findall("PMID")
    assert len(PMIDs) == 1, "ERROR: expected 1 PMID, got %d" % len(PMIDs)
    PMID = PMIDs[0]

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
            print >> sys.stderr, "NOTE: %d 'other' abstracts for PMID %s. Only using first." % (len(otherAbstracts), PMID.text)
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
            sectionTexts = []
            for at in abstractTexts:
                t = ""
                if "Label" in at.attrib:
                    t = at.attrib["Label"] + ": "
                else:
                    print >> sys.stderr, "Warning: missing 'Label' for multiple <AbstractText>s"
                t += at.text
                sectionTexts.append(t)
            textAbstract = " ".join(sectionTexts)
                    
    # OK, we've got all the elements we need. Grab the texts
    textPMID     = PMID.text
    textTitle    = articleTitle.text        
    return (textPMID, textTitle, textAbstract)
