#!/usr/bin/env python

from __future__ import with_statement

import sys
import os
import codecs
import gzip
import logging
import json

from collections import OrderedDict, namedtuple
from logging import error, warn, info

try:
    import xml.etree.ElementTree as ET
except ImportError:
    import cElementTree as ET

utf8_stdout = codecs.getwriter('utf8')(sys.stdout)

output_count, skipped_count = 0, 0

def argparser():
    import argparse

    ap=argparse.ArgumentParser(description='Extract texts from PubMed XML.')
    ap.add_argument('-j', '--json', default=False, action='store_true',
                    help='Output JSON')
    ap.add_argument('-gt', '--PMID-greater-than', metavar='PMID', default=None,
                    help='Only process citations with PMIDs > PMID.')
    ap.add_argument('-lt', '--PMID-lower-than', metavar='PMID', default=None,
                    help='Only process citations with PMIDs < PMID.')
    ap.add_argument('-sa', '--single-line-abstract', default=False,
                    action='store_true', help='Output abstract on single line.')
    ap.add_argument('-mh', '--mesh-headings', default=False,
                    action='store_true', help='Output MeSH headings.')
    ap.add_argument('-mt', '--mesh-trees', default=False, action='store_true',
                    help='Output expanded MeSH trees (implies -mh).')
    ap.add_argument('-na', '--no-abstract', default=False, action='store_true',
                    help='Do not output abstracts.')
    ap.add_argument('-nt', '--no-title', default=False, action='store_true',
                    help='Do not output titles.')
    ap.add_argument('-nc', '--no-colon', default=False, action='store_true',
                    help='Do not add a colon to structured abstract headings.')
    ap.add_argument('-o', '--output-dir', metavar='DIR', default='texts',
                    help='Output directory (default "texts\", "-" for stdout)')
    ap.add_argument('-v', '--verbose', default=False, action='store_true',
                    help='Verbose output.')
    ap.add_argument('files', metavar='FILE', nargs='+',
                    help='Input PubMed distribution XML file(s).')
    return ap

class Citation(object):
    """Represents a PubMed citation."""

    def __init__(self, PMID, title, sections, mesh):
        self.PMID = PMID
        self.title = title
        self.sections = sections
        self.mesh = mesh

    def abstract_text(self, options=None):
        section_texts = [s.text(options) for s in self.sections]
        space = '\n' if not options or not options.single_line_abstract else ' '
        return space.join(section_texts)

    def text(self, options=None):
        lines = []
        if not options or not options.no_title:
            lines.append(self.title)
        if not options or not options.no_abstract:
            lines.append(self.abstract_text(options))
        if options and options.mesh_headings:
            lines.append('\t'.join(['MeSH Terms:'] +
                                   [m.text(options) for m in self.mesh]))
        return '\n'.join(lines)

    def to_dict(self, options=None):
        obj = { '_id': self.PMID }
        if not options or not options.no_title:
            obj['title'] = self.title
        if not options or not options.no_abstract:
            obj['abstract'] = [s.to_dict(options) for s in self.sections]
        if options and options.mesh_headings:
            obj['mesh'] = [m.to_dict(options) for m in self.mesh]
        return obj

    @classmethod
    def from_xml(cls, element):
        PMID = find_only(element, 'PMID').text
        article = find_only(element, 'Article')
        title = find_only(article, 'ArticleTitle').text
        abstract = find_abstract(element, PMID)
        if abstract is None:
            abstractTexts = []
        else:
            abstractTexts = abstract.findall('AbstractText')
            assert len(abstractTexts), 'ERROR: no <AbstractText> for %s' % PMID
        sections = [AbstractSection.from_xml(a, PMID) for a in abstractTexts]
        mesh_headings = find_mesh_headings(element, PMID)
        mesh = [MeshHeading.from_xml(h) for h in mesh_headings]
        return cls(PMID, title, sections, mesh)

class EmptySection(Exception):
    pass

class AbstractSection(object):
    """Represents a section of an abstract with text and an optional label."""

    def __init__(self, text, label=None):
        self._text = text if text is not None else ''
        self.label = label if label is not None else ''

    def _separator(self, options=None):
        if not self.label or not self._text:
            return ''
        colon = ':' if not options or not options.no_colon else ''
        space = '\n' if not options or not options.single_line_abstract else ' '
        return colon + space

    def text(self, options=None):
        return self.label + self._separator(options) + self._text

    def to_dict(self, options=None):
        label = self.label
        label += ':' if label and (not options or not options.no_colon) else ''
        return {
            'label': label,
            'text': self._text
        }

    @classmethod
    def from_xml(cls, element, PMID):
        if element.text and element.text.strip() != '':
            text = element.text
        else:
            warn('empty text for <AbstractText>s in %s' % PMID)
            text = ''
        label = element.attrib.get('Label')
        # The special Label "UNLABELLED" is interpreted as empty.
        if label == 'UNLABELLED':
            label = ''
        # Empty text and label would imply an empty section. Refuse
        # to create such aberrations.
        if not text and not label:
            raise EmptySection('empty unlabelled <AbstractText> in %s' % PMID)
        return cls(text, label)

Descriptor = namedtuple('Descriptor', 'id name major')
Qualifier = namedtuple('Qualifier', 'id name major')

class MeshHeading(object):
    """Represents a MeSH heading with a Descriptor and optional Qualifiers."""

    def __init__(self, descriptor, qualifiers):
        self.descriptor = descriptor
        self.qualifiers = qualifiers

    def descriptor_qualifier_pairs(self):
        quals = self.qualifiers if self.qualifiers else [None]
        return [(self.descriptor, q) for q in quals]

    def heading_texts(self):
        texts = []
        for d, q in self.descriptor_qualifier_pairs():
            did = d.id + ('*' if d.major else '')
            dname = d.name + ('*' if d.major else '')
            if q is None:
                texts.append('%s (%s)' % (did, dname))
            else:
                qid = q.id + ('*' if q.major else '')
                qname = q.name + ('*' if q.major else '')
                texts.append('%s/%s (%s/%s)' % (did, qid, dname, qname))
        return texts

    def tree_numbers(self):
        uid_to_node, treenum_name = get_mesh_data()
        expanded = OrderedDict()
        # TODO: trace major topics through ancestor expansion
        for d, q in self.descriptor_qualifier_pairs():
            for treenum in uid_to_node[d.id]['tnum']:
                for t in mesh_ancestors(treenum):
                    expanded[(t, q)] = True
        return [tree_number_text(t, q, treenum_name) for t, q in expanded]

    def text(self, options=None):
        if not options or not options.mesh_trees:
            return '\t'.join(self.heading_texts())
        else:
            return '\t'.join(self.tree_numbers())

    def to_dict(self, options=None):
        return {
            'descriptor': {
                'id': self.descriptor.id,
                'name': self.descriptor.name,
                'major': self.descriptor.major
            },
            'qualifiers': [
                {
                    'id': qual.id,
                    'name': qual.name,
                    'major': qual.major
                }
                for qual in self.qualifiers
            ]
        }

    @classmethod
    def from_xml(cls, element):
        desc = find_only(element, 'DescriptorName')
        id_ = desc.attrib.get('UI')
        major = (desc.attrib.get('MajorTopicYN') == 'Y')
        descriptor = Descriptor(id_, desc.text, major)
        qualifiers = []
        for qual in element.findall('QualifierName'):
            id_ = qual.attrib.get('UI')
            major = (qual.attrib.get('MajorTopicYN') == 'Y')
            qualifiers.append(Qualifier(id_, qual.text, major))
        return cls(descriptor, qualifiers)

def find_only(element, match):
    """Return the only matching child of the given element.

    Fail on assert if there are no or multiple matches.
    """
    found = element.findall(match)
    assert len(found) == 1, 'Error: expected 1 %s, got %d' % (match, len(found))
    return found[0]

def find_abstract(citation, PMID):
    """Return the Abstract element for given Article, or None if none."""
    # basic case: exactly one <Abstract> in <Article>.
    article = find_only(citation, 'Article')
    abstracts = article.findall('Abstract')
    assert len(abstracts) in (0,1), 'ERROR: %d abstracts for PMID %s' % \
        (len(abstracts, PMID))
    abstract = None if not abstracts else abstracts[0]

    # if there's no <Abstract>, look for <OtherAbstract> in the
    # citation (*not* the article).
    if abstract is None:
        otherAbstracts = citation.findall('OtherAbstract')
        # This happens a few times.
        if len(otherAbstracts) > 1:
            warn('%d "other" abstracts for PMID %s. Only using first.' %
                 (len(otherAbstracts), PMID))
        if otherAbstracts != []:
            abstract = otherAbstracts[0]

    return abstract

def find_mesh_headings(citation, PMID, options=None):
    """Return list of MeshHeading elements in given citation."""
    if options and not options.mesh_headings:
        return []     # avoid unnecessary load
    heading_lists = citation.findall('MeshHeadingList')
    if not heading_lists:
        info('no MeSH headings for %s' % PMID)
        return []
    assert len(heading_lists) == 1, 'Multiple MeshHeadingLists for %s' % PMID
    headings = heading_lists[0]
    return headings.findall('MeshHeading')

def tree_number_text(treenum, qualifier, treenum_to_name):
    """Return human-readable text for MeSH treenumber."""
    num = treenum
    text = treenum_to_name[treenum]
    if qualifier is not None:
        num += '/' +qualifier.id
        text += '/'+ qualifier.name
    return '%s (%s)' % (num, text)

def get_mesh_data():
    from meshdata import mesh
    if get_mesh_data._cache is None:
        treenum_to_name = {}
        for uid, obj in mesh.iteritems():
            name = obj['name']
            for tnum in obj['tnum']:
                treenum_to_name[tnum] = name
        get_mesh_data._cache = mesh, treenum_to_name
    return get_mesh_data._cache
get_mesh_data._cache = None

def mesh_ancestors(treenum):
    """Return ancestor tree numbers for given MeSH tree number."""
    # MeSH tree numbers have dotted forms such as "H02.403.640", and
    # ancestor numbers can be generated by removing later
    # dot-separated substrings. As an exception, the top-level ID
    # consists of just the first letter.
    parts = treenum.split('.')
    return [treenum[0]] + ['.'.join(parts[:i+1]) for i in range(len(parts))]

def skip_pmid(PMID, options):
    """Return True if PMID should be skipped by options, False otherwise."""
    PMID = int(PMID)
    if ((options.PMID_greater_than is not None and
         PMID <= options.PMID_greater_than) or
        (options.PMID_lower_than is not None and
         PMID >= options.PMID_lower_than)):
        warn('skipping %d (limits %d-%d)' %
             (PMID, options.PMID_greater_than, options.PMID_lower_than))
        return True
    else:
        return False

def write_citation(directory, citation, options):
    if not options.json:
        text = citation.text(options)
    else:
        text = json.dumps(citation.to_dict(options), sort_keys=True,
                          indent=2, separators=(',', ': '))
    if directory is None:
        print >> utf8_stdout, text
    else:
        suffix = '.txt' if not options.json else '.json'
        fn = os.path.join(directory, citation.PMID+suffix)
        with codecs.open(fn, 'wt', encoding='utf-8') as out:
            print >> out, text

def strip_extensions(fn):
    """Strip all extensions from file name."""
    while True:
        fn, ext = os.path.splitext(fn)
        if not ext:
            break
    return fn

def make_output_directory(fn, options):
    # create a directory for this package; we don't want to have all
    # the files in a single directory.
    base = strip_extensions(os.path.basename(fn))
    directory = os.path.join(options.output_dir, base)
    try:
        os.makedirs(directory)
    except OSError, e:
        error('Failed to create %s: %s' % (directory, str(e)))
        raise
    return directory

def process_stream(stream, outdir, options):
    global output_count, skipped_count

    for event, element in stream:
        if event != 'end' or element.tag != 'MedlineCitation':
            continue

        PMID = find_only(element, 'PMID').text
        if skip_pmid(PMID, options):
            skipped_count += 1
            element.clear()    # Won't need this
            continue

        citation = Citation.from_xml(element)
        write_citation(outdir, citation, options)
        output_count += 1

        element.clear()

def process(fn, options):
    if options.output_dir == '-':
        outdir = None    # use STDOUT
    else:
        outdir = make_output_directory(fn, options)

    if not fn.endswith('.gz'):
        process_stream(ET.iterparse(fn), outdir, options)
    else:
        with gzip.GzipFile(fn) as stream:
            process_stream(ET.iterparse(stream), outdir, options)

def process_options(argv):
    options = argparser().parse_args(argv[1:])
    if options.verbose:
        logging.getLogger().setLevel(logging.INFO)
    if options.mesh_trees:
        options.mesh_headings = True     # -mt implies -mh
    if options.no_title and options.no_abstract and not options.mesh_headings:
        error('nothing to output (-nt and -na and not -mh)')
        return None
    if options.PMID_greater_than is not None:
        options.PMID_greater_than = int(options.PMID_greater_than)
    if options.PMID_lower_than is not None:
        options.PMID_lower_than = int(options.PMID_lower_than)
    return options

def main(argv):
    global output_count, skipped_count

    options = process_options(argv)
    if options is None:
        return 1

    for fn in options.files:
        process(fn, options)

    print 'Done. Output texts for %d PMIDs, skipped %d.' % (output_count,
                                                            skipped_count)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
