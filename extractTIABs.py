#!/usr/bin/env python

from __future__ import with_statement

import sys
import os
import codecs
import gzip
import logging
import tarfile
import json

from time import time
from StringIO import StringIO
from collections import OrderedDict, namedtuple
from logging import error, warn, info

from gtbtokenize import tokenize

try:
    import xml.etree.ElementTree as ET
except ImportError:
    import cElementTree as ET

utf8_stdout = codecs.getwriter('utf8')(sys.stdout)

output_count, skipped_count = 0, 0

month_abbr_map = {
    'Jan': '01',
    'Feb': '02',
    'Mar': '03',
    'Apr': '04',
    'May': '05',
    'Jun': '06',
    'Jul': '07',
    'Aug': '08',
    'Sep': '09',
    'Oct': '10',
    'Nov': '11',
    'Dec': '12',
}

class FormatError(Exception):
    pass

def argparser():
    import argparse

    ap=argparse.ArgumentParser(description='Extract texts from PubMed XML.')
    ap.add_argument('-a', '--ascii', default=False, action='store_true',
                    help='Map text to ASCII')
    ap.add_argument('-i', '--ids', metavar='FILE', default=None,
                    help='Only process citations with IDs in FILE.')
    ap.add_argument('-j', '--json', default=False, action='store_true',
                    help='Output JSON')
    ap.add_argument('-ha', '--has-abstract', default=False, action='store_true',
                    help='Only process citations with abstracts.')
    ap.add_argument('-gt', '--PMID-greater-than', metavar='PMID', default=None,
                    help='Only process citations with PMIDs > PMID.')
    ap.add_argument('-lt', '--PMID-lower-than', metavar='PMID', default=None,
                    help='Only process citations with PMIDs < PMID.')
    ap.add_argument('-sa', '--single-line-abstract', default=False,
                    action='store_true', help='Output abstract on single line.')
    ap.add_argument('-ii', '--include-id', default=False,
                    action='store_true', help='Include PMID in ouput.')
    ap.add_argument('-m', '--metadata', default=False, action='store_true',
                    help='Output metadata (e.g. publication date)')
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
    ap.add_argument('-ss', '--ssplit', default=False, action='store_true',
                    help='Perform sentence splitting.')
    ap.add_argument('-tt', '--tokenize', default=False, action='store_true',
                    help='Perform tokenization')
    ap.add_argument('-s', '--substances', default=False, action='store_true',
                    help='Output substances (chemicals).')
    ap.add_argument('-o', '--output-dir', metavar='DIR', default='texts',
                    help='Output directory (default "texts\", "-" for stdout)')
    ap.add_argument('-v', '--verbose', default=False, action='store_true',
                    help='Verbose output.')
    ap.add_argument('-z', '--tgz', default=False, action="store_true",
                    help='Output .tar.gz file')
    ap.add_argument('files', metavar='FILE', nargs='+',
                    help='Input PubMed distribution XML file(s).')
    return ap

def inner_text(element):
    """Return the catenated text of element and all subelements."""
    return ''.join(element.itertext())

class Citation(object):
    """Represents a PubMed citation."""

    def __init__(self, PMID, title, sections, mesh, chemicals, metadata):
        self.PMID = PMID
        self.title = title
        self.sections = sections
        self.mesh = mesh
        self.chemicals = chemicals
        self.metadata = metadata

    def abstract_text(self, options=None):
        section_texts = [s.text(options) for s in self.sections]
        space = '\n' if not options or not options.single_line_abstract else ' '
        return space.join(section_texts)

    def text(self, options=None):
        lines = []
        if options and options.include_id:
            lines.append(self.PMID)
        if not options or not options.no_title:
            lines.append(self.title)
        if not options or not options.no_abstract:
            lines.append(self.abstract_text(options))
        if options and options.mesh_headings:
            lines.append('\t'.join(['MeSH Terms:'] +
                                   [m.text(options) for m in self.mesh]))
        if options and options.substances:
            lines.append('\t'.join(['Substances:'] +
                                   [c.text(options) for c in self.chemicals]))
        if options and options.metadata:
            for k, v in self.metadata.items():
                lines.append('{}\t{}'.format(k, v))
        text = '\n'.join(lines)
        if not text.endswith('\n'):
            text = text + '\n'
        return text

    def to_dict(self, options=None):
        obj = { 'id': self.PMID }
        if not options or not options.no_title:
            obj['title'] = self.title
        if not options or not options.no_abstract:
            obj['abstract'] = [s.to_dict(options) for s in self.sections]
        if options and options.mesh_headings:
            obj['mesh'] = [m.to_dict(options) for m in self.mesh]
        if options and options.substances:
            obj['chemicals'] = [c.to_dict(options) for c in self.chemicals]
        if options and options.metadata:
            for k, v in self.metadata.items():
                k = k.lower()    # consistency w/JSON from esummary API
                assert k not in obj, 'internal error: metadata clobbers %s' % k
                obj[k] = v
        return obj

    @classmethod
    def from_xml(cls, element):
        PMID = find_only(element, 'PMID').text
        article = find_only(element, 'Article')
        title = inner_text(find_only(article, 'ArticleTitle'))
        if not title:
            # fallback
            try:
                title = inner_text(find_only(article, 'VernacularTitle'))
            except KeyError:
                pass
        if not title:
            warn('missing title for {}'.format(PMID))
            title = ''
        abstract = find_abstract(element, PMID)
        if abstract is None:
            abstractTexts = []
        else:
            abstractTexts = abstract.findall('AbstractText')
            assert len(abstractTexts), 'ERROR: no <AbstractText> for %s' % PMID
        sections = []
        for a in abstractTexts:
            try:
                sections.append(AbstractSection.from_xml(a, PMID))
            except EmptySection, e:
                warn(str(e))
        mesh_headings = find_mesh_headings(element, PMID)
        mesh = [MeshHeading.from_xml(h) for h in mesh_headings]
        chemical_list = find_chemicals(element, PMID)
        chemicals = [Chemical.from_xml(c) for c in chemical_list]
        metadata = find_metadata(element, PMID)
        return cls(PMID, title, sections, mesh, chemicals, metadata)

class EmptySection(Exception):
    pass

class AbstractSection(object):
    """Represents a section of an abstract with text and an optional label."""

    def __init__(self, text, label=None):
        self._text = text if text is not None else ''
        self.label = label if label is not None else ''

    def _colon(self, options=None):
        if options and options.no_colon:
            return ''
        elif options and options.tokenize:
            return ' :'
        else:
            return ':'

    def _separator(self, options=None):
        if not self.label:
            return ''
        colon = self._colon(options)
        if not self._text or self._text.isspace():
            space = ''
        elif not options or not options.single_line_abstract:
            space = '\n'
        else:
            space = ' '
        return colon + space

    def text(self, options=None):
        return self.label + self._separator(options) + self._text

    def to_dict(self, options=None):
        label = self.label
        label += self._colon(options) if label else ''
        return {
            'label': label,
            'text': self._text
        }

    @classmethod
    def from_xml(cls, element, PMID):
        text = inner_text(element)
        if not (text and text.strip() != ''):
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
            for treenum in uid_to_node[d.id]['treenums']:
                for t in mesh_ancestors(treenum):
                    expanded[(t, q)] = True
        return [tree_number_text(t, q, treenum_name) for t, q in expanded]

    def text(self, options=None):
        if not options or not options.mesh_trees:
            return '\t'.join(self.heading_texts())
        else:
            return '\t'.join(self.tree_numbers())

    def to_dict(self, options=None):
        if options.mesh_trees:
            raise NotImplementedError('-mt not supported in JSON output')
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

class Chemical(object):
    """Represents a Chemical entry with a ID, name, and registry number."""

    def __init__(self, id_, name, regnum):
        self.id = id_
        self.name = name
        self.regnum = regnum

    def text(self, options=None):
        # Note: registry number not represented
        return '%s (%s)' % (self.id, self.name)

    def to_dict(self, options=None):
        return {
            'id': self.id,
            'name': self.name,
            'regnum': self.regnum
        }

    @classmethod
    def from_xml(cls, element):
        id_name = find_only(element, 'NameOfSubstance')
        id_ = id_name.attrib.get('UI')
        name = id_name.text
        regnum = find_only(element, 'RegistryNumber').text
        return cls(id_, name, regnum)

def find_only(element, match):
    """Return the only matching child of the given element.

    Raise KeyError if no match and FormatError if multiple matches.
    """
    found = element.findall(match)
    if not found:
        raise KeyError('Error: expected 1 %s, got %d' % (match, len(found)))
    elif len(found) > 1:
        raise FormatError('Error: expected 1 %s, got %d' % (match, len(found)))
    else:
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

def find_chemicals(citation, PMID, options=None):
    """Return list of Chemical elements in given citation."""
    if options and not options.substances:
        return []    # avoid unnecessary load
    chemical_lists = citation.findall('ChemicalList')
    if not chemical_lists:
        info('No chemical list for %s' % PMID)
        return []
    assert len(chemical_lists) == 1, 'Multiple ChemicalLists for %s' % PMID
    chemicals = chemical_lists[0]
    return chemicals.findall('Chemical')

def date_string(element, PMID, expect_full_date=True):
    """Format <Date*> element content as string."""
    values = {}
    for e in element:
        if e.tag in values:
            warn('duplicate %s in %s in %s' % (e.tag, element.tag, PMID))
        else:
            values[e.tag] = e.text
    date = values.pop('MedlineDate', None)
    year = values.pop('Year', None)
    month = values.pop('Month', None)
    day = values.pop('Day', None)
    values.pop('Season', None)    # ignore Season
    if values:
        warn('extra data in %s in %s: %s' % (
            element.tag, PMID,
            ' '.join(['%s="%s"' % (k, v) for k, v in values.items()])))
    if month is not None and not month.isdigit():
        if month in month_abbr_map:
            month = month_abbr_map[month]    # "Jan" -> "01" etc.
        else:
            warn('unexpected month format %s in %s' % (month, PMID))
    if date is not None:
        if year or month or day:
            warn('extra data w/MedlineDate in %s in %s' % (element.tag, PMID))
        return date
    elif year is None:
        warn('missing <Year> in %s in %s' % (element.tag, PMID))
        return ''
    elif month is None:
        if expect_full_date:
            warn('missing <Month> in %s in %s' % (element.tag, PMID))
        return year
    elif day is None:
        if expect_full_date:
            warn('missing <Day> in %s in %s' % (element.tag, PMID))
        return '%s-%s' % (year, month)
    else:
        return '%s-%s-%s' % (year, month, day)    # https://xkcd.com/1179/

def find_metadata(citation, PMID, options=None):
    """Return dictionary of citation metadata."""
    if options and not options.metadata:
        return {}    # avoid unnecessary load
    metadata = OrderedDict()

    # PubMed citation dates (if present)
    for date in ('DateCreated', 'DateCompleted'):
        try:
            element = find_only(citation, date)
        except KeyError:
            info('missing <%s> in %s' % (date, PMID))
            continue
        metadata[date] = date_string(element, PMID)

    # Article dates
    article = find_only(citation, 'Article')
    pubdate = find_only(article, './/PubDate')    # recursive
    metadata['PubDate'] = date_string(pubdate, PMID, expect_full_date=False)
    for e in article.findall('ArticleDate'):
        datetype = e.attrib.get('DateType')
        if datetype == 'Electronic':
            metadata['EPubDate'] = date_string(e, PMID)
        else:
            warn('unknown <ArticleDate DateType="%s"> in %s' % (datetype, PMID))
    return metadata

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
            for tnum in obj['treenums']:
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
        info('skipping %d (limits %d-%d)' %
             (PMID, options.PMID_greater_than, options.PMID_lower_than))
        return True
    elif options.ids is not None and PMID not in options.ids:
        info('skipping %d (not in given IDs)' % PMID)
        return True
    else:
        return False

def skip_citation(element, options):
    """Return True if citation should be skipped by options, False otherwise."""

    PMID = find_only(element, 'PMID').text
    if skip_pmid(PMID, options):
        return True
    elif options.has_abstract and find_abstract(element, PMID) is None:
        info('skipping %s (no abstract)' % PMID)
        return True
    else:
        return False

def to_ascii(s):
    """Map string to ASCII"""
    import unicode2ascii
    if not s:
        return s
    if to_ascii.mapping is None:
        mapfn = os.path.join(os.path.dirname(__file__), 'entities.dat')
        with codecs.open(mapfn, encoding='utf-8') as f:
            to_ascii.mapping = unicode2ascii.read_mapping(f, mapfn)
    out = StringIO()
    unicode2ascii.process([s], out, to_ascii.mapping)
    return out.getvalue()
to_ascii.mapping = None

def write_to_ascii_statistics(out=sys.stderr):
    from unicode2ascii import missing_mapping
    if not missing_mapping:
        return    # nothing missing
    print >> out, "Characters without mapping\t%d" % sum(missing_mapping.values())
    sk = missing_mapping.keys()
    sk.sort(lambda a,b : cmp(missing_mapping[b],missing_mapping[a]))
    for c in sk:
        try:
            print >> out, "\t%.4X\t%s\t%d" % (ord(c), c.encode("utf-8"), missing_mapping[c])
        except:
            print >> out, "\t%.4X\t?\t%d" % (ord(c), missing_mapping[c])

def citation_to_ascii(citation):
    """Map citation text content to ASCII"""
    citation.title = to_ascii(citation.title)
    for section in citation.sections:
        section._text = to_ascii(section._text)
        section.label = to_ascii(section.label)

def citation_ssplit(citation):
    """Split sentences in citation text content."""
    if citation_ssplit.ssplitter is None:
        from ssplit import ssplitter
        citation_ssplit.ssplitter = ssplitter
    ssplitter = citation_ssplit.ssplitter
    citation.title = ssplitter(citation.title)
    for section in citation.sections:
        section._text = ssplitter(section._text)
        section.label = ssplitter(section.label)
citation_ssplit.ssplitter = None

def tokenize_multiline(text):
    return '\n'.join(tokenize(s) for s in text.split('\n'))

def citation_tokenize(citation):
    citation.title = tokenize_multiline(citation.title)
    for section in citation.sections:
        section._text = tokenize_multiline(section._text)
        section.label = tokenize_multiline(section.label)

def save_in_tar(tar, name, text):
    info = tar.tarinfo(name)
    sio = StringIO(text.encode('utf-8'))
    sio.seek(0)
    info.size = sio.len
    info.mtime = time()
    tar.addfile(info, sio)

def write_citation(directory, name, outfile, citation, options):
    if options.ascii:
        citation_to_ascii(citation)
    if options.ssplit:
        citation_ssplit(citation)
    if options.tokenize:
        citation_tokenize(citation)
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
        if options.tgz:
            fn = os.path.join(os.path.basename(name).split('.')[0],
                              os.path.basename(fn))
            save_in_tar(outfile, fn, text)
        else:
            with codecs.open(fn, 'wt', encoding='utf-8') as out:
                out.write(text)

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
    if not options.tgz:
        directory = os.path.join(options.output_dir, base)
    else:
        # tgz: create output_dir only, no subdirs
        directory = options.output_dir
        if os.path.isdir(directory):
            return directory
    try:
        os.makedirs(directory)
    except OSError, e:
        error('Failed to create %s: %s' % (directory, str(e)))
        raise
    return directory

def tarname(outdir, name):
    base = os.path.basename(name).split('.')[0]
    return os.path.join(outdir, base + '.tar.gz')

def process_stream(stream, name, outdir, options):
    global output_count, skipped_count

    if options.tgz:
        outfile = tarfile.open(tarname(outdir, name), 'w:gz') # TODO use `with`
    else:
        outfile = None

    for event, element in stream:
        if event != 'end' or element.tag != 'MedlineCitation':
            continue

        if skip_citation(element, options):
            skipped_count += 1
            element.clear()    # Won't need this
            continue

        citation = Citation.from_xml(element)
        write_citation(outdir, name, outfile, citation, options)
        output_count += 1

        element.clear()

    if options.tgz:
        outfile.close()

def process(fn, options):
    if options.output_dir == '-':
        outdir = None    # use STDOUT
    else:
        outdir = make_output_directory(fn, options)

    if not fn.endswith('.gz'):
        process_stream(ET.iterparse(fn), fn, outdir, options)
    else:
        with gzip.GzipFile(fn) as stream:
            process_stream(ET.iterparse(stream), fn, outdir, options)

def read_ids(fn):
    ids = set()
    with open(fn) as f:
        for ln, l in enumerate(f, start=1):
            l = l.strip()
            try:
                i = int(l)
            except ValueError:
                raise ValueError('Error on line %d in %s: not an ID: %s' % (
                    ln, fn, l))
            ids.add(i)
    info('read %d IDs from %s' % (len(ids), fn))
    return ids

def process_options(argv):
    options = argparser().parse_args(argv[1:])
    if options.verbose:
        logging.getLogger().setLevel(logging.INFO)
    if options.mesh_trees:
        options.mesh_headings = True     # -mt implies -mh
    if options.tokenize and not options.ssplit:
        # Tokenizer assumes sentence-split input
        warn('--ssplit recommended with --tokenize')
    if (options.no_title and options.no_abstract and
        not (options.mesh_headings or options.include_id or options.metadata)):
        error('nothing to output (-nt and -na without other output options)')
        return None
    if options.PMID_greater_than is not None:
        options.PMID_greater_than = int(options.PMID_greater_than)
    if options.PMID_lower_than is not None:
        options.PMID_lower_than = int(options.PMID_lower_than)
    if options.ids is not None:
        options.ids = read_ids(options.ids)
    return options

def main(argv):
    global output_count, skipped_count

    options = process_options(argv)
    if options is None:
        return 1

    for fn in options.files:
        try:
            process(fn, options)
        except:
            error('Failed to process %s' % fn)
            raise

    if options.ascii:
        write_to_ascii_statistics(sys.stderr)

    print >> sys.stderr, 'Done. Output data for %d PMIDs, skipped %d.' % (
        output_count, skipped_count)

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
