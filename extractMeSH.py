#!/usr/bin/env python

# Extract unique ID, name and tree numbers from MeSH XML data.

import sys
import gzip
import json

import xml.etree.ElementTree as ET

from logging import info, warning, error


# MeSH top-level structure (from https://www.nlm.nih.gov/cgi/mesh/2016/MB_cgi),
# not explicitly found in the XML
meshtop = [
    ['A', 'Anatomy', ['A']],
    ['B', 'Organisms', ['B']],
    ['C', 'Diseases', ['C']],
    ['D', 'Chemicals and Drugs', ['D']],
    ['E', 'Analytical, Diagnostic and Therapeutic Techniques and Equipment', ['E']],
    ['F', 'Psychiatry and Psychology', ['F']],
    ['G', 'Phenomena and Processes', ['G']],
    ['H', 'Disciplines and Occupations', ['H']],
    ['I', 'Anthropology, Education, Sociology and Social Phenomena', ['I']],
    ['J', 'Technology, Industry, Agriculture', ['J']],
    ['K', 'Humanities', ['K']],
    ['L', 'Information Science', ['L']],
    ['M', 'Named Groups', ['M']],
    ['N', 'Health Care', ['N']],
    ['V', 'Publication Characteristics', ['V']],
    ['Z', 'Geographicals', ['Z']],
]


def argparser():
    import argparse
    ap=argparse.ArgumentParser(description="Extract data from MeSH XML")
    ap.add_argument('-t', '--top', default=False, action='store_true',
                    help='Include implicit top-level structure')
    ap.add_argument('-j', '--json', default=False, action='store_true',
                    help='Output JSON (default TSV)')
    ap.add_argument('-d', '--dict', default=False, action='store_true',
                    help='Output python dictionary (default TSV)')
    ap.add_argument('-b', '--brat-norm', default=False, action='store_true',
                    help='Output brat normalization format (default TSV)')
    ap.add_argument("file", metavar="FILE", help="Input MeSH XML.")
    return ap


def pretty_dumps(obj):
    return json.dumps(obj, sort_keys=True, indent=4, separators=(',', ': '))


def find_only(element, match):
    """Return the only matching child of the given element.

    Fail on assert if there are no or multiple matches.
    """
    found = element.findall(match)
    assert len(found) == 1, 'Error: expected 1 %s, got %d' % (match, len(found))
    return found[0]


def find_first(element, match):
    """Return the first matching child of the given element.

    Fail on assert if there are no matches.
    """
    found = element.findall(match)
    assert len(found) > 0, 'Error: did not find %s' % match
    return found[0]


def write_header(options, out=None):
    if out is None:
        out = sys.stdout
    if options.dict:
        out.write('mesh = {\n')
    elif options.json:
        out.write('[\n')


def write_trailer(options, out=None):
    if out is None:
        out = sys.stdout
    if options.dict:
        out.write('}')
    elif options.json:
        out.write('\n]\n')


class Descriptor(object):
    """MeSH descriptor."""

    def __init__(self, id_, name, scope, treenums):
        self.id = id_
        self.name = name
        self.scope = scope
        self.treenums = treenums

    def to_dict(self, no_id=False):
        d = {
            'name': self.name,
            'scope': self.scope,
            'treenums': self.treenums
        }
        if not no_id:
            d.update({'_id': self.id })
        return d

    @classmethod
    def from_xml(cls, element):
        uid = find_only(element, 'DescriptorUI').text
        name_element = find_only(element, 'DescriptorName')
        name = find_only(name_element, 'String').text
        try:
            scope = find_first(element, './/ScopeNote').text
            scope = scope.rstrip()
        except Exception as e:
            info('no scope note for %s (%s)' % (uid, name))
            scope = ''
        try:
            tree_number_list_element = find_only(element, 'TreeNumberList')
            tree_numbers = tree_number_list_element.findall('TreeNumber')
        except Exception as e:
            warning('missing tree numbers for %s (%s)' % (uid, name))
            tree_numbers = []
        return cls(uid, name, scope, [e.text for e in tree_numbers])


def write_data(descriptor, options, out=None):
    if out is None:
        out = sys.stdout
    desc_dict = descriptor.to_dict(no_id=options.dict)
    sorted_kv = sorted((k, v) for k, v in desc_dict.items())
    if options.dict:
        print('    %s: {' % repr(descriptor.id), file=out)
        print('\n'.join([8 * ' ' + '%s: %s,' % (repr(k), repr(v))
                                 for k, v in sorted_kv]), file=out)
        print('    },', file=out)
    elif options.json:
        if not write_data.first:
            out.write(',\n')
        out.write(pretty_dumps(descriptor.to_dict()))
    elif options.brat_norm:
        print('\t'.join([
            descriptor.id,
            'name:Title:{}'.format(descriptor.name),
            'info:Scope:{}'.format(descriptor.scope),
        ]), file=out)
    else:
        print('\t'.join([str(v) for k, v in sorted_kv]), file=out)
    write_data.first = False
write_data.first = True


def process_stream(stream, options):
    for event, element in stream:
        if event != 'end' or element.tag != 'DescriptorRecord':
            continue
        write_data(Descriptor.from_xml(element), options)
        element.clear()


def process(path, options):
    if not path.endswith('.gz'):
        return process_stream(ET.iterparse(path), options)
    else:
        with gzip.GzipFile(path) as stream:
            return process_stream(ET.iterparse(stream), options)


def main(argv):
    args = argparser().parse_args(argv[1:])
    if sum(1 for f in ('json', 'dict', 'brat_norm') if getattr(args, f)) > 1:
        error('at most one of -j, -d and -b arguments allowed.')
        return 1

    write_header(args)
    if args.top:
        for uid, name, treenums in meshtop:
            write_data(Descriptor(uid, name, None, treenums), args)
    process(args.file, args)
    write_trailer(args)
        
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
