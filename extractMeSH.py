#!/usr/bin/env python

# Extract unique ID, name and tree numbers from MeSH XML data.

import sys

import xml.etree.ElementTree as ET

from logging import warn

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
    ap.add_argument('-d', '--dict', default=False, action='store_true',
                    help='Output python dictionary (default TSV)')
    ap.add_argument("file", metavar="FILE", help="Input MeSH XML.")
    return ap
    
def find_only(element, match):
    """Return the only matching child of the given element.

    Fail on assert if there are no or multiple matches.
    """
    found = element.findall(match)
    assert len(found) == 1, 'Error: expected 1 %s, got %d' % (match, len(found))
    return found[0]

def write_header(options, out=None):
    if out is None:
        out = sys.stdout
    if options.dict:
        out.write('mesh = {\n')

def write_trailer(options, out=None):
    if out is None:
        out = sys.stdout
    if options.dict:
        out.write('}')

def get_uid_name_treenums(descriptor):
    """Return text for descriptor UID, name and tree numbers."""
    uid_element = find_only(descriptor, 'DescriptorUI')
    name_element = find_only(descriptor, 'DescriptorName')
    name_str_element = find_only(name_element, 'String')
    try:
        tree_number_list_element = find_only(descriptor, 'TreeNumberList')
        tree_numbers = tree_number_list_element.findall('TreeNumber')
    except Exception, e:
        warn('Failed to get tree numbers for %s (%s)' %
             (uid_element.text, name_str_element.text))
        tree_numbers = []
    return (uid_element.text,
            name_str_element.text,
            [e.text for e in tree_numbers])

def write_data(uid, name, treenums, options, out=None):
    if out is None:
        out = sys.stdout
    if not options.dict:
        print >> out, '\t'.join([uid, name] + treenums)
    else:
        print >> out, '''    %s: {
        'name': %s,
        'tnum': %s
    },''' % (repr(uid), repr(name), repr(treenums))
        
def main(argv):
    args = argparser().parse_args(argv[1:])

    write_header(args)
    if args.top:
        for uid, name, treenums in meshtop:
            write_data(uid, name, treenums, args)
    for event, element in ET.iterparse(args.file):
        if event != 'end' or element.tag != 'DescriptorRecord':
            continue
        uid, name, treenums = get_uid_name_treenums(element)
        write_data(uid, name, treenums, args)        
        element.clear()
    write_trailer(args)
        
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
