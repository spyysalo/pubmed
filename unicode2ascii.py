#!/usr/bin/env python

# Replaces Unicode characters in input text with ASCII
# approximations based on file with mappings between the two.

import sys
import os
import codecs
import re

from itertools import tee, chain

verbose = True

# The name of the file from which to read the replacement. Each line
# should contain the hex code for the unicode character, TAB, and
# the replacement string.

MAPPING_FILE_NAME = "entities.dat"

# For statistics and summary of missing mappings in verbose mode
map_count = {}
missing_mapping = {}


# Support wide unichr on narrow python builds. From @marcovzla, see
# https://github.com/spyysalo/nxml2txt/pull/4.
def wide_unichr(i):
    try:
        return chr(i)
    except ValueError:
        return (r'\U' + hex(i)[2:].zfill(8)).decode('unicode-escape')


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ... (sN, None)"
    a, b = tee(iterable)
    next(b, None)
    return zip(a, chain(b, (None, )))


def read_mapping(f, fn="mapping data"):
    """
    Reads in mapping from Unicode to ASCII from the given input stream
    and returns a dictionary keyed by Unicode characters with the
    corresponding ASCII characters as values. The expected mapping
    format defines a single mapping per line, each with the format
    CODE\tASC where CODE is the Unicode code point as a hex number and
    ASC is the replacement ASCII string ("\t" is the literal tab
    character). Any lines beginning with "#" are skipped as comments.
    """

    # read in the replacement data
    linere = re.compile(r'^([0-9A-Za-z]{4,})\t(.*)$')
    mapping = {}

    for i, l in enumerate(f):
        # ignore lines starting with "#" as comments
        if len(l) != 0 and l[0] == "#":
            continue

        m = linere.match(l)
        assert m, "Format error in %s line %s: '%s'" % (
            fn, i+1, l.replace("\n","").encode("utf-8"))
        c, r = m.groups()

        c = wide_unichr(int(c, 16))
        assert c not in mapping or mapping[c] == r, "ERROR: conflicting mappings for %.4X: '%s' and '%s'" % (ord(c), mapping[c], r)

        # exceptions: literal '\n' maps to newline, initial and terminal
        # '\b' map to backspace (see map_character)
        if r == '\\n':
            r = '\n'
        if r[:2] == '\\b':
            r = '\b'+r[2:]
        if r[-2:] == '\\b':
            r = r[:-2]+'\b'

        mapping[c] = r

    return mapping


word_char = re.compile(r'^\w$', re.U)


def map_character(prev, char, next_, mapping):
    """
    Maps character char to mapping where prev and next_ are the characters
    preceding and following c, handling boundary space insertion.
    """
    if not mapping:
        return mapping
    # Boundary space insertion: if mapping begins (or ends) with '\b',
    # that character is replaced with a space if the previous (next)
    # character is a word character and removed otherwise.
    if mapping[0] == '\b':
        space = ' ' if prev and word_char.match(prev) else ''
        mapping = space + mapping[1:]
    if mapping and mapping[-1] == '\b':
        space = ' ' if next_ and word_char.match(next_) else ''
        mapping = mapping[:-1] + space
    return mapping


def process(f, out, mapping):
    """
    Applies the given mapping to replace characters other than 7-bit
    ASCII from the given input stream f, writing the mapped text to
    the given output stream out.
    """
    global map_count, missing_mapping

    missing_count = 0
    for line in f:
        prev = None
        for c, next_ in pairwise(line):
            curr = c
            if ord(c) >= 128:
                # higher than 7-bit ASCII, might wish to map
                if c in mapping:
                    map_count[c] = map_count.get(c,0)+1
                    c = map_character(prev, c, next_, mapping[c])
                else:
                    missing_mapping[c] = missing_mapping.get(c,0)+1
                    # escape into numeric Unicode codepoint
                    c = "<%.4X>" % ord(c)
                    missing_count += 1
            out.write(c)
            prev = curr
    return missing_count


def print_summary(out, mapping):
    """
    Prints human-readable summary of statistics and missing mappings
    for the input into the given output stream.
    """

    global map_count, missing_mapping

    print("Characters replaced       \t%d" % sum(map_count.values()), file=out)
    for c in sorted(map_count.keys(), key=lambda k: map_count[k]):
        try:
            print("\t%.4X\t%s\t'%s'\t%d" % (ord(c), c.encode("utf-8"),
                                            mapping[c], map_count[c]), file=out)
        except:
            print("\t%.4X\t'%s'\t%d" % (ord(c), mapping[c], map_count[c]),
                  file=out)
    print("Characters without mapping\t%d" % sum(missing_mapping.values()),
          file=out)
    for c in sorted(missing_mapping.keys(), key=lambda k: missing_mapping[k]):
        try:
            print("\t%.4X\t%s\t%d" % (ord(c), c.encode("utf-8"),
                                      missing_mapping[c]), file=out)
        except:
            print("\t%.4X\t?\t%d" % (ord(c), missing_mapping[c]), file=out)


def argparser():
    """
    Returns an argument parser for the script.
    """
    import argparse
    ap=argparse.ArgumentParser(description="Replaces Unicode characters in input text with ASCII approximations.")
    ap.add_argument('-d', '--directory', default=None, help="Directory for output (stdout by default)")
    ap.add_argument('-v', '--verbose', default=False, action='store_true', help="Verbose output")
    ap.add_argument('file', nargs='+', help='Input text file')
    return ap


def main(argv):
    global options

    # argument processing
    options = argparser().parse_args(argv[1:])

    # read in mapping
    try:
        mapfn = MAPPING_FILE_NAME

        if not os.path.exists(mapfn):
            # fall back to trying in script dir
            mapfn = os.path.join(os.path.dirname(__file__), 
                                 os.path.basename(MAPPING_FILE_NAME))

        with codecs.open(mapfn, encoding="utf-8") as f:
            mapping = read_mapping(f, mapfn)
    except IOError as e:
        print("Error reading mapping from %s: %s" % (MAPPING_FILE_NAME, e),
              file=sys.stderr)
        return 1

    # primary processing
    for fn in options.file:
        try:
            if fn == '-':
                fn = '/dev/stdin' # TODO: make portable
            with codecs.open(fn, encoding="utf-8") as f:
                if options.directory is None:
                    process(f, sys.stdout, mapping)
                else:
                    bfn = os.path.basename(fn)
                    ofn = os.path.join(options.directory, bfn)
                    with codecs.open(ofn, 'wt', encoding="utf-8") as out:
                        process(f, out, mapping)
        except IOError as e:
            print("Error processing %s: %s" % (fn, e), file=sys.stderr)

    # optionally print summary of mappings
    if options.verbose:
        print_summary(sys.stderr, mapping)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
