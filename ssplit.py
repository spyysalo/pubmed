#!/usr/bin/env python

import re
import nltk.data

# rarely followed by a sentence split
NS_STRING = [
    'a.k.a.', 'approx.', 'ca.', 'cf.', 'e.g.', 'et al.', 'f.c.', 'i.e.',
    'lit.', 'vol.'
]

# rarely followed by a sentence split if the next character is a number
NS_NUM_STRING = [
    'fig.', 'ib.', 'no.',
]

NS_RE = re.compile(r'(?i)\b((?:' +
                   '|'.join(re.escape(s) for s in NS_STRING) +
                   ')\s*)\n')
NS_NUM_RE = re.compile(r'(?i)\b((?:' +
                       '|'.join(re.escape(s) for s in NS_NUM_STRING) +
                       ')\s*)\n(\d)')

nltk_splitter = nltk.data.load('tokenizers/punkt/english.pickle')

def ssplitter(s):
    split = '\n'.join(nltk_splitter.tokenize(s.strip()))
    split = NS_RE.sub(r'\1 ', split)
    split = NS_NUM_RE.sub(r'\1 \2', split)
    return split
