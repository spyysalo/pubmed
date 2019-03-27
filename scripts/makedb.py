#!/usr/bin/env python3

import os
import sys
import tarfile

from logging import warning, error

try:
    import sqlitedict
except ImportError:
    error('failed to import sqlitedict; try `pip3 install sqlitedict`')
    raise


def argparser():
    from argparse import ArgumentParser
    ap = ArgumentParser()
    ap.add_argument('-s', '--suffix', default=None,
                    help='suffix of files to insert (default any)')
    ap.add_argument('-p', '--keep-path', default=False, action='store_true',
                    help='keep paths as part of db keys')
    ap.add_argument('db', help='database name')
    ap.add_argument('path', nargs='+', help='file or dir to insert into DB')
    return ap


def is_tar_gzip(path):
    return path.endswith('.gz') or path.endswith('.tar.gz')


def get_key(name, options):
    if options.keep_path:
        return name
    else:
        return os.path.basename(name)
    

def process_tgz(db, path, options):
    insert_count = 0
    print('Processing {}'.format(path), file=sys.stderr)
    with tarfile.open(path, 'r:gz') as tar:
        for m in tar.getmembers():
            suffix = os.path.splitext(m.name)
            if options.suffix is not None and options.suffix != suffix:
                continue
            f = tar.extractfile(m)
            if f is None:
                warning('failed to extract {} from {}'.format(m.name, path))
                continue
            content = f.read().decode('utf-8')
            key = get_key(m.name, options)
            db[key] = content
            insert_count += 1
            if insert_count % 1000 == 0:
                print('Inserted {} ...'.format(insert_count), end='\r',
                      file=sys.stderr, flush=True)
    print('Done, inserted {}, committing...'.format(insert_count),
          end='', file=sys.stderr, flush=True)
    db.commit()
    print('done.', file=sys.stderr)
    return insert_count
    

def process_file(db, path, options):
    with open(path) as f:
        content = f.read()
        key = get_key(path, options)
        db[key] = content
    print('Inserted {}, committing...'.format(key),
          end='', file=sys.stderr, flush=True)
    db.commit()
    print('done.', file=sys.stderr)
    return 1
        

def process(dbname, path, options):
    with sqlitedict.SqliteDict(dbname, autocommit=True) as db:
        if is_tar_gzip(path):
            return process_tgz(db, path, options)
        elif os.path.isfile(path):
            return process_file(db, path, options)
        elif os.path.isdir(path):
            raise NotImplementedError()
    


def main(argv):
    args = argparser().parse_args(argv[1:])
    total = 0
    for path in args.path:
        total += process(args.db, path, args)
    print('Finished, inserted {}.'.format(total), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
