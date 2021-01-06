#!/usr/bin/env python3
'''
@file: text_gio.py
@auth: Sprax Lines
@date: 2020.11.22
DNA pattern matching functions and some text file utilities.
Written with Python version >= 3.8.5
'''
import argparse
import errno
# import fnmatch
import glob
import ntpath
import os
import os.path
import pickle
import random
# import re
import sys
import time
from pdb import set_trace
from typing import Deque, Dict, List, Set, Tuple

# from pprint import pprint


def cwd():
    '''current working directory'''
    return os.path.dirname(os.path.realpath('.'))


def path_base(path):
    '''
    Returns only the basename (no parent path)
    for Unix or Windows style paths.
    Logic:  'C:\\tmp/some\\file.txt' => 'file'
    '''
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)


def path_name(path):
    '''
    Returns only the basename (no parent path, no file extension)
    for Unix or Windows style paths.
    Logic:  'C:\\tmp/some\\file.txt' => 'file'
    '''
    return os.path.splitext(path_base(path))[0]


def get_abs_path(path):
    '''
    Convert the specified path to an absolute path (if it isn't already one).
    Returns the corresponding absolute path.
    '''
    return os.path.abspath(path)


def make_abs_path(dirpath, filepath):
    '''
    Given a directory path and a filename or relative file path,
    get the absolute path for the specified file under that directory.
    Returns this absolute path as a string suitable as an argument to open().
    '''
    return os.path.abspath(os.path.join(dirpath, filepath))


def print_stdout_stderr(text):
    ''' print text to stdout and stderr '''
    print("sys.stdout: ", text, file=sys.stdout)
    print("sys.stderr: ", text, file=sys.stderr)


def open_out_file(file_spec, label='text'):
    '''
    returns a file handle open for writing,
    to be closed by the caller, else None
    '''
    if file_spec:
        if file_spec in ['-', 'stdout']:
            return sys.stdout
        else:
            try:
                out_file = open(file_spec, 'w')
            except IOError as ex:
                if ex.errno != errno.ENOENT:
                    raise
                print("IOError opening {} file [{}]:".format(
                    label, file_spec), ex)
                out_file = sys.stdout
            return out_file
    else:
        return None


def read_lines(file_spec, charset='utf8'):
    '''read and yield all lines of a text file as a iter of str'''
    with open(file_spec, 'r', encoding=charset) as text:
        for line in text:
            yield line.rstrip()


def read_text_lines(file_spec, charset='utf8'):
    '''read and yield non-empty lines of a text file as a iter of str'''
    with open(file_spec, 'r', encoding=charset) as text:
        for line in text:
            line = line.strip()
            if line:
                yield line


def pickle_file(in_path, out_path, data_struct, data_adder, charset='utf8'):
    '''read in_file into data_struct via data_adder then save to out_path'''
    lines_in = 0
    lines = read_text_lines(in_path, charset)
    for line in lines:
        data_adder(data_struct, line)
        lines_in += 1
    with open(out_path, 'wb') as out_file:
        pickle.dump(data_struct, out_file)
    return (lines_in, len(data_struct))


def pickle_word_list(in_path, out_path, word_set=None, adder=set.add, charset='utf8'):
    '''
    read single words/strings per line from in_file
    and save them as a set to out_path as a pickle file
    '''
    if word_set is None:
        word_set = set()
    return pickle_file(in_path, out_path, word_set, adder, charset)


def read_file(file_spec, charset='us-ascii'):   # Not: 'utf8'
    '''
    read and return all contents of file as one str
    '''
    with open(file_spec, 'r', encoding=charset) as src:
        return src.read()


def read_file_eafp(file_spec, charset='us-ascii'):   # Not: 'utf8'
    '''
    Read contents of file_spec.
    Easier to Ask for Forgiveness than ask Permission.
    '''
    try:
        src = open(file_spec, 'r', encoding=charset)
    except IOError as ex:
        if ex.errno != errno.ENOENT:
            raise
        print("WARNING: {} does not exist".format(file_spec))
        return None
    else:
        text = src.read()
        src.close()
        return text


def read_text_file(file_spec):
    '''
    Read and return all contents of file as one str
    Try to read ascii or utf-8 and failover to iso-8859-1, etc.
    '''
    try:
        return read_file(file_spec, 'utf-8')
    except UnicodeDecodeError:
        return read_file(file_spec, 'iso-8859-1')


def glob_files(dir_path: str, end_pat: str, recursive=False, verbose=0):
    '''
    Returns list of (relative) paths of files maching file_pat in dir_path.
    Uses glob for *nix-like file name matching.
    Recursive is OFF by default.
    '''
    glob_pat = dir_path + "/*" + end_pat
    if verbose > 3:
        print("find_file_paths: glob_pat(%s)" % glob_pat)
    return filter(os.path.isfile, glob.iglob(glob_pat))


def scan_files(dir_path: str, end_pat: str) -> Dict:
    '''
    Returns list of file names (only) maching file_pat in dir_path.
    Uses os.scandir for Posix-like file name matching.
    Recursive is always ON.
    '''
    return [f.name for f in os.scandir(dir_path) if f.name.endswith(end_pat)]


def name_seq_map_from_dir(dna_dir: str,
                          end_pat: str,
                          charset: str = 'us-ascii') -> Dict:
    '''
    TODO: separate file finding from loading.
    '''
    name_acgt_map = {}
    glob_pat = dna_dir + "/*" + end_pat
    for path in filter(os.path.isfile, glob.iglob(glob_pat)):
        name = path_name(path)
        name_acgt_map[name] = read_file(path, charset)
    return name_acgt_map


def load_dna_map(dna_dir: str,
                 file_pat: str,
                 charset: str,
                 verbose: int = 1) -> Dict:
    '''
    Load NCBI-named Proteins as DNA sequences (acgt) from raw text files.
    Uses glob instead of scandir for the simplicity of one flat data dir,
    not a tree of subdirectories.
    '''
    dna_files = glob_files(dna_dir, file_pat, recursive=False, verbose=verbose)
    if verbose > 3:
        print("glob DNA files:", list(dna_files))
    return name_seq_map_from_dir(dna_dir, file_pat, charset)


def find_proteins(name_acgt_map: Dict,
                  acgt_str: str,
                  max_find: int = 1,
                  verbose: int = 0) -> int:
    '''
    Find up to max_find proteins matching the search string acgt_str.
    If max_find == 0, all proteins are searched.
    If max_find == 1, the proteins are searched in random order;
    Otherwise, the proteins are searched in sorted key order.
    Threaded searching means the found order may vary.
    By default, the names of matching proteins are re-sorted.

    NOTE: str.find uses The Fast Search Algorithm (aka “BMHBNFS”)
    as in Boyer-Moore-Horspool-..., expected to run in sub-linear
    time, and slightly faster than KMP, for packed strings with
    small alphabet / many repetitions.  The extra space for BMH
    is amortized.
    @see: http://effbot.org/zone/stringlib.htm

    One SQLite3 query similar to the use of str.find here might be:
        sqlite> .timer ON
        sqlite> SELECT ncbi_name, INSTR(acgt_text, 'tattatttttatat') - 1 AS idx
                FROM pfind_protein
                WHERE idx > 0
                LIMIT 0, $(max_find);
    when max_find > 0; if max_find == 0, omit the LIMIT part.

    Or, using the Django model with parameters passed into raw query:
        offset = 0
        INT_MAX = 2**31 - 1
        limit = max_find if max_find > 0 else INT_MAX
        Protein.objects.raw("SELECT ncbi_name,
                INSTR(acgt_text, %s) - 1 AS idx
                WHERE idx >= 0
                ORDER BY ncbi_name ASC
                LIMIT %s, %s", [acgt_str, offset, limit])

    Or omit the LIMIT part if max_find == 0.
    '''
    found = 0
    names = list(name_acgt_map.keys())
    if max_find == 1:
        random.shuffle(names)
    else:
        names.sort()

    for name in names:
        raw = name_acgt_map[name]
        idx = raw.find(acgt_str)
        if idx < 0:
            if verbose > 2:
                print("%s excludes %s" % (name, acgt_str))
        else:
            found += 1
            if verbose > 1:
                print("%s CONTAINS %s at %d" % (name, acgt_str, idx))
            if 0 < max_find <= found:
                break
    return found


def unit_test(opt):
    '''
    Test the functions above with one set of options.
    '''
    verbose = opt.verbose

    prot_map = load_dna_map(opt.dna_dir, opt.file_pat, opt.charset, verbose)

    beg_time = time.perf_counter()
    found_ct = find_proteins(prot_map, opt.acgt_str, opt.max_find, verbose)
    dur_time = time.perf_counter() - beg_time

    print("END unit_test: Found %d / %d matches for %s in %7.4g seconds."
          % (found_ct, opt.max_find, opt.acgt_str, dur_time))

###############################################################################


def main():
    '''
    Test driver for finding DNA fragments in sequenced proteins.
    '''
    default_acgt_str = "tattatttttatat"
    default_max_find = 1

    parser = argparse.ArgumentParser(
        # usage='%(prog)s [options]',
        description="Test driver for DNA alignment/Protein search")

    parser.add_argument('acgt_str', type=str, nargs='?',
                        default=default_acgt_str,
                        help=('DNA search string.  Example: catattaggaatttt. '
                              'Default: %s.' % default_acgt_str
                              ))

    parser.add_argument('max_find', type=int, nargs='?',
                        default=default_max_find,
                        help=('Maximum matches to find and show. '
                              ' 0 means no limit. '
                              ' Default: %d (stop at first match)'
                              % default_max_find
                              ))

    parser.add_argument('-c', '--charset', dest='charset', type=str,
                        default='us-ascii',
                        # default='iso-8859-1',
                        help='charset encoding of input text')
    parser.add_argument('-d', '--dna_dir', dest='dna_dir', type=str,
                        default='DNA',
                        help='path to dir containing DNA sequence files')
    parser.add_argument('-file_pat', type=str,
                        default='.raw',
                        help='Pattern matching DNA sequence file names')
    parser.add_argument('-name_order', action='store_true',
                        help=("Find and show matches in name-sorted order v. "
                              "random (maybe parallelized) order.  "
                              "NOT IMPLEMENTED"))

    parser.add_argument('-out_file', type=str, nargs='?', const='-',
                        help='output file for search results (default: None)')
    parser.add_argument('-repr', action='store_true',
                        help='output repr of data, not raw data')
    parser.add_argument('-verbose', type=int, nargs='?', const=1, default=1,
                        help='verbosity of output (default: 1)')
    args = parser.parse_args()

    if args.verbose > 5:
        print("outfile: <{}>".format(args.out_file))
        print("args:", args)
        print(__doc__)
        exit(0)

    unit_test(args)


if __name__ == '__main__':
    main()
