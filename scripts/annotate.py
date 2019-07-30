#!/usr/bin/env python3

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>

import argparse
import os
import re
import sqlite3
import sys


def get_args():
    parser = argparse.ArgumentParser(description='Annotate BLAST for muSCOPE')

    parser.add_argument(
        '-b',
        '--blast_out',
        help='BLAST out file',
        type=str,
        metavar='FILE',
        required=True)

    parser.add_argument(
        '-a',
        '--annot_dir',
        help='Annotation directory',
        type=str,
        metavar='FILE',
        default='/work/03137/kyclark/ohana/sqlite')

    parser.add_argument(
        '-o',
        '--out_dir',
        help='Output directory',
        type=str,
        metavar='DIR',
        default='blast-annotated')

    parser.add_argument(
        '-v', '--verbose', help='Say more stuff', action='store_true')

    return parser.parse_args()


def main():
    args = get_args()
    out_dir = args.out_dir
    blast_out = args.blast_out
    annot_dir = args.annot_dir
    verbose = args.verbose

    if not os.path.isfile(blast_out):
        print('--blast_out file "{}" is not a file'.format(blast_out))
        exit(1)

    if not os.path.isdir(annot_dir):
        print('--annot_dir "{}" is not valid'.format(annot_dir))
        exit(1)

    dbs = list(filter(lambda x: x.endswith('.db'), os.listdir(annot_dir)))
    if len(dbs) == 0:
        print('Cannot find SQLite dbs in annot_dir "{}"'.format(annot_dir))
        exit(1)

    dbhs = dict()
    for db in dbs:
        base, ext = os.path.splitext(db)
        dbhs[base] = sqlite3.connect(os.path.join(annot_dir, db))

    os.makedirs(out_dir, exist_ok=True)

    blast_fields = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    gene_fields = [
        'gene_id', 'gene_name', 'cog_id', 'source', 'evalue', 'desc',
        'cog_categories'
    ]
    sql = 'select ' + ', '.join(gene_fields) + ' from gene where gene_name=?'

    # print headers for output
    out_file = os.path.join(out_dir, os.path.basename(blast_out))
    out_fh = open(out_file, 'wt')
    out_fh.write('\t'.join(['qseqid', 'sample'] + gene_fields) + '\n')

    def err(msg):
        if verbose:
            sys.stderr.write(msg + '\n')

    # BLAST output will contain gene ids like "HOT233_1_0770m_c4_1"
    # The sample name here would be "HOT233_1_0770m"
    # But we may have no annotations for that sample, so skip
    with open(blast_out, 'rt') as fh:
        for i, line in enumerate(fh):
            rec = dict(zip(blast_fields, line.rstrip().split('\t')))
            seqid = rec['sseqid']
            match = re.match('^(HOT\d{3}_(?:\d*[a-z]?_)?\d*m)', seqid)

            if not match:
                err('Failed to extract sample name from seqid "{}"'.format(
                    seqid))
                continue

            sample = match.group(0)
            if not sample in dbhs:
                err('Unknown sample "{}"'.format(sample))
                continue

            dbh = dbhs[sample]
            for row in dbh.execute(sql, (seqid, )):
                out_fh.write('\t'.join([rec['qseqid'], sample] +
                                       list(map(str, row))) + '\n')

    out_fh.close()
    print('Done, see output file "{}"'.format(out_file))


if __name__ == '__main__':
    main()
