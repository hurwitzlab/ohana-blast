#!/usr/bin/env python3
"""
This script inserts a header row into a BLAST output file. The header looks like this:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

usage:
    python3 inserthdr.py test_HOT224_1_0025m.fa-contigs.tab
"""
import argparse
import shutil


def inserthdr(target_fp):

    target_with_header_fp = target_fp + '.add-hdr'
    with open(target_fp, 'rt') as input_file:
        if input_file.readline().startswith('qseqid'):
            print('file "{}" already has a header row'.format(target_fp))
        else:
            input_file.seek(0)
            with open(target_with_header_fp, 'wt') as output_file:

                output_file.write(
                    'qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n')
                for line in input_file:
                    output_file.write(line)

    # this function will not copy file metadata
    shutil.move(target_with_header_fp, target_fp)


def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('target_fp', metavar='FILE', help='path to file to which BLAST header will be added')
    args = arg_parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    inserthdr(**args.__dict__)
