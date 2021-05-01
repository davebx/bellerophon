#!/usr/bin/env python
import argparse
import os
from bellerophon import filter_reads, merge_bams, __version__

def main():
    long_description = 'Filter two single-end BAM, SAM, or CRAM files for reads where there is high-quality mapping on both sides of a ligation junction, retaining the 5´ side of that mapping, then merge them into one paired-end BAM file. '
    parser = argparse.ArgumentParser(description='Filter chimeric reads.', epilog=long_description)
    parser.add_argument('--forward', '-f', dest='forward', action='store', required=True, help='SAM/BAM/CRAM file with the first set of reads.')
    parser.add_argument('--reverse', '-r', dest='reverse', action='store', required=True, help='SAM/BAM/CRAM file with the second set of reads.')
    parser.add_argument('--output', '-o', dest='outfile', action='store', required=True, help='Output BAM file for filtered and paired reads.')
    parser.add_argument('--quality', '-q', dest='mq', type=int, action='store', required=True, help='Minimum mapping quality.')
    parser.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))
    args = parser.parse_args()

    filtered_forward, filtered_reverse = filter_reads(args)
    # print(filtered_forward, filtered_reverse)
    return merge_bams(filtered_forward, filtered_reverse, args.outfile, args.mq)

if __name__ == '__main__':
    exit(main())
