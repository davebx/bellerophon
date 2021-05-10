#!/usr/bin/env python
import argparse
from bellerophon import filter_reads, merge_bams, __version__, __description__


def main():
    parser = argparse.ArgumentParser(description='Filter chimeric reads.', epilog=__description__)
    parser.add_argument('--forward', '-f', dest='forward', action='store', required=True, help='SAM/BAM/CRAM file with the first set of reads.')
    parser.add_argument('--reverse', '-r', dest='reverse', action='store', required=True, help='SAM/BAM/CRAM file with the second set of reads.')
    parser.add_argument('--output', '-o', dest='output', action='store', required=True, help='Output BAM file for filtered and paired reads.')
    parser.add_argument('--quality', '-q', dest='quality', type=int, action='store', required=False, default=20, help='Minimum mapping quality.')
    parser.add_argument('--threads', '-t', dest='threads', type=int, action='store', required=False, default=1, help='Threads.')
    parser.add_argument('--log-level', '-l', dest='log_level', type=str, action='store', required=False, default='ERROR', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG'], help='Log level.')
    parser.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))
    args = parser.parse_args()
    filtered_forward, filtered_reverse = filter_reads(args)
    return merge_bams(args, filtered_forward, filtered_reverse)


if __name__ == '__main__':
    exit(main())
