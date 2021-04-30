import argparse
import os
import pysam
import re
import tempfile


def filter_reads(args):
    retval = []
    for strand in [args.forward, args.reverse]:
        previous_read = None
        all_reads = []
        unmapped_reads = []
        five_reads = []
        three_reads = []
        mid_reads = []
        counter = 0
        come_in_here = re.compile(r'^[0-9]*M')
        dear_boy = re.compile(r'.*M$')
        have_a_cigar = re.compile(r'^[0-9]*[HS].*M.*[HS]$') # You're gonna go far, you're gonna fly
        infile = pysam.AlignmentFile(strand, 'r')
        outfile = tempfile.NamedTemporaryFile(prefix='filtered_', suffix='.bam', delete=False, dir=os.getcwd())
        retval.append(outfile.name)
        outfile.close()
        outfh = pysam.AlignmentFile(outfile.name, 'wb', header=infile.header)
        for read in infile:
            # If this is 1. Not the first read, and 2. Not the previous read again:
            if previous_read is not None and read.query_name != previous_read:
                if counter in [1, 2] and len(five_reads) == 1:
                    outfh.write(five_reads[0])
                else:
                    # Get the most recent read, set the unmapped flag, and send it
                    # to the output file.
                    new_read = all_reads[0]
                    new_read.is_unmapped = 1
                    outfh.write(new_read)
                counter = 0
                all_reads = []
                unmapped_reads = []
                five_reads = []
                three_reads = []
                mid_reads = []
            counter += 1
            all_reads.append(read)
            previous_read = read.query_name
            # Determine whether read is unmapped, or has mapped reads spanning a junction
            if read.is_unmapped:
                unmapped_reads.append(read)
            # If the read is aligned - and has mapped reads at the end, or it is +
            # and has mapped reads at the beginning, it goes in fives.
            elif (read.is_reverse and dear_boy.match(read.cigarstring) is not None) or (not read.is_reverse and come_in_here.match(read.cigarstring) is not None):
                five_reads.append(read)
            # If the read is aligned + and has mapped reads at the end, or it is -
            # and has mapped reads at the beginning, it goes in threes and gets discarded.
            elif (read.is_reverse and come_in_here.match(read.cigarstring) is not None) or (not read.is_reverse and dear_boy.match(read.cigarstring) is not None):
                three_reads.append(read)
            # If it has mapped reads in the middle, put it in that list.
            elif have_a_cigar.match(read.cigarstring):
                mid_reads.append(read)
        if counter == 1:
            if len(five_reads) == 1:
                outfh.write(five_reads[0])
            else:
                new_read = all_reads[0]
                new_read.is_unmapped = 1
                outfh.write(new_read)
        elif counter == 2 and len(five_reads) == 1:
            outfh.write(five_reads[0])
        else:
            new_read = all_reads[0]
            new_read.is_unmapped = 1
            outfh.write(new_read)
    return retval

def merge_bams(filtered_forward, filtered_reverse, output_file, quality):
    forward = pysam.AlignmentFile(filtered_forward, 'r')
    reverse = pysam.AlignmentFile(filtered_reverse, 'r')
    outfile = pysam.AlignmentFile(output_file, 'wb', header=forward.header)
    forward_check = str(forward.header)
    if str(forward.header) != str(reverse.header):
        print('Error: The input SAM headers do not match.')
        return 1
    input_line = 0
    for forward_read, reverse_read in zip(forward, reverse):
        proper_pairs = 0
        input_line += 1
        if input_line % 100000 == 0:
            print('Processed and merged %d of %d reads' % (input_line, len(forward)))
        if forward_read.query_name != reverse_read.query_name:
            print('Error at input line %d, forward read name %s does not match reverse read name %s.' % (input_line, forward_read.query_name, reverse_read.query_name))
            return 1
        if (forward_read.is_unmapped or forward_read.mapping_quality < quality) or (reverse_read.is_unmapped or reverse_read.mapping_quality < quality):
            continue
        if False not in [forward_read.is_unmapped, reverse_read.is_unmapped]:
            proper_pairs = 1
            if forward_read.reference_id == reverse_read.reference_id:
                distance = abs(forward_read.reference_start - reverse_read.reference_start)
                if forward_read.reference_start >= reverse_read.reference_start:
                    forward_distance = -1 * distance
                    reverse_distance = distance
                else:
                    forward_distance = distance
                    reverse_distance = -1 * distance
            else:
                forward_distance = 0
                reverse_distance = 0

        else:
            proper_pairs = 0
            forward_distance = 0
            reverse_distance = 0

        forward_read.flag = forward_read.flag ^ (1 << 8) | (1 << 7)
        forward_read.flag = forward_read.flag ^ (1 << 8) | (1 << 6)
        forward_bitmask = (forward_read.flag & (1 << 10) | forward_read.flag & (1 << 9) | \
                           reverse_read.flag & (1 << 4) | forward_read.flag & (1 << 4) | \
                           reverse_read.flag & (1 << 2) | forward_read.flag & (1 << 2) | \
                           (proper_pairs << 1) | (1 << 0))
        reverse_bitmask = (reverse_read.flag & (1 << 10) | reverse_read.flag & (1 << 9) | \
                           forward_read.flag & (1 << 4) | reverse_read.flag & (1 << 4) | \
                           forward_read.flag & (1 << 2) | reverse_read.flag & (1 << 2) | \
                           (proper_pairs << 1) | (1 << 0))
        forward_read.flag = forward_read.flag | forward_bitmask
        reverse_read.flag = reverse_read.flag | reverse_bitmask
        forward_start = forward_read.reference_start
        forward_read.next_reference_start = reverse_read.reference_start
        reverse_read.next_reference_start = forward_read.reference_start
        outfile.write(forward_read)
        outfile.write(reverse_read)
    for filename in [filtered_forward, filtered_reverse]:
        os.unlink(filename)
    return 0