import argparse
import os
import pysam
import re
import tempfile

__version__ = '1.0rc0'

def filter_reads(args):
    retval = []
    forward_ids = []
    reverse_ids = []
    save = pysam.set_verbosity(0)
    ffh = pysam.AlignmentFile(args.forward, 'r', threads=args.threads)
    rfh = pysam.AlignmentFile(args.reverse, 'r', threads=args.threads)
    pysam.set_verbosity(save)
    if ffh.header.references != rfh.header.references or ffh.header.lengths != rfh.header.lengths:
        print('Error: The input files do not have the same sequence names or lengths.')
        return 1
    for handle in [ffh, rfh]:
        print('Loading reads from %s...' % str(handle.filename))
        previous_read = None
        all_reads = []
        unmapped_reads = []
        five_reads = []
        three_reads = []
        mid_reads = []
        counter = 0
        written_reads = 0
        processed_reads = 0
        come_in_here = re.compile(r'^[0-9]*M')
        dear_boy = re.compile(r'.*M$')
        have_a_cigar = re.compile(r'^[0-9]*[HS].*M.*[HS]$') # You're gonna go far, you're gonna fly
        output_tempfile = tempfile.NamedTemporaryFile(prefix='filtered_', suffix='.bam', delete=False, dir=os.getcwd())
        retval.append(output_tempfile.name)
        output_tempfile.close()
        output_fh = pysam.AlignmentFile(output_tempfile.name, 'wb', header=handle.header)
        for read in handle:
            processed_reads += 1
            # If this is 1. Not the first read, and 2. Not the previous read again:
            if previous_read is not None and read.query_name != previous_read:
                # If we have more than one read in the current batch and one
                # read is on the 5´ side of a ligation junction.
                if counter in [1, 2] and len(five_reads) == 1:
                    # Serve it forth.
                    output_fh.write(five_reads[0])
                    written_reads += 1
                else:
                    # Get the most recent read, set the unmapped flag, and send it
                    # to the output file.
                    new_read = all_reads[0]
                    new_read.is_unmapped = 1
                    output_fh.write(new_read)
                    written_reads += 1
                # Reset these variables to their original values.
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
            # If the read is aligned - and has mapped reads at the end, or it is
            # aligned + and has mapped reads at the beginning, it goes in the 5´
            # bin and is retained.
            elif (read.is_reverse and dear_boy.match(read.cigarstring) is not None) or (not read.is_reverse and come_in_here.match(read.cigarstring) is not None):
                five_reads.append(read)
            # If the read is aligned + and has mapped reads at the end, or it is
            # aligned - and has mapped reads at the beginning, it goes in the 3´
            # bin and is discarded.
            elif (read.is_reverse and come_in_here.match(read.cigarstring) is not None) or (not read.is_reverse and dear_boy.match(read.cigarstring) is not None):
                three_reads.append(read)
            # If it has mapped reads in the middle, put it in that list.
            elif have_a_cigar.match(read.cigarstring):
                mid_reads.append(read)
        # If we have a read.
        if counter == 1:
            # And it is on the 5´ side of a ligation junction
            if len(five_reads) == 1:
                # We send it to the output
                output_fh.write(five_reads[0])
            else:
                # Otherwise we flag it unmapped and push it out.
                new_read = all_reads[0]
                new_read.is_unmapped = 1
                output_fh.write(new_read)
                written_reads += 1
        # Or if we have two reads and one of them is on the 5´ side of a junction.
        elif counter == 2 and len(five_reads) == 1:
            # We do.
            output_fh.write(five_reads[0])
            written_reads += 1
        else:
            # The same kind of thing.
            new_read = all_reads[0]
            new_read.is_unmapped = 1
            output_fh.write(new_read)
            written_reads += 1
        print('Processed %d reads and output %d.' % (processed_reads, written_reads))
    # Send the filenames of the filtered alignments back to the caller.
    return retval

def merge_bams(args, filtered_forward, filtered_reverse):
    save = pysam.set_verbosity(0)
    forward = pysam.AlignmentFile(filtered_forward, 'r', threads=args.threads)
    reverse = pysam.AlignmentFile(filtered_reverse, 'r', threads=args.threads)
    pysam.set_verbosity(save)
    output_fh = pysam.AlignmentFile(args.output, 'wb', header=forward.header)
    processed_reads = 0
    skipped_reads = 0
    for forward_read, reverse_read in zip(forward, reverse):
        proper_pairs = 0
        if forward_read.query_name != reverse_read.query_name:
            read_tups = forward_read.query_name, reverse_read.query_name
            skipped_reads += 1
            continue
        if (forward_read.is_unmapped or forward_read.mapping_quality < args.quality) or (reverse_read.is_unmapped or reverse_read.mapping_quality < args.quality):
            continue
        if not forward_read.is_unmapped or reverse_read.is_unmapped:
            proper_pairs = 1
            if forward_read.reference_id == reverse_read.reference_id:
                distance = abs(forward_read.reference_start - reverse_read.reference_start)
                if forward_read.reference_start >= reverse_read.reference_start:
                    forward_length = -1 * distance
                    reverse_length = distance
                else:
                    forward_length = distance
                    reverse_length = -1 * distance
            else:
                forward_length = 0
                reverse_length = 0

        else:
            proper_pairs = 0
            forward_length = 0
            reverse_length = 0
        # Update some of the AlignedSegment attributes.
        forward_read.is_secondary = 0
        reverse_read.is_secondary = 0
        forward_read.is_unmapped = 0
        reverse_read.is_unmapped = 0
        forward_read.is_reverse = 0
        reverse_read.is_reverse = 0
        forward_read.is_supplementary = 0
        reverse_read.is_supplementary = 0
        forward_read.is_read1 = 1
        reverse_read.is_read2 = 1
        reverse_read.is_read1 = 0
        forward_read.is_read2 = 0
        reverse_read.is_unmapped = forward_read.is_unmapped
        forward_read.is_unmapped = reverse_read.is_unmapped
        forward_read.mate_is_unmapped = forward_read.is_unmapped
        reverse_read.mate_is_unmapped = reverse_read.is_unmapped
        forward_read.is_reverse = 0
        reverse_read.is_reverse = 1
        forward_read.mate_is_reverse = 1
        reverse_read.mate_is_reverse = 0
        forward_read.is_proper_pair = proper_pairs
        reverse_read.is_proper_pair = proper_pairs
        forward_read.is_paired = 1
        reverse_read.is_paired = 1
        forward_read.next_reference_start = reverse_read.reference_start
        reverse_read.next_reference_start = forward_read.reference_start
        forward_read.next_reference_name = reverse_read.reference_name
        reverse_read.next_reference_name = forward_read.reference_name
        forward_read.template_length = forward_length
        reverse_read.template_length = reverse_length
        output_fh.write(forward_read)
        output_fh.write(reverse_read)
        processed_reads += 1
    print('Processed and merged %d read pairs, skipped %d with mismatched read names.' % (processed_reads, skipped_reads))
    for filename in [filtered_forward, filtered_reverse]:
        os.unlink(filename)
    return 0