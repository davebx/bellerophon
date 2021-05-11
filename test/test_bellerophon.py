import argparse
import hashlib
import os
import pysam
import tempfile

from unittest import TestCase

from bellerophon import filter_reads, merge_bams, __version__, __description__


class TestBellerophon(TestCase):
    def test_bellerophon(self):
        parser = argparse.ArgumentParser(description='Filter chimeric reads.', epilog=__description__)
        parser.add_argument('--forward', '-f', dest='forward', action='store', required=True, help='SAM/BAM/CRAM file with the first set of reads.')
        parser.add_argument('--reverse', '-r', dest='reverse', action='store', required=True, help='SAM/BAM/CRAM file with the second set of reads.')
        parser.add_argument('--output', '-o', dest='output', action='store', required=True, help='Output BAM file for filtered and paired reads.')
        parser.add_argument('--quality', '-q', dest='quality', type=int, action='store', required=False, default=20, help='Minimum mapping quality.')
        parser.add_argument('--threads', '-t', dest='threads', type=int, action='store', required=False, default=1, help='Threads.')
        parser.add_argument('--log-level', '-l', dest='log_level', type=str, action='store', required=False, default='ERROR', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG'], help='Log level.')
        parser.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))
        input_forward_reads = os.path.abspath('test_data/test_1500_forward.bam')
        input_reverse_reads = os.path.abspath('test_data/test_1500_reverse.bam')
        test_file = os.path.abspath('test_data/test_1500_merged_reads.bam')
        output_tempfile = tempfile.NamedTemporaryFile(prefix='test_filtered_merged_', suffix='.bam', delete=False, dir=os.getcwd())
        output_tempfile.close()
        args = ['--forward', input_forward_reads, '--reverse', input_reverse_reads, '--output', os.path.abspath(output_tempfile.name), '--quality', '10', '--log-level', 'DEBUG']
        arguments = parser.parse_args(args)
        output_filtered_forward, output_filtered_reverse = filter_reads(arguments)
        merged_output = merge_bams(arguments, output_filtered_forward, output_filtered_reverse)
        self.assertEqual(merged_output, 0)
        save = pysam.set_verbosity(0)
        test_out_fh = pysam.AlignmentFile(arguments.output, 'r')
        test_cmp_fh = pysam.AlignmentFile(test_file, 'r')
        pysam.set_verbosity(save)
        # Compare each read individually since bellerophon.merge_bams()
        # adds a row to the @PG section of the SAM header, per the SAM specification,
        # making the checksums differ.
        for output_read, test_read in zip(test_out_fh, test_cmp_fh):
            self.assertEqual(output_read, test_read)
        test_out_fh.close()
        test_cmp_fh.close()
        os.unlink(output_tempfile.name)
        self.assertFalse(os.path.exists(output_tempfile.name))

    def _hash_files(self, out_file, cmp_file):
        with open(out_file, 'rb') as fh:
            out_hash = hashlib.sha1(fh.read()).hexdigest()
        with open(cmp_file, 'rb') as fh:
            cmp_hash = hashlib.sha1(fh.read()).hexdigest()
        return out_hash, cmp_hash
