import argparse
import hashlib
import os
import pysam
import tempfile

from unittest import TestCase

from bellerophon import filter_reads, merge_bams, __version__, __description__


class TestBellerophon(TestCase):
    def setUp(self):
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
        self.test_filtered_forward_reads = os.path.abspath('test_data/test_1500_filtered_forward.bam')
        self.test_filtered_reverse_reads = os.path.abspath('test_data/test_1500_filtered_reverse.bam')
        self.output_tempfile = tempfile.NamedTemporaryFile(prefix='test_filtered_merged_', suffix='.bam', delete=False, dir=os.getcwd())
        args = ['--forward', input_forward_reads, '--reverse', input_reverse_reads, '--output', self.output_tempfile.name, '--quality', '10', '--log-level', 'DEBUG']
        self.arguments = parser.parse_args(args)

    def test_filter(self):
        output_hash = None
        test_hash = None
        output_filtered_forward, output_filtered_reverse = filter_reads(self.arguments)
        self.assertTrue(os.path.exists(output_filtered_forward))
        self.assertTrue(os.path.exists(output_filtered_reverse))
        # The output files should be identical including headers after filtering.
        for output_file, test_file in zip([output_filtered_forward, output_filtered_reverse], [self.test_filtered_forward_reads, self.test_filtered_reverse_reads]):
            output_hash, test_hash = self._hash_files(output_file, test_file)
        self.assertTrue(output_hash == test_hash, 'SHA1 checksum mismatch in test output:\n[OUT]\t%s:\t"%s" !=\n[CMP]\t%s:\t"%s"' % (output_file, output_hash, test_file, test_hash))
        for temporary_file in [output_filtered_forward, output_filtered_reverse]:
            os.unlink(temporary_file)

    def test_merge(self):
        test_file = os.path.abspath('test_data/test_1500_merged_reads.bam')
        output_filtered_forward, output_filtered_reverse = filter_reads(self.arguments)
        merged_output = merge_bams(self.arguments, output_filtered_forward, output_filtered_reverse)
        self.assertEqual(merged_output, 0)
        save = pysam.set_verbosity(0)
        test_out_fh = pysam.AlignmentFile(self.arguments.output, 'r')
        test_cmp_fh = pysam.AlignmentFile(test_file, 'r')
        pysam.set_verbosity(save)
        # Compare each read individually since bellerophon.merge_bams()
        # adds a row to the @PG section of the SAM header, making the checksums differ.
        for output_read, test_read in zip(test_out_fh, test_cmp_fh):
            self.assertEqual(output_read, test_read)
        os.unlink(self.arguments.output)

    def _hash_files(self, out_file, cmp_file):
        with open(out_file, 'rb') as fh:
            out_hash = hashlib.sha1(fh.read()).hexdigest()
        with open(cmp_file, 'rb') as fh:
            cmp_hash = hashlib.sha1(fh.read()).hexdigest()
        return out_hash, cmp_hash
