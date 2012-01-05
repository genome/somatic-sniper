#!/usr/bin/env python

from integrationtest import IntegrationTest, main
import unittest

class TestSniper(IntegrationTest, unittest.TestCase):

    def test_sniper(self):
        normal_bam = self.inputFiles("n-small.bam")[0]
        tumor_bam = self.inputFiles("t-small.bam")[0]
        reference_fasta = self.inputFiles("small.fa")[0]
        expected_file = self.inputFiles("expected.vcf")[0]
        output_file = self.tempFile("output.vcf")
        params = [
            "-F vcf", "-f", reference_fasta, tumor_bam, normal_bam, output_file
        ]
        rv, err = self.execute(params)
        self.assertEqual(0, rv)
        self.assertFilesEqual(expected_file, output_file, filter_regex="##fileDate|##reference=")

if __name__ == "__main__":
    main()
