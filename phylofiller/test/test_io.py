from unittest import TestCase, main

from skbio.util import get_data_path

from phylofiller.io import get_assembly_checksum


class IOTests(TestCase):
    def setUp(self):
        self.fp_assembly = get_data_path('assembly.fasta')
        self.exp_md5 = '1b89479f6bbae9a684b75b71e0cb492e'

    def tearDown(self):
        pass

    def test_get_assembly_checksum(self):
        obs_md5 = get_assembly_checksum(self.fp_assembly)
        self.assertEqual(
            self.exp_md5, obs_md5,
            "md5 for assembly differs from expectation")

        obs_md5 = get_assembly_checksum(get_data_path(
            'assembly_smallchars.fasta'))
        self.assertEqual(
            self.exp_md5, obs_md5,
            "Nucleotids in seq 5 are now lower case, "
            "should NOT result in different MD5")

        obs_md5 = get_assembly_checksum(get_data_path(
            'assembly_switchOrder.fasta'))
        self.assertEqual(
            self.exp_md5, obs_md5,
            "seqs order changed, should NOT result in different MD5")

        obs_md5 = get_assembly_checksum(get_data_path(
            'assembly_headerChanged.fasta'))
        self.assertEqual(
            self.exp_md5, obs_md5,
            "header changed, should NOT result in different MD5")

        obs_md5 = get_assembly_checksum(get_data_path(
            'assembly_layout.fasta'))
        self.assertEqual(
            self.exp_md5, obs_md5,
            "sequence layout changed, should NOT result in different MD5")

        obs_md5 = get_assembly_checksum(get_data_path(
            'assembly_seqdel.fasta'))
        self.assertNotEqual(
            self.exp_md5, obs_md5,
            "deletion of one nucleotide, should result in different MD5")

        obs_md5 = get_assembly_checksum(get_data_path(
            'assembly_newSeq.fasta'))
        self.assertNotEqual(
            self.exp_md5, obs_md5,
            "one additional sequence, should result in different MD5")


if __name__ == '__main__':
    main()
