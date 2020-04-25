from unittest import TestCase, main

from skbio.util import get_data_path

from phylofiller.io import get_assembly_checksum, read_metadata


class IOTests(TestCase):
    def setUp(self):
        self.fp_assembly = get_data_path('assembly.fasta')
        self.exp_md5 = '1b89479f6bbae9a684b75b71e0cb492e'
        self.fp_metadata = get_data_path('meta.tsv')

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

    def test_read_metadata(self):
        with self.assertRaisesRegex(ValueError, 'Conflicting column name\(s\) in your metadata: "__line_number" please rename!'):
            read_metadata(get_data_path("meta_sysColName.tsv"))

        with self.assertRaisesRegex(ValueError, 'Header of first column must be "organism"'):
            read_metadata(get_data_path("meta_wrongIndexName.tsv"))

        with self.assertRaisesRegex(ValueError, '2 organisms have been re-defined in your metadata. Please fix!'):
            read_metadata(get_data_path("meta_ambigOrganisms.tsv"))

        # config = {'projects': {'fungi'}}
        # with self.assertRaisesRegex(ValueError, "No 'assemblies' defined in configuration for project 'fungi'"):
        #     validate_input_configuration(config)


if __name__ == '__main__':
    main()
