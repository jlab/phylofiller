from unittest import TestCase, main

from skbio.util import get_data_path

from phylofiller.metadata import read_metadata, get_augustus_reference_species


class MetadataTests(TestCase):
    def setUp(self):
        self.fp_metadata = get_data_path('meta.tsv')

    def tearDown(self):
        pass

    def test_read_metadata(self):
        with self.assertRaisesRegex(
                ValueError,
                'Conflicting column name\\(s\\) in your metadata: '
                '"__line_number" please rename!'):
            read_metadata(get_data_path("meta_sysColName.tsv"))

        with self.assertRaisesRegex(
                ValueError,
                'Header of first column must be "organism"'):
            read_metadata(get_data_path("meta_wrongIndexName.tsv"))

        with self.assertRaisesRegex(
                ValueError,
                '2 organisms have been re-defined in your metadata. '
                'Please fix!'):
            read_metadata(get_data_path("meta_ambigOrganisms.tsv"))

        obs = read_metadata(get_data_path("meta_localfiles.tsv"),
                            fp_assemblyprefix="./")
        self.assertTrue(obs['__file_exists'].all(),
                        'Not all assembly files could be located.')

        with self.assertRaisesRegex(
                ValueError,
                'Cannot find files for following assemblies'):
            read_metadata(get_data_path('meta_onemissingfile.tsv'),
                          fp_assemblyprefix="./")

    def test_get_augustus_reference_species(self):
        meta = read_metadata(
            get_data_path('meta.tsv'), skip_file_exists_test=True)
        organism = "Aspergillus sydowii"
        with self.assertRaisesRegex(
                ValueError,
                "Could not determine reference species for organism '%s'" %
                organism):
            get_augustus_reference_species(organism, meta, None)

        exp = 'kurt'
        obs = get_augustus_reference_species(organism, meta, {'augustus': {
            'default_reference_species': 'kurt'}})
        self.assertEqual(
            exp, obs, 'Project wide reference species not recovered.')

        meta = read_metadata(
            get_data_path('meta_refspec.tsv'), skip_file_exists_test=True)
        exp = 'Aspergillus nidulans'
        obs = get_augustus_reference_species(organism, meta, None)
        self.assertEqual(exp, obs, 'Metadata reference species not recovered.')

        exp = 'Aspergillus nidulans'
        obs = get_augustus_reference_species(organism, meta,  {'augustus': {
            'default_reference_species': 'kurt'}})
        self.assertEqual(
            exp, obs, "Metadata does not have precedence over config")


if __name__ == '__main__':
    main()
