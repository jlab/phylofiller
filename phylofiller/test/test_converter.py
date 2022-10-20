from unittest import TestCase, main

from skbio.util import get_data_path
from io import StringIO
import pandas as pd
from pandas.testing import assert_frame_equal

from phylofiller.converter import (
    easel_table2pd, parse_easel_output, create_CIGAR, easle2sam,
    kreport2feature)


def _str2pd(input):
    table = pd.read_csv(StringIO(input), sep=";", index_col=0, dtype=str)
    table = table.rename(columns={
        c: '' for c in table.columns if c.startswith('Unnamed: ')})
    table = table.fillna("")
    return table


class IOTests(TestCase):
    def setUp(self):
        self.fp_infernal = get_data_path('easel2sam/fw16_fssc.cm.out')
        self.fp_sam = get_data_path('easel2sam/fw16_fssc.sam')
        self.fp_nohit = get_data_path('easel2sam/nohit.cm_out')

    def tearDown(self):
        pass

    def test_parse_easel_output(self):
        obs = parse_easel_output(self.fp_infernal)
        self.assertEqual('INFERNAL', obs['software'].iloc[0])
        self.assertEqual('1.1.2', obs['software version'].iloc[0])
        self.assertEqual('Markergenes/FSSC/allFSSC.cm',
                         obs['fp_query'].iloc[0])
        self.assertEqual('Sequences_Fusarium/taxid_99000016/FW16.genome.fna',
                         obs['fp_target'].iloc[0])

        with open(get_data_path('easel2sam/exp_parse_easle_output.txt')) as f:
            exp = ''.join(f.readlines())
        assert_frame_equal(_str2pd(exp), obs)

        obs = parse_easel_output(self.fp_nohit)
        self.assertEqual(obs.shape[0], 0)

    def test_easel_table2pd(self):
        obs = easel_table2pd(
            [" rank     E-value  score  bias mdl mdl from   mdl to       seq f"
             "rom      seq to       acc trunc   gc",
             " ----   --------- ------ ----- --- -------- --------    --------"
             "--- -----------      ---- ----- ----",
             "  (1) !   4.2e-23   84.6   0.1 hmm      710      930 .]      284"
             "427      284630 + .. 0.86     - 0.44"])
        exp = (';rank;;E-value;score;bias;mdl;mdl from;mdl to;;seq from;seq to'
               ';;acc;trunc;gc\n0;(1);!;4.2e-23;84.6;0.1;hmm;710;930;.];284427'
               ';284630;+ ..;0.86;-;0.44\n')
        assert_frame_equal(_str2pd(exp), obs)

        obs = easel_table2pd(
            [" rank     E-value  score  bias  sequence    start     end   mdl "
             "trunc   gc  description",
             " ----   --------- ------ -----  --------- ------- -------   --- "
             "----- ----  -----------",
             "  (1) !  2.7e-265  885.0   0.0  Scaffold3 2966345 2965521 - hmm "
             "    - 0.54  -"])
        exp = (';rank;;E-value;score;bias;;sequence;start;end;;mdl;trunc;gc;;d'
               'escription\n0;(1);!;2.7e-265;885.0;0.0;;Scaffold3;2966345;2965'
               '521;-;hmm;-;0.54;;-\n')
        assert_frame_equal(_str2pd(exp), obs)

        obs = easel_table2pd(
            [" rank     E-value  score  bias  sequence    start     end   mdl "
             "trunc   gc  description",
             " ----   --------- ------ -----  --------- ------- -------   --- "
             "----- ----  -----------",
             "  (1) !  2.7e-265  885.0   0.0  Scaffold3 2966345 2965521 - hmm "
             "    - 0.54  -",
             "  (2) !  1.1e-243  813.4   0.0  Scaffold3 2965391 2964591 - hmm "
             "    - 0.53  -",
             "  (3) !   1.7e-34  120.9   0.0  Scaffold4 2315121 2315839 + hmm "
             "    - 0.51  -",
             "  (4) !   1.6e-26   94.6   0.0  Scaffold6 1337459 1336642 - hmm "
             "    - 0.55  -",
             "  (5) !   5.9e-08   33.1   0.0  Scaffold4 2314349 2314665 + hmm "
             "    - 0.53  -",
             "  (6) !    0.0059   16.5   0.0  Scaffold1 5920205 5920315 + hmm "
             "    - 0.60  -",
             "  (7) ?     0.037   13.9   0.1  Scaffold4  478876  478690 - hmm "
             "    - 0.55  -",
             "  (8) ?     0.071   12.9   0.0  Scaffold4 2314870 2315072 + hmm "
             "    - 0.51  -",
             "  (9) ?      0.31   10.8   0.0  Scaffold4 3659261 3659088 - hmm "
             "    - 0.56  -",
             " (10) ?         1    9.1   0.0  Scaffold6 1338065 1337981 - hmm "
             "    - 0.53  -",
             " (11) ?       1.7    8.3   0.8  Scaffold2 3674449 3674606 + hmm "
             "    - 0.59  -"])
        exp = (
            ';rank;;E-value;score;bias;;sequence;start;end;;mdl;trunc;gc;;desc'
            'ription\n0;(1);!;2.7e-265;885.0;0.0;;Scaffold3;2966345;2965521;-;'
            'hmm;-;0.54;;-\n1;(2);!;1.1e-243;813.4;0.0;;Scaffold3;2965391;2964'
            '591;-;hmm;-;0.53;;-\n2;(3);!;1.7e-34;120.9;0.0;;Scaffold4;2315121'
            ';2315839;+;hmm;-;0.51;;-\n3;(4);!;1.6e-26;94.6;0.0;;Scaffold6;133'
            '7459;1336642;-;hmm;-;0.55;;-\n4;(5);!;5.9e-08;33.1;0.0;;Scaffold4'
            ';2314349;2314665;+;hmm;-;0.53;;-\n5;(6);!;0.0059;16.5;0.0;;Scaffo'
            'ld1;5920205;5920315;+;hmm;-;0.60;;-\n6;(7);?;0.037;13.9;0.1;;Scaf'
            'fold4;478876;478690;-;hmm;-;0.55;;-\n7;(8);?;0.071;12.9;0.0;;Scaf'
            'fold4;2314870;2315072;+;hmm;-;0.51;;-\n8;(9);?;0.31;10.8;0.0;;Sca'
            'ffold4;3659261;3659088;-;hmm;-;0.56;;-\n9;(10);?;1;9.1;0.0;;Scaff'
            'old6;1338065;1337981;-;hmm;-;0.53;;-\n10;(11);?;1.7;8.3;0.8;;Scaf'
            'fold2;3674449;3674606;+;hmm;-;0.59;;-\n')
        assert_frame_equal(_str2pd(exp), obs)

    def test_create_CIGAR(self):
        with self.assertRaisesRegex(
                AssertionError,
                'Reference and read string have different length.'):
            create_CIGAR("", " ")

        with self.assertRaisesRegex(
                AssertionError,
                'Reference is not of type str.'):
            create_CIGAR([""], " ")

        with self.assertRaisesRegex(
                AssertionError,
                'Read is not of type str.'):
            create_CIGAR(" ", [""])

        exp = (
            '1=1X5=3D1=1X1=1X3=7X4=1X2=1X3=1X3=2X2=3X1=5X1=1X1=1X5=2X1=1X1=5X3'
            '=1X1=1X2=2X3=1X1=2X3=2X2=2X3I1X5=3X1=5X1=1X2=1X1=6I5=1X3=4I2=1X1='
            '5X2I1X3=1X2=2X1=3X5=1X3=')
        obs = create_CIGAR(
            ("ACCAGGAuauGCACGACAGCCUGACUGGCAUCGCCGAGGAGGCCCAGAAGACGGCUAAGGAUCU"
             "GGAGUCGAAGCUUGCCGAGGCGUUGGCCAAGGU---UGAGGAUGGAGAGAAGCAGGU------G"
             "GAGGUUCU----CGAGGCCCA--AAUCAAGGUCAAGGAUGCCGAG"),
            ("AGCAGGA...GGAu"
             "GACCAcGAgUCUGGUAUuGCCAAGGGAGCuuUGGUuCUGACCAAGGACcUuGUCAACAAGcUuG"
             "CuAAGGAGCAGGCGGAGCCaCCAGAGGACCCAUCaauGAAGAUUGGAUGGGAGGGUCUGAUCCG"
             "aGCcGGaACcAUCGAGUACCUCGAUGCuGAG"))
        self.assertEqual(exp, obs)

    def test_easle2sam(self):
        with open(self.fp_sam, 'r') as f:
            exp = ''.join(f.readlines())
            obs = easle2sam(parse_easel_output(self.fp_infernal))
            self.assertEqual(exp, obs)

    def test_parse_easel_output_format113(self):
        obs = parse_easel_output(
            get_data_path('easel2sam/SLCC2753.RF014787.cmout'))
        # just test if the new format (since 1.1.3) can be parsed to 17 hits
        self.assertEqual(obs.shape[0], 17)


class KrakenTests(TestCase):
    def setUp(self):
        self.fp_report = get_data_path('KB0060_c.k2report')

    def tearDown(self):
        pass

    def test_kreport2feature(self):
        obs = kreport2feature(self.fp_report, rank="Genus")
        self.assertEqual(obs.shape[0], 145)
        self.assertEqual(obs['#reads_clade'].sum(), 12538)

        obs = kreport2feature(self.fp_report, rank=None)
        self.assertEqual(obs.shape[0], 778)


if __name__ == '__main__':
    main()
