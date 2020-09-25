from unittest import TestCase, main

from skbio.util import get_data_path
from io import StringIO
from pandas.testing import assert_frame_equal

from phylofiller.converter import *


def _str2pd(input):
    table = pd.read_csv(StringIO(input), sep=";", index_col=0, dtype=str)
    table = table.rename(columns={c: '' for c in table.columns if c.startswith('Unnamed: ')})
    table = table.fillna("")
    return table

class IOTests(TestCase):
    def setUp(self):
        self.fp_infernal = get_data_path('easel2sam/fw16_fssc.cm.out')
        #self.exp_md5 = '1b89479f6bbae9a684b75b71e0cb492e'


    def tearDown(self):
        pass

    def test_easel2sam(self):
        obs = easel2sam(self.fp_infernal)
        self.assertEqual('INFERNAL', obs['software'].iloc[0])
        self.assertEqual('1.1.2', obs['software version'].iloc[0])
        self.assertEqual('Markergenes/FSSC/allFSSC.cm', obs['fp_query'].iloc[0])
        self.assertEqual('Sequences_Fusarium/taxid_99000016/FW16.genome.fna', obs['fp_target'].iloc[0])

        exp = ';rank;;E-value;score;bias;mdl;mdl from;mdl to;;seq from;seq to;;acc;trunc;gc;target name;query sequence;target sequence;model name;model clen;software;software version;fp_query;fp_target\n0;(1);!;4.2e-23;84.6;0.1;hmm;710;930;.];284427;284630;+ ..;0.86;-;0.44;Scaffold6;AGAGACCGAUAGCGCACAAGUAGAGUGAUCGAAAGAUGAAAAGAACUUUGAAAAGAGAGUUAAAaAGUACGUGAAAUUGUUGAAAGGGAAGCGCUUGUGACCAGACUUGGGCUuGGUUGAUCAUCCgGGGUUCUCCCuGGUGCACUCUUCCGGCccAGGCCAGCAUCAGUUcGCCcuGGGGGAuAAAGGCuuCGGGAAuGUGGCUCuCUCCGGGGAGUGUU;AAAUGCCGAUGGUGUACGAUUA--CUGAUCAAAAGAUGAA-AGAACUUAAGAAAGAGAGCUAA-AAGUACGUGAG---GUUGAAAG-GAAACGCUUGGCACUAGAUUUGGGCUUGGUUGAUU-----GAGGUCCCUUCAGUGCACUCU---GGCCUAUGCUGGUAUUAACUUGCCUAAGGGAAUAAAGGCUUG-GGGAUGUAGUUCUCUCCUGGGAGUGUU;FSSC_ITS+28S-rDNA--956-bp.fna;930;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(2);!;7e-17;64.0;0.0;hmm;802;928;..;4687313;4687190;- ..;0.91;-;0.54;Scaffold3;GCUUGUGACCAGACUUGGGCUuGGUUGAUCAUCCgGGGUUCUCCCuGGUGCACUCUUCCGGCccAGGCCAGCAUCAGUUcGCCcuGGGGGAuAAAGGCuuCGGGAAuGUGGCUCuCUCCGGGGAGUG;GCUCGAGACCAGACCUGGGCUUGGUUGAUCAUCCAG-GUUCUUCCUGGUGCACUCUUCUGGCUUAUGCCAGCAUCGAUUUGUCUUGAGGAACAAUGGUUUCCGGGAUGUGACUA--CCCGGGCAGUG;FSSC_ITS+28S-rDNA--956-bp.fna;930;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(3);?;0.016;16.6;0.0;hmm;752;808;..;1013790;1013850;+ ..;0.78;-;0.39;Scaffold1;GAACUUUGAAAAGAGAGUUAAAaAGUACGUGAAAUUGUUGA....AAGGGAAGCGCUUGUG;AAAUUUUAAAAAGAGAGUUAAGCAGUACGUGAAAUUGUCGUucgcAUUGGACCCGUUGGUG;FSSC_ITS+28S-rDNA--956-bp.fna;930;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(4);?;1.1;10.6;1.0;hmm;693;821;..;2222270;2222153;- ..;0.73;-;0.46;Scaffold5;AAAGCUAAAUACCGGCCAGAGACCGAUAGCGCACAAGUAGAGUGAUCGAAAGAUGAAAAGAACUUUGAAAAGAGAGUUAAAaAGUACGUGAAAUUGUUGAAAGGGAAGCGCUUGUGAC...CAGACUUGGGC;AAGGCAAGCUACCACCUUCAAAGAAAUAGUGUCCGGCAAGAAGGAUCGAAAGG--------------AGAAGGAGGUUAUAGAGGACGGGAAAUUGGUGUAAUAAAAGCGUCAAUGACagaCAGACGAGGUC;FSSC_ITS+28S-rDNA--956-bp.fna;930;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(5);?;5.7;8.2;0.7;hmm;756;861;..;3020655;3020767;+ ..;0.79;-;0.58;Scaffold2;UUUGAAAAGAGAGUUAAAaAGUACGUGAAAUUGUUGAAAGGGAAGCGCUUGUGACCAG.......ACUUGGGCUuGGUUGAUCAUCCgGGGUUCUCCCuGGUGCACUCUUCCG;UCUGGAGAGGGAAUAAAAAAGCUGGAGAAAUGGCGCGGGGGAAAGGGCUCUUUGCCAGgggggagUGUGGAGCUUGGGCGAUGAUGACGGGUGAUACCUGGGGCAUGAUGCGG;FSSC_ITS+28S-rDNA--956-bp.fna;930;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(1);!;2.7e-265;885.0;0.0;hmm;2;826;..;2966345;2965521;- ..;0.99;-;0.54;Scaffold3;ACcCCuAUuGGACGAGAUGGAAAGCUCGCCAAGCCCCGUCAGCUaCACAACACCCAUUGGGGuCUGGUcUGuCCAGCCGAGACgCCcGAGGGUCAGGCUUGUGGuCUGGUCAAGAACUUGUCcCUGAUGUGcUAcGUCAGUGUCGGcUCUCCcUCcGAACCucUGAUuGAGUUCAUGAUCAACCGAGGUAUGGAaGUCGUGGAAGAGUACGAgCCCCUGAGAUAcCCgCAUGCuACCAAGAUcUUuGUCAAuGGUGUcUGGUGuGGUGUcCACUCgGACCCcAAGCAUCUCGUCAGcCAGGUccUGGACACaCGACGAAAGUCGUAccUGCAGUAcGAGGUGUCgCUuGUUCGUGACAUUCGAGAuCGAGAGUUCAAGGUCUUCUCCGACGCuGGCCGAGUCAUGAGgCCGGUCUUUACGGUuCAGCAGGAGGAuGACCAcGAgUCUGGUAUuGCCAAGGGAGCuuUGGUuCUGACCAAGGACcUuGUCAACAAGcUuGCuAAGGAGCAGGCGGAGCCaCCAGAGGACCCAUCaauGAAGAUUGGAUGGGAGGGUCUGAUCCGaGCcGGaACcAUCGAGUACCUCGAUGCuGAGGAAGAGGAGaCGGCuAUGAUUUGCAUGACuCCUGAGGAuCUuGAcCUCUAuCGcAUGCAAAAGGCuGGUUACGUcGUaGAuGAcGAUAACACGGACGACCCcAACAGGAGaUUGAAGACcAAGACgAACCCCACAACUCACAUGUACACUCAUUGUGAGAUUCACCCcAGuAUGAUUCUuGGCAUuUGuGCCAGUAUCAUUCCcUUCCCcGAUCACAACCAGGUAuACgAC;ACCCCUAUUGGACGAGAUGGAAAGCUCGCCAAGCCCCGUCAGCUACACAACACCCAUUGGGGUCUGGUCUGUCCAGCCGAGACACCCGAGGGUCAGGCUUGUGGUCUGGUCAAGAACUUGUCCCUGAUGUGUUACGUCAGUGUCGGUUCUCCUUCUGAACCUCUGAUCGAGUUCAUGAUCAACCGAGGCAUGGAAGUCGUGGAAGAGUACGAGCCCCUGAGAUACCCGCAUGCUACCAAGAUCUUUGUCAAUGGUGUCUGGUGCGGUGUCCACUCGGACCCCAAGCAUCUCGUCAGCCAGGUCCUGGACACACGACGAAAGUCGUACCUGCAAUACGAGGUGUCGCUCGUCCGUGACAUUCGAGAUCGAGAGUUCAAGGUCUUCUCCGACGCUGGCCGAGUCAUGAGGCCAGUCUUUACGGUUCAGCAGGAGGAUGACCACGAGUCUGGUAUCGCCAAGGGAGCUUUGGUUCUGACCAAGGACCUUGUCAACAAGCUUGCUAAGGAGCAGGCGGAGCCACCAGAGGACCCAUCAAUGAAGAUUGGAUGGGAGGGUCUGAUCCGAGCCGGAACCAUCGAGUACCUCGAUGCUGAGGAAGAGGAGACGGCUAUGAUUUGCAUGACUCCUGAGGAUCUUGACCUCUAUCGCAUGCAAAAGGCUGGUUACGUCGUAGAUGACGAUAACACGGACGACCCCAACAGGAGAUUGAAGACCAAGACGAACCCCACAACUCACAUGUACACUCAUUGUGAGAUUCACCCCAGUAUGAUUCUUGGCAUUUGUGCUAGUAUCAUCCCCUUCCCCGAUCACAACCAGGUAUGCGCC;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(2);!;1.1e-243;813.4;0.0;hmm;787;1588;.];2965391;2964591;- ..;0.96;-;0.53;Scaffold3;CAGUAUCAU.UCCcUUCCCcGAUCACAACCAGGUAuACgACaCGaUCCAUGGAGUUCCUCAAGUUCCGUGAaCUGCCaGCcGGuCAGAACGCCAUcGUCGCuAUCGCUUGcUACUCuGGUUAcAACCAGGAAGAUUCCGUCAUuAUGAACCAGAGuAGUAUCGAuCGAGGCuUGUUCCGcAGUcUGUUCUUCaGAUCcUACUCuGAcCAGGAGAAGAAGGUCGGuCUgAACUACACGGAAGUGUUuGAGAAGCCCUUCCAGCAGUCGACgCUUCGuAUGAAGCAcGGUACcUACGACAAGCUGGAcGAGGAuGGUAUcGUGGCcCCcGGuGUGCGAGUGUCgGGUGAaGAUAUCAUcAUCGGCAAGACuGCGCCGAUuGAuCAAGAGAACCAGGAUCUgGGuACCAGGACaACgGUGCACCAGCGUCGUGAUAUCUCcACgCCGCUGCGAAGuACcGAgAACGGuAUCGUuGAuucGGUCAUUGUGACuGUCAAuGCcGACAACGUCAAGUAcGUCAAGGUCCGUGUGAGGACGACCAAGAUUCCuCAGAUUGGuGACAAGUUuGCCUCUCGUCACGGACAGAAGGGUACCAUUGGUGUuACcUACcGaCAGGAGGAuAUGCCcUUcagCAGGGAGGGuGUGACaCCaGACAUuAUCAUuAACCCCCACGCCAUUCCgUCGCGAAUGACAAUUGCcCAUUUGAUUGAAUGCCUCCUCAGUAAGGUGUCAACgCUcGAAGGCAUGGAGGGUGAuGCaACaCCuUUCACcGAuGUCACuGUcGACUCcGUuUCGGAGCUGCUCCG;CAACAUUCUuUACUAUCCGCAAA--AACCUCUGGCGACAACACGAUCCAUGGAGUUCCUCAAGUUCCGUGAACUACCAGCCGGUCAGAAUGCCAUCGUCGCUAUCGCUUGCUACUCUGGUUACAACCAGGAAGAUUCCGUCAUUAUGAACCAGAGUAGUAUCGAUCGAGGCCUGUUCCGCAGUCUGUUCUUCAGAUCCUACUCUGACCAGGAGAAGAAGGUCGGACUGAACUACACGGAAGUGUUCGAGAAGCCCUUCCAGCAGUCGACGCUUCGUAUGAAGCACGGUACUUACGACAAGCUGGACGAGGAUGGUAUCGUGGCCCCCGGUGUCCGAGUGUCAGGUGAAGACAUCAUCAUCGGCAAGACUGCGCCGAUUGACCAAGAGAACCAGGAUCUGGGUACCAGGACGACGGUGCACCAGCGUCGUGAUAUCUCCACGCCGCUGCGAAGUACCGAGAACGGUAUCGUCGAUUCGGUCAUUGUGACUGUCAAUGCCGACAACGUCAAGUACGUCAAGGUCCGUGUGAGGACGACCAAGAUUCCUCAGAUUGGUGACAAGUUCGCCUCUCGUCACGGACAGAAGGGUACAAUUGGUGUCACCUACCGGCAGGAGGAUAUGCCCUUCAGCAGGGAGGGUGUGACACCAGACAUUAUCAUUAACCCCCACGCCAUUCCCUCGCGAAUGACAAUUGCCCAUUUGAUUGAAUGCCUCCUCAGUAAGGUGUCAACGCUCGAAGGCAUGGAGGGUGAUGCAACGCCUUUCACCGAUGUCACCGUCGACUCCGUUUCGGAGCUGCUCCG;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(3);!;1.7e-34;120.9;0.0;hmm;828;1525;..;2315121;2315839;+ ..;0.66;-;0.51;Scaffold4;CGaUCCAUGGAGUUCCUCAAGUUCC..................GUGAaCUGCCaGCcGGuCAGAACGCCAUcGUCGCuAUCGCUUGcUACUCuGGUUAcAACCAGGAAGAUUCCGUCAUuAUGAACCAGAGuAGUAUCGAuCGAGGCuUGUUCCGcAGUcUGUUCUUCaGAUCcUACUCuGAcCAGGAGAAGAAGGUCGGuCUgAACUACACGGAA.GUGUUuGAGAAGCCCUUCCAGCAGUCGACgCUUCGuAUGAAG....CAcGGU..........ACc......UACGACAAGCUGGAcGAGGAuGGUAUcGUGGCcCCcGGuGUGCGAGUGUCgGGUGAaGAUAUCAUcAU..............CGGCAAGACuGCGCCGAUuGAuCAAGAGAACCAGGAUCUgGGuACCAGGAC.aACgGU.GCACCAGCG....UCGUGAUAUCUCcACgCCGCUGCGAAGuACcGAgAACGGuAUCGUuGAuucGGUCAUUGUGACuGUCAAuGCcGACAACGUCA...AGUAcGUCAAGGUCCGUGUGAGGACGACCAAGAUUCCuCAGAUUGGuGACAAGUUuGCCUCUCGUCACGGACAGAAGGGUACCAUUGGUGUuACcUACcGaCAGGAGGAuAUGCCcUUcagCAGGGAGGGuGUGACaCCaGACAUuAUCAUuAACCCCCACGCCAUUCCgUCGCGAAUGACAAUUGCcCAUUUGAUUGAAUGCCUCCUCAGUAAGGUGUCAACgCUcGAAGG;CGACCCAUGGUUAUCUCCAAGACCAuccagcuuauuggcuaugAUAAGCUCCCCGCAGGCCAGAACGCGACUGUCGUCGUCAUGUCGUACUCUGGAUAUGAUAUCGAAGAUGCUUUGGUUCUGAACAAGGCGUCGAUCGACAGAGGAUUUGGACGCUGCCAGGUCUUCCGCAAAUACACG-ACCGAGCUACA-----------AAAAUACCCUAAUgGUCGCCGAGA-GCGCAUC-GG----CGAUCCCCAAAAUGAAGaaggCAAGGUcaagcaaagaAUCaagaagCAUGAAGGUCUCGACGACGAUGGUUUGG--------------CCAUUGUGGGAUACAGAAUUCAUAAUggcgaggccaugauCAAGAAGGAAACACCCCUUGACCA-GACAACCACC-----GGCAUCGGAAUgGAUCGUgGACCCAGCGaauaCCGUGAUUCUUCAGUCUCGUACCGUAUUGCUGAUCCGGCAUACAUUGACAAGGUGAUGGUAUC---CCAAACAGAGAAGGACAcuaCAGUUAUCAAGGUCCAGACCAGACAGACUCGUCGUCCAGAGCUUGGAGACAAGUUUUCAUCUCGUCACGGUCAAAAGGGUGUGGUUGGCAUCAUUGUUGAUCAGGAGGAUCUGCCCUUCUCUGACAAGGGCCUGACCCCUGAUAUCAUCAUGAACCCUCACGGUUUCCCGUCUCGAAUGACGGUCGGCAAGCUGCUCGAGUGCUUGACGGGCAAGGCUUCCAUUAUCCAUGG;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(4);!;1.6e-26;94.6;0.0;hmm;761;1564;..;1337459;1336642;- ..;0.63;-;0.55;Scaffold6;CCcAGuAUGAUUCUuGGCAUuUGuGCCAGUAUCAUUCCcUUCCCcGAUCACAACCAGGUAuACgACaCGaUCCAUGGAGUUCCUC...A..A........GUUCCGUGAaCUGC.CaGCcGGuCAGAACGCCAUcGUCGCuAUCGCUUGcUACUCuGGUUAcAACCAGGAAGAUUCCGUCAUuAUGAACCAGAGuAGUAUCGAuCGAGGCuUGUUCCGcAGUcUGUUCUUCaGAUCcUACUCuG.AcCAGGAGAAGAAGGUCGGuCUgAACUA...CACG.....G.........A..AGUGUUuGAGAAGCCCUUC...CAGCA.GUCGACgCUUCGuAUGAAGCAcGGUAC.cUACGACAAGCUGGAcGAGGAuGGUAUcGUGGCcCCcGGuG..UGCGAGUGUCgGGU.GAaGAUAUCAUcAUCG...GCAAGACuGCGCCGAUuGAuCAAGAGAACCAGGAUCUgG.GuACCAGGACaACgGUGCACCAGCGUCGUGAUAUCUCcACgCCGCUGCG..AAGuAC........cGAgAACGGuAUCGUuGAuucGGUCAUUGUGACuG...UCAAuGCcGACAACG............UCAAGUAcGUCAAGGUCCGUGUGAGGACGACCAAGAUUCCuCAG.AUUGGuGACAAGUUuGCCUCUCGUCACGGACAGAAGGGUACCAUUGGUGUuACcUACcGaCAGGAGGAuAUGCCcUUcagCAGGGAGGGuGUGACaCCaGACAUuAUCAUuAACCCCCACGCCAUUCCgUCGCGAAUGACAAUUGCcCAUUUGAUUGAAUGCCUCCUCAGUAAGGUGUCAACgCUcGAAGGCAUGGAGGGUGAuGCaACaCCuUUCACcGAuGUCACuGU;CCCAGUACUAAUCUUCGCUACCGAACCGACAACAAGUCAUACC--GCAUACAAACUGGACAAACACCCGU-CGUGCGAGCACCUCuucAcaAcacauaugGUUUCGACAACUUCcCCAACGGCAUGAACGCAGUUGUUGCCGUCAUCUCGUAUACUGGAUAUGACAUGGACGACGCCAUGAUUCUCAACAAGAGCGCCCAUGAGCGUGGUUUUGGGCA----------------UGGUACUAUCuACAAGACAAAGAAGAUCUCUCUCAAGGAugaCUCGcgaacAaaggccaccAagAGUGUUACCAAGGCUUUUGgcuUCGCAcCUCACAGCUACGUG---AGCGCGUCAUaCCAGGGAAUGCUGGAUGACGACGGCCU---GCCUCACGUGGgcCGCAUGAUCCAGGAgGGAGAUGUAAUCUGCGcguGGCACACAGUGACGCCCGACUACAAUGGCAAG---CUGGuGAACCUGGAUGGC-------------------AU---CACU-CACUACGaaAAGUACaaggauggCGAGACGGGCUUUGUCGAAGAGGUCCGCCUGAUCGgagCCGAAACCGGCAACGagccucuccagaCCAUCUCUGUCAAGUUCCGUGUU-------CCUCGAUCCCCCAUcAUCGGUGACAAGUUCUCGUCCAGACACGGACAGAAGGGUGUCGCCUCGCAAAAGUGGCCAACGGUGGAUCUGCCCUUCUCGGAGACUGGCAUCCAACCGGACAUCAUCAUUAACCCCCACGCUUUCCCUUCCCGUAUGACGAUAGGCAUGUUUGUGGAAUCGCUCGCCGGCAAGGCUGGCGCGCUGCAUGGCCUUGCGCAAGACUCGACACCGUUCAAGUUUGACGAGGA;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(5);!;5.9e-08;33.1;0.0;hmm;17;320;..;2314349;2314665;+ ..;0.69;-;0.53;Scaffold4;GAUG.GAAAGCUCGCCAAGCCCCGUCAGCUaCACAACACCCAUUGGGGuCUGGUcUGuCCAGCCGAGACgCCcGAGGGUCAGGCUUGUGGuCUGGUCAAGAACUUGUCcCUGAUGUGcUAcGUCA......GUGUCGGcUCUCCcUCcGAACCucUGAUuGAGUUCAUGAUCAACCGAGGUAUGGAaGUCGUGG......AAGAG...UACGAgCCCCUGAGAUAcCCgCAUGCuACCAAG...AUcUUuGUCAAuGGUGUcUGGUGuGGUGUcCACUCgGACCCcAAGCAUCUCGUCAGcCAGGUccUGGACAC.aCGACGAA;GACGcGAAAGGUCUCAGGCCCUCGUGCGCUCCAACCGUCGCAAUGGGGUAUGCUGUGCACUUCGGAUACACCUGAAGGAGAAGCCUGCGGUCUGGUGAAGAACUUGGCUUUGAUGACACAUAUCAccacaaAUGUCGAG------GAAGGACCCGUGAAGGAGACUAUCCUGACGAUCGACAAGGAAGUUGAGGcuaucgAGAAAuucUCUGGCUCAAUGAUGCACCGGGAAGGGAGCUACgucAUCCAUGUCAAUGGUACGCCAUUC-GCAUUGACACGGAAUCCGAAGAGGUUCGCACAAAGGUUCAGGACAUuGCGACGAA;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(6);!;0.0059;16.5;0.0;hmm;438;542;..;5920205;5920315;+ ..;0.80;-;0.60;Scaffold1;ACCAcGAgUCUGGUAUuGCCAAGGGAGCuuUGGUuCUGACCAAGGACcU......uGUCAACAAGcUuGCuAAGGAGCAGGCGGAGCCaCCAGAGGACCCAUCaauGAAGA;CCCUGAAGGCUGGCGUUGACAAGGGCAUCUUUGAGCAGCCCAAGGGUCCuuccggUGGCACCAAGCUUGCCAAGAAGCAGCCUGAGCCCAAGAAGGCUGCCGCUGAGAAGA;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(7);?;0.037;13.9;0.1;hmm;551;715;..;478876;478690;- ..;0.62;-;0.55;Scaffold4;GAGGGUCUGAUCCGaGCcGGaACcAUCGAGUACCUCGAUGCuGAGGAAGAGGAG.aCGGCuAUGAUUUGCAUGACuCCUGAGG...........AuCUu.......GA..cCUCUAuCGcAUGCAAAA.....GGCuGGUUACGUcGUaGAuGAcGAUAACACGGACGACCCcAACAG....GAGaUUGAAGACc;GAGGUUGAGAUUCCAGAGGAGAUUGAAGAGGAAGAAGAGGUCGAAGAGGAAGAGgAUGACGAGGAUA-------CUCUUGAGGucgagggcgagAUCAUcgcugaaGAgcCAUCCGGCGAA-GAAGAAgagccGGCCGUAGAAGUUGUGGAGGAAGAUGACCCGGACGACACCAUCAGuccgGAAGCUGAAGAGG;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(8);?;0.071;12.9;0.0;hmm;669;899;..;2314870;2315072;+ ..;0.59;-;0.51;Scaffold4;UcGUaGAuGAcGAUAACACGGACGACCCcAACAGGAGaU.UGAAGACcAAGACgAACCCCACAACUCACAUGUACACUCAUUGUGAGAUUCACCCcAGuAUGAUUCUuGGCAUuUGuGCCAGUAUCAUUCCcUUCCCcGAUCACAACCAGGUAuACgACaCGaUCCAUGGAGUUCCUCAAGUUCCGUGAaCUGCCaGCcGGuCAGAACGCCAUcGUCGCuAUCGCUUGcUAC;UGGAUGUGAACGAGGAGAAUGACGCCCUCAUCACCAUCUaCGAAGACCAAG----------UGACGCAAAGUACCACGCAUCUGGAGAUUGAGCCGUUCACCAUCCUGGGUGCAGUUGCAGGACUGAUUCCAUUCCCUCACCACAACCAAU----CGCCU-------CGUAAUACCUAC----CAAUGUGC----UAUGGGUAAACAAGCCAUUGGUGCCAUCGCCUACAAC;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(9);?;0.31;10.8;0.0;hmm;421;600;..;3659261;3659088;- ..;0.62;-;0.56;Scaffold4;GGUuCAGCAGGAGGAuGACCAcG...AgUCUGGUAUuGCCAAGGGAGCuuUGGUuCUGACCAAGGACcUuGUCAACAAGcUuGCuAAGGAGCAGGCGGAGCCaCCAGAGGACCCAUCaauGAAGAUUGGAUGGGAGGGUCUGAUCCGaGCcGGaACcAUCGAGUACCUCGAUGCuGAGGAAGA;GGACCUGGAUGACGAUGCCCGUGcccAGCGCGACGAGGAGCAGGGAGUCAAGGACAUUAUCAAUGACCUCAAAAAGCAGCUCGCCAACAAGCAGGCUGAGCUGGCC-CGCACCCAGAAC--AAGAUCUU----GCGUACUCGUAUCGAGCAG--ACCAUCAAGCAGCUCAAGUCGGAGAUUGA;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(10);?;1;9.1;0.0;hmm;56;140;..;1338065;1337981;- ..;0.74;-;0.53;Scaffold6;CAUUGGGGuCUGGUcUGuCCAGCCGAGACgCCcGAGGGUCAGGCUUGUGGuCUGGUCAAGAACUUGUCcCUGAUGUGcUAcGUCA;UCAUGGGGCUUCAUGUGUCCUGUCCACACGCCUGAUGGUGCGCCUUGUGGUCUUCUGAACCAUCUUGCGCACAAGUGCAAGAUUA;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(11);?;1.7;8.3;0.8;hmm;426;595;..;3674449;3674606;+ ..;0.64;-;0.59;Scaffold2;AGCAGGA...GGAuGACCAcGAgUCUGGUAUuGCCAAGGGAGCuuUGGUuCUGACCAAGGACcUuGUCAACAAGcUuGCuAAGGAGCAGGCGGAGCCaCCAGAGGACCCAUCaauGAAGAUUGGAUGGGAGGGUCUGAUCCGaGCcGGaACcAUCGAGUACCUCGAUGCuGAG;ACCAGGAuauGCACGACAGCCUGACUGGCAUCGCCGAGGAGGCCCAGAAGACGGCUAAGGAUCUGGAGUCGAAGCUUGCCGAGGCGUUGGCCAAGGU---UGAGGAUGGAGAGAAGCAGGU------GGAGGUUCU----CGAGGCCCA--AAUCAAGGUCAAGGAUGCCGAG;FSSC_RPB2--1588-bp.fna;1588;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(1);!;2.6e-192;642.5;0.1;hmm;1;632;[];6427210;6427837;+ ..;0.99;-;0.54;Scaffold2;AGUCGACCACCGUAAGUCaaaCCCUCaUcGCGAuCUGCUuAUCUCGGGUCGUGGAACCCCGCCuGguauCUCGGGCGGGGUauuCAUCAgUCACUucaUGCUGACAAUCAucuACAGACCGGUCACUUGAUCUACCAGUGCGGUGGUAUCGACAAGCGAACCAUCGAGAAGUUCGAGAAGGUUGGUGACAUCUCCCCCGAUCGCGCCUUGCUaUuccaCAUCGAAUUCCcCGUCGAAUUCCCUCCcuCGCGAuaCGCuCUGCGCCCGCUuCUCCCGAGUCcCAAAAuUUUUGCGGUuCGACCGuaAuUUUUUUuGGUGGGGCAUuUACCCCGCCaCUCGgGcGACGUUGGACAAAGCCCUGAUCCCUGCACAC.AAAAACACCAAACCCUCUUGGCGCGCAUCaUCACGUGGUucaCaaCAGAcgCUaACcGguuCAACAAuAGGAAGCCGCUGAGCUCGGuAAGGGUUCCUUCAAGUACGCCUGGGUCCUUGACAAGCUCAAGGCCGAGCGUGAGCGUGGUAUCACCAUCGAcAUUGCcCUCUGGAAGUUCGAGACUCCCCGCUACUAUGUCACCGUCAUUGGUAUGUcGCuGucgucuCucUCacuCauGUCUCacCaCUAACaAUCaACA;AGUCGACCACCGUAAGUCAAACCCCCAUCGCGAUCUGCUUAUCUCGGGUCGUGGAACCCCGCCUGGUAUCUCGGGCGGGGUACUCAUCAGUCACUUCAUGCUGACAAUCAUCUACAGACCGGUCACUUGAUCUACCAGUGCGGUGGUAUCGACAAGCGAACCAUCGAGAAGUUCGAGAAGGUUGGUGACAUCUCCCCCGAUCGCGCCUUGCUAUUCCACAUCGAAUUCCCCGUCGAAUUCCCUCUUCCGCGACACGCUCUGCGCCCGCUUCUCCCGAGUCCCAAAAAUUUUGCGGUUCGACCGUAAUUUUUUU-GGUGGGGCAUUUACCCCGCCACUCGGGCGACGUUGGACG-AGCCCUGAACCCUGCACACaAAAAACACCAAACCCUCUUGGCGCG---CAUCACGUGGUUCACAACAGACACUAACAGGUUCAACAAUAGGAAGCCGCUGAGCUCGGUAAGGGUUCCUUCAAGUACGCCUGGGUCCUUGACAAGCUCAAGGCCGAGCGUGAGCGUGGUAUCACCAUCGAUAUUGCUCUCUGGAAGUUCGAGACUCCCCGCUACUAUGUCACCGUCAUUGGUAUGUCGCCCUCAUCUCUCUCAAUCACGUCUCAUCAUUAACAAUCAACA;FSSC_TEF1-665-bp.fna;632;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(2);?;0.042;16.1;0.0;hmm;468;568;..;4896877;4896777;- ..;0.66;-;0.53;Scaffold3;UCCUUCAAGUACGCCUGGGUCCUUGACAAGCUCAAGGCCGAGCGUGAGCGUGGUAUCACCAUCGAcAUUGCcCUCUGGAAGUUC......GAGACUCCCCGCUACUA;UUCCUCGAGUAUG---GUGCCAUUGACAAGGCUCCCGAGGAGCGAAAGCGUGGUAUCACCAUUUCCACUGCUCACAUCGAGUACucuaccGACAA---CCGCCACUA;FSSC_TEF1-665-bp.fna;632;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(3);?;0.048;15.9;0.0;hmm;438;541;..;3187913;3188016;+ ..;0.85;-;0.54;Scaffold2;CAAuAGGAAGCCGCUGAGCUCGGuAAGGGUUCCUUCAAGUACGCCUGGGUCCUUGACAAGCUCAAGGCCGAGCGUGAGCGUGGUAUCACCAUCGAcAUUGCcCU;CGAAGGCAAGCUGAAAAGUCUGGCAAGCAGUCUUUUGCUCUUGCAUGGGUCAUGGACCAGAGAAGUGAAGAGCGAGAACGCGGCGUGACCAUUGACAUCGCCAC;FSSC_TEF1-665-bp.fna;632;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(4);?;0.17;14.1;0.0;hmm;78;309;..;2145448;2145215;- ..;0.65;-;0.54;Scaffold7;GGGUauuCAUC....AgUCACUucaUGCUGACAAUCAucuACAGACCGGUCACUUGAUCUACCAGUGCGGUGGUAUCGACAAGCGAACCAUCGAGAAGUUCGAGAAGG........UUGGUGACAUCUCCCCCGAUCGCGCCUU....GCUaUuccaCAUCGAAUUCCcCGUCGAAUUCCCUCCcuCGCGAuaCGCuCUGCGCCCGCUuCUCCCGAGUCcCAAAAuUUUUGCGGUuCGACCGuaAuUU;GGGGAUCCGUCggggAGUCCCAUCGUGCGUAGAGUCCGUUAGAGAUUGAU-----GAGCUGACAGUG-------AAAGACUUGGGAA-CAUGGAGAGGGAAGAGAGACuccuacucUCCGAAUCAUCCGCCGCAACUGCGUUAUcuccGCUGGUCAGCA-CCAAUCAUGGGUCGCUUGACCCUCCGCACCACCCGCUCAUCGUCCGCUCCUCUCAAGUCUCAAUUCAUACGGGCUACGUAUAUAAAUU;FSSC_TEF1-665-bp.fna;632;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(5);?;0.56;12.4;0.0;hmm;298;505;..;573210;573407;+ ..;0.69;-;0.51;Scaffold7;CGACCGuaAuUUUUUUuGGUGGGGCAUuUACCCCGCCaCUCGgGcGACGUUGGACAAAGCCCUGAUCCCUGCACACAAAAACACCAAACCC......UCUUGGCGCGCAUCaUCACGUGGUucaCaaCAGAcgCUaACcGguuCAAC...AAuAGGAAGCCGCUGAGCUCGGuAAGGGUUCCUUCAAGUACGCCUGGGUCCUUGACAAGCUCAAGGC;CCACCGUCAUAUU----GGUC-AGCUUUUACCGAGCCGUUAUAGAG-UGUCGUA------CCCGAUCG--GGUGACAAUUCCAUCAAUGCCgaggauUAUUGAC--GAAUCGGAACGUGAGUGUCGGAACAAGA--UCCAGUUCCACgacCAAUGAGAG-CCCUGAGCUCGGCGAGGGCCACUUAACGUUAGGCUUGCACUUUGUCUAUUUCAAGGU;FSSC_TEF1-665-bp.fna;632;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(6);?;4.9;9.3;0.0;hmm;150;293;..;6008452;6008603;+ ..;0.71;-;0.57;Scaffold1;CGACAAGCGAACCAUCGAGAAGUUCGAGAAGGUUGGUGACAUCUCCCCCGAUCGCGCCUUGCUaUuccaCAUCGAAUUCCcCGUCGAAUUC.....CCUCCcuCGCGAuaCGCuCU...GCGCCCGCUuCUCCCGAGUCcCAAAAuUUUUGC;CUUCAAGCCUGCCUAUGAGAAGUUCGCCAGCCUGGGUGACAAGAUCCCCAAGCUCGUCUUCACCACCUACUUCGGUGACAUCGUCCACAACcucgaCCUCCUCCCCAAGGAGGUCUaugGUGUCCACAUCGACCUUGUCCGCAACCCUGAGC;FSSC_TEF1-665-bp.fna;632;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n0;(7);?;9.7;8.3;0.1;hmm;55;249;..;6008632;6008801;+ ..;0.66;-;0.62;Scaffold1;AACCCCGCCuGguauCUCGGGCGGGGUauuCAUCAgUCACUucaUGCUGACAAUCAucuACAGACCGGUCACUUGAUCUACCAGUGCGGUGGUAUCGACAAGCGAACCAUCGAGAAGUUCGAGAAGGUUGGUGACAUCUCCCCCGAUCGCGCCUUGCUaUuccaCAUCGAAUUCCcCGUCGAAUUCCCUCCcuCG;CCCCAAGACUGUCCUCUCUGCCGGUGUCGUCGACGGCCGCAACAUCUGGAAGACCAACCUGAAGCGCGC------------CAUCG-----AGAUUGUCGAGACUGCCAUCCAGAAGCUCGGCAAGGACCGUGUCAUUGCCGCCACCUCCUCUUCCCUCCUCCACA--------CCCCUCACACUCUCGCCAGCG;FSSC_TEF1-665-bp.fna;632;INFERNAL;1.1.2;Markergenes/FSSC/allFSSC.cm;Sequences_Fusarium/taxid_99000016/FW16.genome.fna\n'
        assert_frame_equal(_str2pd(exp), obs)


    def test_easel_table2pd(self):
        obs = easel_table2pd(
            [" rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       acc trunc   gc",
             " ----   --------- ------ ----- --- -------- --------    ----------- -----------      ---- ----- ----",
             "  (1) !   4.2e-23   84.6   0.1 hmm      710      930 .]      284427      284630 + .. 0.86     - 0.44",])
        exp = ';rank;;E-value;score;bias;mdl;mdl from;mdl to;;seq from;seq to;;acc;trunc;gc\n0;(1);!;4.2e-23;84.6;0.1;hmm;710;930;.];284427;284630;+ ..;0.86;-;0.44\n'
        assert_frame_equal(_str2pd(exp), obs)

        obs = easel_table2pd(
            [" rank     E-value  score  bias  sequence    start     end   mdl trunc   gc  description",
             " ----   --------- ------ -----  --------- ------- -------   --- ----- ----  -----------",
             "  (1) !  2.7e-265  885.0   0.0  Scaffold3 2966345 2965521 - hmm     - 0.54  -"])
        exp = ';rank;;E-value;score;bias;;sequence;start;end;;mdl;trunc;gc;;description\n0;(1);!;2.7e-265;885.0;0.0;;Scaffold3;2966345;2965521;-;hmm;-;0.54;;-\n'
        assert_frame_equal(_str2pd(exp), obs)

        obs = easel_table2pd(
            [" rank     E-value  score  bias  sequence    start     end   mdl trunc   gc  description",
             " ----   --------- ------ -----  --------- ------- -------   --- ----- ----  -----------",
             "  (1) !  2.7e-265  885.0   0.0  Scaffold3 2966345 2965521 - hmm     - 0.54  -",
             "  (2) !  1.1e-243  813.4   0.0  Scaffold3 2965391 2964591 - hmm     - 0.53  -",
             "  (3) !   1.7e-34  120.9   0.0  Scaffold4 2315121 2315839 + hmm     - 0.51  -",
             "  (4) !   1.6e-26   94.6   0.0  Scaffold6 1337459 1336642 - hmm     - 0.55  -",
             "  (5) !   5.9e-08   33.1   0.0  Scaffold4 2314349 2314665 + hmm     - 0.53  -",
             "  (6) !    0.0059   16.5   0.0  Scaffold1 5920205 5920315 + hmm     - 0.60  -",
             "  (7) ?     0.037   13.9   0.1  Scaffold4  478876  478690 - hmm     - 0.55  -",
             "  (8) ?     0.071   12.9   0.0  Scaffold4 2314870 2315072 + hmm     - 0.51  -",
             "  (9) ?      0.31   10.8   0.0  Scaffold4 3659261 3659088 - hmm     - 0.56  -",
             " (10) ?         1    9.1   0.0  Scaffold6 1338065 1337981 - hmm     - 0.53  -",
             " (11) ?       1.7    8.3   0.8  Scaffold2 3674449 3674606 + hmm     - 0.59  -"])
        exp = ';rank;;E-value;score;bias;;sequence;start;end;;mdl;trunc;gc;;description\n0;(1);!;2.7e-265;885.0;0.0;;Scaffold3;2966345;2965521;-;hmm;-;0.54;;-\n1;(2);!;1.1e-243;813.4;0.0;;Scaffold3;2965391;2964591;-;hmm;-;0.53;;-\n2;(3);!;1.7e-34;120.9;0.0;;Scaffold4;2315121;2315839;+;hmm;-;0.51;;-\n3;(4);!;1.6e-26;94.6;0.0;;Scaffold6;1337459;1336642;-;hmm;-;0.55;;-\n4;(5);!;5.9e-08;33.1;0.0;;Scaffold4;2314349;2314665;+;hmm;-;0.53;;-\n5;(6);!;0.0059;16.5;0.0;;Scaffold1;5920205;5920315;+;hmm;-;0.60;;-\n6;(7);?;0.037;13.9;0.1;;Scaffold4;478876;478690;-;hmm;-;0.55;;-\n7;(8);?;0.071;12.9;0.0;;Scaffold4;2314870;2315072;+;hmm;-;0.51;;-\n8;(9);?;0.31;10.8;0.0;;Scaffold4;3659261;3659088;-;hmm;-;0.56;;-\n9;(10);?;1;9.1;0.0;;Scaffold6;1338065;1337981;-;hmm;-;0.53;;-\n10;(11);?;1.7;8.3;0.8;;Scaffold2;3674449;3674606;+;hmm;-;0.59;;-\n'
        assert_frame_equal(_str2pd(exp), obs)

if __name__ == '__main__':
    main()
