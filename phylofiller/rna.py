import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_dotplot(fp_dotplot, omit_MEA=True):
    """Read a dot postscript file from Vienna.

    Parameters
    ----------
    fp_dotplot : str
        Filepath to dot.ps file, procuded by e.g. RNAfold
    omit_MEA : bool
        If True, only report base pair probabilities, i.e. upper right
        triangle but NO base pairs used in the Maximum Expected Accurracy
        secondary structure, i.e. lower left triangle (= lbox)

    Returns
    -------
    A tuple of a pandas.DataFrame and a nucleotide sequence as a str.
    The dataframe has a two level index of type int to indicate start and
    end positions of all read base pairs and one column calles "probabilities"
    that holds the actual base pair probabilities, i.e. read data from the ps
    file CONVERTED from sqrt(prob) to prob.
                       probability
        start end
        1     17      0.177101
              19      0.172224
              21      0.277901
              24      0.084927
              38      0.212154
        ...                ...
        38    46      0.194564
        39    44      0.130276
        40    45      0.347719
        41    45      0.090910
        42    46      0.103339
    If omit_MEA == False, the DataFrame holds an additional column "box"
    flagging each probability as either ubox (uppler right triangle) or
    lbox (lower left triangle = MEA).
    """
    data = None
    seq = []

    recordseq = False
    breakNext = False
    peek = None
    with open(fp_dotplot, 'r') as f:
        for lnr, line in enumerate(f.readlines()):
            if line.startswith(') } def'):
                recordseq = False
            if recordseq:
                seq.append(line.strip().replace('\\', ''))
            if line.startswith('/sequence { ('):
                recordseq = True
            if line.startswith('%start of base pair probability data'):
                breakNext = True
                continue
            if breakNext:
                peek = line
                break
        # rnafold and rnaalifold produce different data blocks:
        # RNAfold:
        # %start of base pair probability data
        # 1 8 0.003375481 ubox
        # 1 9 0.003326045 ubox
        #
        # RNAalifold:
        # %start of base pair probability data
        # 0.99 1.00 hsb 1 73 0.999728 ubox
        # 0.60 1.00 hsb 1 73 0.399728 ubox
        #
        # we peek into the first data row and decide on the number of fields
        # which program might be the source
        fields = ['start', 'end', 'probability', 'box']  # default to RNAfold
        if len(peek.strip().split()) > 4:
            fields = ['field1', 'field2', 'hsb', 'start', 'end', 'probability',
                      'box']

        data = pd.read_csv(fp_dotplot, sep=" ", skiprows=lnr+1, header=None,
                           names=fields).iloc[:-3]
        # convert sqrt(p(i,j)) back to p(i,j)
        data['probability'] = data['probability'].apply(lambda x: np.sqrt(x))

        for field in ['start', 'end']:
            data[field] = data[field].astype(int)
        for field in ['probability']:
            data[field] = data[field].astype(float)

        if omit_MEA:
            data = data[data['box'] != 'lbox']
            del data['box']

    return data.set_index(['start', 'end']), ''.join(seq)


def dotPlot(bp_ref, seq_ref, bp_mut=None, seq_mut=None, plotsize=10, title="",
            drawSequence=True, color_ref='darkgreen', color_mut='darkred',
            color_both='yellow'):
    """Produces an RNA 'dot-plot' to visualize base-pair probabilities.

    Parameters
    ----------
    bp_ref : pandas.DataFrame
        Base pair probabilities for a reference folding space.
        The datatype is a pandas DataFrame with a two level index of type
        int that address positions of the base pair and one column that
        holds base pair probabilities, e.g.

                   probability
        start end
        1     17      0.177101
              19      0.172224
              21      0.277901
              24      0.084927
              38      0.212154
        ...                ...
        38    46      0.194564
        39    44      0.130276
        40    45      0.347719
        41    45      0.090910
        42    46      0.103339
    seq_ref : str
        The nucleotide sequence of the RNA molecule.
        (It is used to determine the size of the dot-plot.)
    bp_mut : pandas.DataFrame
        Base pair probabilities of an alternative folding space,
        e.g. a mutated RNA molecule. Datatype is the same as in
        bp_ref.
    seq_mut : str
        The nucleotide sequence of a mutated RNA molecule.
    plotsize : int
        The relative size of the figure, which also depends on the
        seq_ref length.
    title : str
        A title for the dot-plot.
    drawSequence : bool
        seq_ref and/or seq_mut won't we used as tick labels if set
        to False.
    color_ref : str
        The color to plot base-pair probabilities for the reference
        folding space.
    color_mut : str
        The color to plot base-pair probabilities for the mutated
        folding space.
    color_both : str
        The color to plot base-pair probabilities that overlap between
        the reference and the mutated folding space.

    Returns
    -------
    plt.figure
    """
    length = len(seq_ref)
    plotsize = 0.15 * length

    fig, ax = plt.subplots(figsize=(plotsize, plotsize))

    # grid at every 10th base
    for pos in range(10, length, 10):
        ax.axvline(x=pos, color='gray', linewidth=0.5, zorder=-5.0)
        ax.axhline(y=pos, color='gray', linewidth=0.5, zorder=-5.0)

    # diagonal line
    ax.plot([-0.5, length], [-0.5, length], color='black', linewidth=0.5)

    # plot overlay of two base pair probability sets
    for ((start, end), probs) in pd.concat(
            [bp_ref, bp_mut], axis=1).fillna(0).iterrows():
        bpprobA = probs.iloc[0]
        bpprobB = 0
        if probs.shape[0] == 2:
            bpprobB = probs.iloc[1]

        # scale basepair probability with square area -> take square root
        # (same as in Vienna)
        (sizeA, sizeB) = map(lambda x: x**2, (bpprobA, bpprobB))

        def _getCoords(start, end, size):
            return start - (size / 2) - 1, end - (size / 2) - 1
        size_larger, (x_larger, y_larger) = sizeA, _getCoords(
            start, end, sizeA)
        size_smaller, (x_smaller, y_smaller) = sizeB, _getCoords(
            start, end, sizeB)
        color = color_ref
        if bpprobB > bpprobA:
            size_larger, (x_larger, y_larger), \
                size_smaller, (x_smaller, y_smaller) = \
                size_smaller, (x_smaller, y_smaller), \
                size_larger, (x_larger, y_larger)
            color = color_mut
        ax.add_patch(plt.Rectangle((y_larger, x_larger),
                                   size_larger, size_larger, color=color))
        if size_smaller > 0:
            ax.add_patch(plt.Rectangle((y_smaller, x_smaller),
                                       size_smaller, size_smaller,
                                       color=color_both))

    ax.set_xlim((-0.5, length - 0.5))
    ax.set_ylim((-0.5, length - 0.5))

    if seq_mut is None:
        seq_mut = seq_ref
    top_x = list(map(lambda x: '\n'.join(x) if x[0] != x[1] else "\n" + x[0],
                 zip(seq_mut, seq_ref)))
    bottom_x = list(
        map(lambda x: '\n'.join(x) if x[0] != x[1] else x[0] + "\n",
            zip(seq_ref, seq_mut)))
    left_y = list(map(lambda x: ' '.join(x) if x[0] != x[1] else x[0],
                  zip(seq_mut, seq_ref)))
    right_y = list(map(lambda x: ' '.join(x) if x[0] != x[1] else x[0],
                   zip(seq_ref, seq_mut)))
    # sequence as tick labels at all four corners
    if drawSequence:
        ax.set_xticks(range(0, length), )
        ax.set_xticklabels(bottom_x, fontsize=6)
        ax.tick_params(axis=u'both', which=u'both', length=0)

        ax.set_yticks(ax.get_xticks())
        ax.set_yticklabels(left_y, fontsize=6, ha="right")

        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(ax.get_xticks())
        ax2.set_xticklabels(top_x, fontsize=6)
        ax2.tick_params(axis=u'both', which=u'both', length=0)

        ax3 = ax.twinx()
        ax3.set_ylim(ax.get_ylim())
        ax3.set_yticks(ax.get_xticks())
        ax3.set_yticklabels(reversed(right_y), fontsize=6, ha='left')
        ax3.tick_params(axis=u'both', which=u'both', length=0)

    ax.invert_yaxis()

    ax.set_title(title, fontsize=30)

    return fig
