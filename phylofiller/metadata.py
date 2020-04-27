import pandas as pd

COL_ORGANISM = 'organism'
COL_FP_ASSEMBLY = 'fp_assembly'
COL_AUGUSTUS_REF_SPECIES = 'augustus_reference_species'

CFG_KEY_AUGUSTUS = 'augustus'
CFG_KEY_AUGUSTUS_REFERENCE_SPECIES = 'default_reference_species'


def read_metadata(fp_metadata: str, fp_assemblyprefix: str = None, skip_file_exists_test=False):
    """.

    Parameters
    ----------
    fp_metadata : str
    fp_assemblyprefix : str

    Returns
    -------

    Raises
    ------
    Diverse erros if misconfigurations are present.
    """
    RESERVED_COL_NAMES = ['__line_number', '__file_exists']

    meta = pd.read_csv(fp_metadata, sep="\t", index_col=0)

    conflicting_col_names = set(RESERVED_COL_NAMES) & set(meta.columns)
    if len(conflicting_col_names) > 0:
        raise ValueError(
            'Conflicting column name(s) in your metadata: %s please rename!' %
            ', '.join(map(lambda x: '"%s"' % x, conflicting_col_names)))

    if meta.index.name != COL_ORGANISM:
        raise ValueError(
            'Header of first column must be "%s"' % COL_ORGANISM)

    # temporarily add line numbers to metadata
    meta['__line_number'] = range(2, meta.shape[0]+2)

    # check if organisms have been defined multiple times
    ambig_organisms = meta.loc[meta.index.value_counts() > 1][
        '__line_number'].reset_index().groupby(COL_ORGANISM).apply(
            lambda x: x['__line_number'].values)
    if ambig_organisms.shape[0] > 0:
        raise ValueError(
            '%i organisms have been re-defined in your '
            'metadata. Please fix!\n%s' %
            (ambig_organisms.shape[0], str(ambig_organisms)))

    # add assembly prefix if given
    if fp_assemblyprefix is not None:
        meta[COL_FP_ASSEMBLY] = meta[COL_FP_ASSEMBLY].apply(
            lambda x: os.path.abspath('%s%s' % (fp_assemblyprefix, x)))

    # check if given assembly file paths exists
    if skip_file_exists_test is False:
        meta['__file_exists'] = meta[COL_FP_ASSEMBLY].apply(
            lambda x: os.path.exists(x))
        if meta[~meta['__file_exists']].shape[0] > 0:
            raise ValueError(
                "Cannot find files for following assemblies:\n%s" %
                str(meta[~meta['__file_exists']]))

    del meta['__line_number']

    return meta


def get_augustus_reference_species(organism: str, meta: pd.DataFrame,
                                   prj_config: dict=None):
    """Return reference species for organism as used in augustus.

    Species can be defined in two ways:
    a) project wide in config.yaml
    b) organism specific via metadata column

    Parameters
    ----------
    organism : str
        Organism for which reference species shall be returned.
    meta : pd.DataFrame
        Organism metadata table.
    prj_config : dict
        Project related snakemake configuration object.
    """
    ref_species = None

    # if defined in config.yaml, obtain project wise default species
    if prj_config is not None:
        ref_species = prj_config.get(CFG_KEY_AUGUSTUS, None).get(
            CFG_KEY_AUGUSTUS_REFERENCE_SPECIES, None)

    if (COL_AUGUSTUS_REF_SPECIES in meta.columns) and (organism in meta.index):
        ref_species = meta.loc[organism, COL_AUGUSTUS_REF_SPECIES]

    if ref_species is None:
        raise ValueError(
            "Could not determine reference species for organism '%s'.\n"
            "Please either define as '%s: %s: xxx' key for your project in "
            "config.yaml\nor as column '%s' in your project metadata file." % (
                organism, CFG_KEY_AUGUSTUS, CFG_KEY_AUGUSTUS_REFERENCE_SPECIES,
                COL_AUGUSTUS_REF_SPECIES))
    else:
        return ref_species
