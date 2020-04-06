import yaml
import pandas as pd
import numpy as np
from os.path import join
from os import makedirs
import glob
import sys
import re


def parse_samplesheet(fp_samplesheet):
    #print(fp_samplesheet.split('/')[-1])
    # in a first iteration, open the file, read line by line and determine start
    # of sample information by looking for a line starting with "[Data]".
    # the following lines will be sample information, about lines are header infos.
    row_sampleinformation = None
    row_reads = None
    with open(fp_samplesheet, "r") as f:
        for linenumber, line in enumerate(f.readlines()):
            if line.startswith("[Data]"):
                row_sampleinformation = linenumber+1
            elif line.startswith("[Reads]"):
                row_reads = linenumber+1
    if row_sampleinformation is None:
        raise ValueError("Could not find [Data] line in file '%s'." % fp_samplesheet)
    if row_reads is None:
        raise ValueError("Could not find [Reads] line in file '%s'." % fp_samplesheet)

    header = pd.read_csv(fp_samplesheet, sep=",", nrows=row_reads-2, index_col=0).dropna(axis=1, how="all").dropna(axis=0, how="all")
    #header = header.set_index(header.columns[0])
    header.index = list(map(lambda x: 'header_%s' % x, header.index))
    header = header.dropna(axis=0, how="any")
    header = header.T.reset_index()
    del header['index']

    # a xxx iteration parses sample information via pandas
    ss = pd.read_csv(fp_samplesheet, sep=",", skiprows=row_sampleinformation, dtype={'Sample_Name': str, 'Sample_ID': str, 'spike_entity_id': str})

    # bcl2fasta automatically changes - into _ char in output filenames
    idx_rawilluminainput = ss[pd.notnull(ss['Lane'])].index
    for f in ['Sample_ID', 'Sample_Name', 'Sample_Project']:
        ss.loc[idx_rawilluminainput, f] = ss.loc[idx_rawilluminainput, f].apply(lambda x: x.replace('-', '_') if type(x) != float else x)

    # bcl2fastq uses a S%03i index to address samples.
    # They are numbered as occuring in the samplesheet order starting with 1.
    # However, number is not increased if Sample_ID was already seen.
    uidx = dict()
    for _, sample_id in ss['Sample_ID'].iteritems():
        if sample_id not in uidx:
            uidx[sample_id] = len(uidx) + 1
    ss['s-idx'] = ss['Sample_ID'].apply(lambda x: uidx[x])
    ss['run'] = fp_samplesheet.split('/')[-1].replace('_spike.csv', '')

    # TODO: ensure that sample names do not clash when not considering s-idx!

    # fastq-prefix
    fp_fastqs = []
    for idx, row in ss.iterrows():
        fp_fastq = ''
        if pd.notnull(row['Sample_Project']):
            fp_fastq = row['Sample_Project']
        if pd.notnull(row['Sample_Name']):
            fp_fastq = join(fp_fastq, row['Sample_ID'])
        fp_fastqs.append(join(fp_fastq,
            '%s' % (
                row['Sample_Name'] if pd.notnull(
                    row['Sample_Name']) else row['Sample_ID'])))
    ss['fastq-prefix'] = fp_fastqs

    # remove samples that are marked to be ignored
    if 'spike_ignore_sample' in ss.columns:
        ss = ss[pd.isnull(ss['spike_ignore_sample'])]

    if 'spike_notes' not in ss.columns:
        ss['spike_notes'] = None

    # merge with header information
    if not all([c not in ss.columns for c in header.columns]):
        raise ValueError("Header name conflicts with sample column in '%s'." % fp_samplesheet)
    for c in header.columns:
        ss[c] = header[c].iloc[0]

    return ss


def validate_samplesheet(ss: pd.DataFrame, config, line_offset: int=22, err=sys.stderr):
    """Checks if sample sheet is valid.

    Parameters
    ----------
    ss : pd.DataFrame
        Samplesheet to be validated.
    config : dict from YAML
        Snakemake configuration file holding information about projects.
    line_offset : int
        Default: 22.
        To give user information about problematic lines, we need to go back
        to the file (not the DataFrame) to address the correct line.
    err : IO.stream
        Default: sys.stderr
        Stream onto which warnings are written.

    Returns
    -------
    [str] : List of warnings

    Raises
    ------
    ValueError if errors are found in the sample sheet.
    """

    errors = []
    warnings = []

    # ensure all needed columns are in the table
    exp_columns = {'Lane', 'Sample_ID', 'Sample_Name', 'I7_Index_ID', 'index',
                   'Sample_Project', 'spike_entity_id', 'spike_entity_role'}
    if len(exp_columns - set(ss.columns)) > 0:
        errors.append(
            'Samplesheet is missing column(s): "%s".' %
            '", "'.join(sorted(exp_columns - set(ss.columns))))

    # ensure to only use [A-z0-9_] in identifiers
    allowedChars = re.compile("^[A-z0-9_]*$")
    for field in ['Sample_ID', 'Sample_Name', 'Sample_Project',
                  'spike_entity_id', 'spike_entity_role']:
        if field in ss:
            for idx, x in ss[field].iteritems():
                if pd.notnull(x):
                    if allowedChars.fullmatch(x) is None:
                        errors.append(
                            ('%s in line %i contains a restricted char'
                             'acter: "%s". Only a-z A-Z 0-9 and _ are al'
                             'lowed!') % (field, line_offset+idx, x))

    # ensure Sample_Project is not empty
    if 'Sample_Project' in ss:
        for idx, x in ss['Sample_Project'].iteritems():
            if pd.isnull(x) or x.strip() == "":
                errors.append('Line %i has an empty Sample_Project.' %
                              (line_offset+idx))

    if len(errors) > 0:
        raise ValueError('The following %i errors(s) were found in your sample sheet:\n%s\n' % (len(errors), '\n'.join(['ERROR %i: %s' % (i+1, error) for i, error in enumerate(errors)])))


    # check that sample project is describes in config.yaml
    for prj in ss['Sample_Project'].unique():
        if prj not in config['projects']:
            warnings.append(('Sample_Project "%s" is not described in config.'
                             'yaml. No processing other than demultiplexing w'
                             'ill be applied.') % (prj))

    # check that spike_entity_role is a defined one
    exp_roles = {         'patient',       'father',       'mother',       'sibling',     'healthy',
                 'tumor', 'tumor_patient', 'tumor_father', 'tumor_mother', 'tumor_sibling'}
    for idx, row in ss.iterrows():
        if pd.notnull(row['spike_entity_role']):
            if row['spike_entity_role'] not in exp_roles:
                warnings.append('spike_entity_role "%s" in line %i for Sample_Project "%s" is unknown!' % (row['spike_entity_role'], line_offset+idx, row['Sample_Project']))

    # test that entity name is infix of sample name
    for idx, row in ss.iterrows():
        if pd.notnull(row['spike_entity_id']):
            if row['spike_entity_id'] not in row['Sample_ID']:
                warnings.append('spike_entity_id "%s" is not part of the Sample_ID "%s" in line %i.' % (row['spike_entity_id'], row['Sample_ID'], line_offset+idx))

    # check assumptions about naming schemas per project
    exp_names = {'Keimbahn': re.compile("^KB\d{4}"),
                 'Alps': re.compile("^ALPS")}
    for idx, row in ss.iterrows():
        if row['Sample_Project'] in exp_names:
            if exp_names[row['Sample_Project']].match(row['Sample_ID']) is None:
                warnings.append('Sample_ID "%s" does not follow expected naming schema "%s" in line %i.' % (row['Sample_ID'], exp_names[row['Sample_Project']].pattern, line_offset+idx))

    # check assumptions about name suffices
    exp_suffices = {'Keimbahn': {'patient': {'_c'},
                                 'father': {'_f'},
                                 'mother': {'_m'}},
                    'Alps': {'patient': {''},
                             'father': {'_a', 'a'},
                             'mother': {'_b', 'b'}},
                    'Maus_Hauer': {'healthy': {'_c', 'c', '_n', 'n'},
                                   'tumor': {'_t', 't'}}}
    for idx, row in ss.iterrows():
        if pd.isnull(row['spike_entity_id']):
            continue
        suffix = row['Sample_ID'][len(row['spike_entity_id']):].lower()
        if row['Sample_Project'] in exp_suffices:
            for role in exp_suffices[row['Sample_Project']].keys():
                if (row['spike_entity_role'] == role) and (suffix not in exp_suffices[row['Sample_Project']][role]):
                    warnings.append('Sample_ID "%s" does not match expected spike_entity_role "%s" for Sample_Project "%s" in line %i.' % (row['Sample_ID'], row['spike_entity_role'], row['Sample_Project'], line_offset+idx))

    # check assumptions about barcodes used by specific wet lab members:
    exp_barcodes = {
        'Keimbahn': {
            'A01': 'ATGCCTAA',
            'B01': 'GAATCTGA',
            'C01': 'AACGTGAT',
            'D01': 'CACTTCGA',
            'E01': 'GCCAAGAC',
            'F01': 'GACTAGTA',
            'G01': 'ATTGGCTC',
            'H01': 'GATGAATC'},
        'Alps': {
            'A01': 'ATGCCTAA',
            'B01': 'GAATCTGA',
            'C01': 'AACGTGAT',
            'D01': 'CACTTCGA',
            'E01': 'GCCAAGAC',
            'F01': 'GACTAGTA',
            'G01': 'ATTGGCTC',
            'H01': 'GATGAATC'},
        'Maus_Hauer': {
            'A02': 'AGCAGGAA',
            'B02': 'GAGCTGAA',
            'C02': 'AAACATCG',
            'D02': 'GAGTTAGC',
            'E02': 'CGAACTTA',
            'F02': 'GATAGACA',
            'G02': 'AAGGACAC',
            'H02': 'GACAGTGC'},
    }
    for idx, row in ss.iterrows():
        if row['Sample_Project'] in exp_barcodes:
            if row['I7_Index_ID'] not in exp_barcodes[row['Sample_Project']]:
                warnings.append('Sample_ID "%s" for Sample_Project "%s" in line %i uses unexpected demultiplexing barcode %s: "%s"' % (row['Sample_ID'], row['Sample_Project'], line_offset+idx, row['I7_Index_ID'], row['index']))
            elif (exp_barcodes[row['Sample_Project']][row['I7_Index_ID']] != row['index']):
                warnings.append('Sample_ID "%s" for Sample_Project "%s" in line %i uses unexpected combination of index and I7_index_ID %s: "%s"' % (row['Sample_ID'], row['Sample_Project'], line_offset+idx, row['I7_Index_ID'], row['index']))

    if len(warnings) > 0:
        err.write('The following %i warning(s) are issued by your sample sheet:\n' % len(warnings))
        for i, warning in enumerate(warnings):
            err.write('warning %i: %s\n' % (i+1, warning))

    return warnings


def get_global_samplesheets(dir_samplesheets, config):
    # parse all available sample sheets
    fps_samplesheets = glob.glob('%s*_spike.csv' % dir_samplesheets)

    global_samplesheet = []
    for fp_samplesheet in fps_samplesheets:
        ss = parse_samplesheet(fp_samplesheet)
        global_samplesheet.append(ss)
    if len(global_samplesheet) <= 0:
        raise ValueError("Could not find a single demultiplexing sample sheet in directory '%s'." % dir_samplesheets)
    global_samplesheet = pd.concat(global_samplesheet, sort=False)

    return add_aliassamples(global_samplesheet, config)


def write_samplesheet(fp_output, samplesheet):
    with open(fp_output, 'w') as f:
        # obtain number of data colums, which dictates number of "," in each line
        data_cols = ['Lane',
                     'Sample_ID',
                     'Sample_Name',
                     'Sample_Plate',
                     'Sample_Well',
                     'I7_Index_ID',
                     'index']
        if 'I5_Index_ID' in samplesheet.columns:
            data_cols += ['I5_Index_ID',
                          'index2']
        data_cols += ['Sample_Project',
                      'Description',
                      'spike_notes']

        # header
        f.write('[Header]\n')
        header_cols = []
        for col in sorted(samplesheet.columns):
            if col.startswith('header_'):
                if samplesheet[col].dropna().unique().shape[0] <= 0:
                    continue
                if samplesheet[col].dropna().unique().shape[0] > 1:
                    raise ValueError("Conflicting header information!")
                header_cols.append(col)
        f.write(samplesheet[header_cols].drop_duplicates().rename(columns={col: col[len('header_'):] for col in header_cols}).T.to_csv(header=False))

        # reads & settings
        if 'header_kind_of_run' in samplesheet.columns:
            f.write('\n')
            pattern = re.compile("^(\d+)x(\d+)bp$")
            match = pattern.fullmatch(samplesheet['header_kind_of_run'].dropna().unique()[0])
            MAP_ADAPTERS = {0: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
                            1: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
                            'miseq': 'CTGTCTCTTATACACATCT'}
            if match is not None:
                f.write('[Reads]\n')
                for r in range(int(match.group(1))):
                    f.write('%s\n' % match.group(2))
                f.write('\n')
                f.write('[Settings]\n')
                f.write('ReverseComplement,0\n')
                for r in range(int(match.group(1))):
                    if r > 0:
                        f.write('AdapterRead%i' % (r+1))
                    else:
                        f.write('Adapter')
                    f.write(',%s\n' % MAP_ADAPTERS['miseq' if '_000000000-' in samplesheet['run'].dropna().unique()[0] else r])
            f.write('\n')

        # data
        f.write('[Data]')
        f.write('\n')
        f.write(samplesheet[data_cols].sort_values('Lane').fillna('').to_csv(index=False, float_format='%i'))


def split_samplesheets(samples, config, dry=False):
    """Creates (multiple) samplesheets for bcl2fastq according to barcode length.

    Parameters
    ----------
    samples : pd.DataFrame
        Sample metadata from parsing global samplesheets.
    config : dict
        Snakemakes config objects
    dry : Bool
        Default: False
        If True, only return filepaths without creating any dirs or files.

    Returns
    -------
        List of created samplesheet filepaths.
    """
    if len(samples['run'].unique()) != 1:
        raise ValueError('Not all samples belong to one unique run.')

    results = []
    ss_split = samples.copy()
    ss_split['barcode_len'] = ss_split['index'].fillna("").apply(len)
    split_by = ['barcode_len']
    if ss_split['index'].dropna().shape[0] > 1:
        split_by.append('Lane')
    for i, (grp, ss_barcode) in enumerate(ss_split.sort_values(by=['barcode_len', 'Lane']).groupby(split_by)):
        fp_dir = join(config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['split_demultiplex'], ss_barcode['run'].unique()[0])
        fp_samplesheet = join(fp_dir, 'samplesheet_part_%i.csv' % (i+1))
        results.append(fp_dir)
        if dry is not True:
            makedirs(fp_dir, exist_ok=True)
            write_samplesheet(fp_samplesheet, ss_barcode)

    if dry is True:
        return len(results)
    else:
        return results


def get_role(spike_project, spike_entity_id, spike_entity_role, samplesheets):
    """Returns file path for bam, given project, entity and role (for trio).

    Parameters
    ----------
    spike_project : str
        Name of project, to avoid entity ID clashes across projects.
    spike_entity_id : str
        Entity ID for which role needs to be obtained.
    spike_entity_role : str
        Role of entity ID whose bam filepath shall be returned.
    config : snakemake.config
        Snakemakes config object to obtain file path of sample sheets.

    Returns
    -------
    str: Filepath of bam file for given entity role.
    """
    samples = samplesheets

    if spike_entity_role in ['patient', 'mother', 'father', 'sibling']:
        # edge case: trios shall be computed not for patient, but for e.g. siblings
        # Usecase in Keimbahn project, e.g. KB0164
        # 1) check that no regular sample can be found, because of e.g. suffix _s1
        if samples[(samples['Sample_Project'] == spike_project) &
                   (samples['spike_entity_id'] == spike_entity_id)].shape[0] == 0:
            alt_samples = samples[(samples['Sample_Project'] == spike_project) &
                                  (samples['Sample_ID'] == spike_entity_id) &
                                  (samples['spike_entity_role'] == 'sibling')]
            # 2) test that excatly ONE alternative sample can be found (might be merged across runs/lanes)
            if alt_samples[['Sample_ID', 'Sample_Name', 'Sample_Project', 'spike_entity_id', 'spike_entity_role', 'fastq-prefix']].drop_duplicates().shape[0] != 1:
                raise ValueError('Alternative entity name leads to none or ambiguous sample information!')
            if spike_entity_role == 'patient':
                return alt_samples['fastq-prefix'].unique()[0]
            else:
                return get_role(spike_project, alt_samples['spike_entity_id'].unique()[0], spike_entity_role, samplesheets)
    elif spike_entity_role in ['tumor', 'healthy']:
        # edge case 2: trios might have additional tumor samples (patient, mother, father, siblings are healthy, i.e. germline)
        # Usecase in Keimbahn project, e.g. KB0049
        if samples[(samples['Sample_Project'] == spike_project) &
                   (samples['spike_entity_id'] == spike_entity_id)].shape[0] == 0:
            alt_samples = samples[(samples['Sample_Project'] == spike_project) &
                                  (samples['Sample_ID'] == spike_entity_id) &
                                  (samples['spike_entity_role'].apply(lambda x: x.split('_')[0] if pd.notnull(x) else "") == 'tumor')]
            if alt_samples[['Sample_ID', 'Sample_Name', 'Sample_Project', 'spike_entity_id', 'spike_entity_role', 'fastq-prefix']].drop_duplicates().shape[0] != 1:
                raise ValueError('Alternative entity name leads to none or ambiguous sample information!')
            if spike_entity_role == 'tumor':
                return alt_samples['fastq-prefix'].unique()[0]
            elif spike_entity_role == 'healthy':
                return get_role(spike_project, alt_samples['spike_entity_id'].unique()[0], alt_samples['spike_entity_role'].unique()[0].split('_')[-1], samplesheets)

    # select correct project
    try:
        x = samples[samples['Sample_Project'] == spike_project]
        x.iloc[0]
    except IndexError:
        raise ValueError('Could not find a spike project with name "%s". Available projects are:\n\t%s\n' % (spike_project, '\n\t'.join(sorted(samples['Sample_Project'].unique()))))
    else:
        samples = x

    # select correct entity
    try:
        x = samples[samples['spike_entity_id'] == spike_entity_id]
        x.iloc[0]
    except IndexError:
        raise ValueError('Could not find a spike entity group with name "%s". Available entities for projects "%s" are:\n\t%s\n' % (spike_entity_id, spike_project, '\n\t'.join(sorted(samples['spike_entity_id'].unique()))))
    else:
        samples = x

    # select correct role
    try:
        x = samples[samples['spike_entity_role'] == spike_entity_role]
        x.iloc[0]
    except IndexError:
        raise ValueError('Could not find a role "%s" for spike entity group with name "%s". Available roles are:\n\t%s\n' % (spike_entity_role, spike_entity_id, '\n\t'.join(sorted(samples['spike_entity_role'].unique()))))
    else:
        samples = x

    res = {sample['fastq-prefix'] for idx, sample in samples.iterrows()}
    if len(res) > 1:
        raise ValueError("Stefan, check if use cases can occour with more than one result!\nspike_project: %s, spike_entity_id: %s, spike_entity_role: %s\n%s" % (spike_project, spike_entity_id, spike_entity_role, res))

    return list(res)[0]


def get_species(sample, samplesheets, config):
    # sample can be a single sample ...
    projects = samplesheets[
        (samplesheets['fastq-prefix'] == sample) &
        (samplesheets['is_alias'] != True)]['Sample_Project'].unique()

    # ... or an entity
    if len(projects) == 0:
        projects = samplesheets[
            (samplesheets['Sample_Project'] == sample.split('/')[0]) &
            ((samplesheets['spike_entity_id'] == sample.split('/')[-1]) |
             (samplesheets['Sample_ID'] == sample.split('/')[-1]))]['Sample_Project'].unique()

    if len(projects) > 1:
        raise ValueError("Ambiguous projects: '%s' for sample '%s'" % (projects, sample))

    if projects[0] not in config['projects']:
        raise ValueError('Project "%s" is not specified in config.yaml!' % projects[0])
    if 'species' not in config['projects'][projects[0]]:
        raise ValueError('"species" is not specified for project "%s" in config.yaml!' % projects[0])
    return config['projects'][projects[0]]['species']


def get_reference_genome(sample, samplesheets, config):
    return config['references']['genomes'][get_species(sample, samplesheets, config)]


def get_reference_knowns(sample, samplesheets, config, _key):
    return [k for k in config['references']['knowns'][get_species(sample, samplesheets, config)] if _key in k]


def get_reference_exometrack(sample, samplesheets, config, returnfield='file', debug=False):
    # there are three ways to define the capture kit for a sample:
    # 1. by adding a column "capture_kit" to the input samplesheet and set the cell value to a capture kit name defined in config.yaml
    # 2. by specifing a key "capture_kit" at the project in config.yaml. Value must match one of the defined capture kits in config.yaml
    # 3. by specifing the species for the project, which than select the default kit defined in config.yaml for the species.
    # 1,2,3 are listed in order of their precedence, i.e. of 1 and 3 is defined, 1 wins.

    if 'references' not in config:
        raise ValueError("Key 'references' is not defined in config.yaml")
    if 'capture_kits' not in config['references']:
        raise ValueError('Definition of "capture_kits" in config.yaml is missing.')

    capture_kit = None
    # use-case 1
    if capture_kit is None:
        if 'capture_kit' in samplesheets.columns:
            sample_capture_kits = samplesheets[(samplesheets['fastq-prefix'] == sample) & (samplesheets['is_alias'] != True)]['capture_kit'].dropna().unique()
            if sample_capture_kits.shape[0] > 1:
                raise ValueError("Ambiguous per-sample capture-kit definition in samplesheet for sample '%s': '%s'" % (sample, "', '".join(sample_capture_kits)))
            if sample_capture_kits.shape[0] == 1:
                capture_kit = sample_capture_kits[0]
                if debug:
                    print('per-sample:', sample, capture_kit)

    # use-case 2
    if capture_kit is None:
        sample_project = samplesheets[(samplesheets['fastq-prefix'] == sample) & (samplesheets['is_alias'] != True)]['Sample_Project'].dropna().unique()
        if sample_project.shape[0] > 1:
            raise ValueError("Ambigious projects for '%s': %s" % (sample, "', '".join(sample_project)))
        if sample_project.shape[0] == 0:
            raise ValueError("Missing project definition for sample '%s'." % sample)
        sample_project = sample_project[0]
        if ('projects' in config) and (sample_project in config['projects']) and ('capture_kit' in config['projects'][sample_project]):
            capture_kit = config['projects'][sample_project]['capture_kit']
            if debug:
                print('per-project', sample, sample_project, capture_kit)

    # use-case 3
    if capture_kit is None:
        species = get_species(sample, samplesheets, config)
        for kit_name in config['references']['capture_kits']:
            if 'default_for_species' in config['references']['capture_kits'][kit_name]:
                if config['references']['capture_kits'][kit_name]['default_for_species'] == species:
                    capture_kit = kit_name
                    if debug:
                        print('per-species', sample, kit_name)

    if capture_kit is None:
        raise ValueError("Could not determine capture kit for sample '%s'." % sample)

    # finally, return found capture_kit information
    return config['references']['capture_kits'][capture_kit][returnfield]


def get_reference_varscan_somatic(sample, samplesheets, config):
    return config['references']['varscan_somatic'][get_species(sample, samplesheets, config)]


######## avoid run
def _run2date(run):
    date='%04i/%02i/%02i' % (
        int(run.split('_')[0][:2])+2000,
        int(run.split('_')[0][2:4]),
        int(run.split('_')[0][4:6]))
    return date


def get_bwa_mem_header(sample, samplesheets, config):
    samples = samplesheets[samplesheets['fastq-prefix'] == sample]
    res = []
    for run in sorted(samples['run'].dropna().unique()):
        res.append(' -R "@RG\\tID:%s\\tCN:Department_of_Pediatric_Oncology_Dusseldorf\\tPU:%s\\tDT:%s\\tPL:ILLUMINA\\tLB:%s\\tSM:%s" ' % (
            run,
            run.split('_')[-1][1:],
            _run2date(run),
            get_reference_exometrack(sample, samplesheets, config, returnfield='protocol_name'),
            samples['Sample_ID'].dropna().unique()[0],
        ))
    return "".join(res)


def get_demux_samples(samplesheets, config):
    # ignore samples aliases
    samples = samplesheets[samplesheets['is_alias'] != True]

    # remove samples that stem from per sample fastq sources, like Macrogen sequencing
    samples = samples[pd.notnull(samples['Lane'])]

    return list(samples['run'].unique())


def get_samples(samplesheets, config):
    # only consider samples that have some spike_entity_role defined
    samples_with_role = samplesheets[pd.notnull(samplesheets['spike_entity_role'])]

    samples = []
    for sample, g in samples_with_role.groupby(['Sample_Project', 'fastq-prefix']):
        samples.append({'Sample_Project': sample[0],
                        'sample': sample[1],
                        'spike_entity_id': g['spike_entity_id'].iloc[0]})

    return samples


def get_tumorNormalPairs(samplesheets, config, species=None):
    # only consider samples that have some spike_entity_role defined
    samples_with_role = samplesheets[pd.notnull(samplesheets['spike_entity_role'])]

    pairs = []
    for pair, g in samples_with_role.groupby(['Sample_Project', 'spike_entity_id']):
        # only choose comlete pairs
        if len(set(g['spike_entity_role'].unique()) & {'healthy','tumor'}) == 2:
            if species is not None:
                if get_species(g['fastq-prefix'].iloc[0], samplesheets, config) != species:
                    continue
            pairs.append({'Sample_Project': pair[0],
                          'spike_entity_id': pair[1]})

    # add tumor/normal computations for trio-like projects, i.e. Keimbahn, where special samples stem from tumor tissue and
    # need to be compared to the normal germline (i.e. healthy) samples, e.g. KB0049
    for (project, tumor), g in samples_with_role[samples_with_role['spike_entity_role'].apply(lambda x: x.startswith('tumor_'))].groupby(['Sample_Project', 'Sample_ID']):
        if species is not None:
            if get_species(g['fastq-prefix'].iloc[0], samplesheets, config) != species:
                continue
        pairs.append({'Sample_Project': project,
                      'spike_entity_id': tumor})

    return pairs


def get_genepanels(samplesheets, config, prefix):
    """Returns list of gene panel result files.

    Parameters
    ----------
    samplesheets : pd.DataFrame
        Global samplesheets.
    config : dict
        Snakemakes config object.
    prefix : str
        Filepath to prefix directory.

    Returns
    -------
    [str] : List of filepaths for gene panel results that shall be computed.
    """
    # collect which panels should be computed for which projects
    project_panels = dict()
    if 'projects' not in config:
        raise ValueError('config.yaml does not contain any projects!')
    for project in samplesheets['Sample_Project'].unique():
        if project in config['projects']:
            if (config['projects'][project] is not None) and ('genepanels' in config['projects'][project]):
                project_panels[project] = config['projects'][project]['genepanels']

    # for every sample, check which panels have to be computed
    to_be_created = []
    for project in project_panels.keys():
        for panel in project_panels[project]:
            to_be_created.extend(
                ['%s%s%s/%s.yaml/%s/%s.tsv' % (prefix, config['dirs']['intermediate'], config['stepnames']['genepanel_coverage'], panel, project, sample)
                 for sample
                 in samplesheets[(samplesheets['Sample_Project'] == project) & (samplesheets['is_alias'] != True)]['Sample_ID'].unique()])

    # in addition to the above, also add samples used as aliases
    if 'sample_aliases' in config:
        for sample in config['sample_aliases']:
            if ('roles' in sample) and ('real_id' in sample):
                for role in sample['roles']:
                    if ('Sample_Project' in role) and (role['Sample_Project'] in project_panels):
                        for panel in project_panels[role['Sample_Project']]:
                            if ('Sample_Project' in sample['real_id']) and ('Sample_ID' in sample['real_id']):
                                if samplesheets[(samplesheets['Sample_Project'] == sample['real_id']['Sample_Project']) & (samplesheets['Sample_ID'] == sample['real_id']['Sample_ID'])].shape[0] > 0:
                                    to_be_created.append('%s%s%s/%s.yaml/%s/%s.tsv' % (prefix, config['dirs']['intermediate'], config['stepnames']['genepanel_coverage'], panel, sample['real_id']['Sample_Project'], sample['real_id']['Sample_ID']))

    return to_be_created


def add_aliassamples(samplesheets, config):
    aliases = []
    if (config is not None) and ('sample_aliases' in config):
        for alias in config['sample_aliases']:
            if not (('roles' in alias) and ('real_id' in alias)):
                raise ValueError('Sample alias is not properly defined (missing keys "roles" or "real_id"), check config.yaml file!')
            for role in alias['roles']:
                if not (('Sample_Project' in alias['real_id']) and ('Sample_ID' in alias['real_id'])):
                    raise ValueError('Sample alias is not properly defined (missing keys "Sample_Project" or "Sample_ID"), check config.yaml file!')
                role['fastq-prefix'] = '%s/%s' % (alias['real_id']['Sample_Project'], alias['real_id']['Sample_ID'])
                if 'run' in samplesheets.columns:
                    role['run'] = '+'.join(samplesheets[(samplesheets['Sample_Project'] == alias['real_id']['Sample_Project']) &
                                                        (samplesheets['Sample_ID'] == alias['real_id']['Sample_ID'])]['run'].unique())
                role['is_alias'] = True
                lanes = samplesheets[(samplesheets['Sample_Project'] == alias['real_id']['Sample_Project']) &
                                     (samplesheets['Sample_ID'] == alias['real_id']['Sample_ID'])]['Lane']
                if lanes.shape[0] > 0:
                    role['Lane'] = lanes.iloc[0]
                aliases.append(role)
    if len(aliases) > 0:
        return pd.concat([samplesheets, pd.DataFrame(aliases)], sort=False)
    else:
        samplesheets['is_alias'] = np.nan
    return samplesheets


def get_trios(samplesheets, config):
    trios = []
    for trio, g in samplesheets[pd.notnull(samplesheets['spike_entity_role'])].groupby(['Sample_Project', 'spike_entity_id']):
        # only choose comlete trios, i.e. entities that have at least patient, mother and father
        # there might also be siblings, but we ignore those samples for now.
        if len(set(g['spike_entity_role'].unique()) & {'patient', 'mother', 'father'}) == 3:
            trios.append({'Sample_Project': trio[0],
                          'spike_entity_id': trio[1]})
        # add extra trios for siblings
        if len(set(g['spike_entity_role'].unique()) & {'sibling', 'mother', 'father'}) == 3:
            for sibling, g_sibling in g[g['spike_entity_role'] == 'sibling'].groupby('Sample_ID'):
                trios.append({'Sample_Project': trio[0],
                              'spike_entity_id': sibling})
    return trios


def get_rejoin_input(prefix, sample, direction, samplesheets, config, _type='files'):
    """Returns list of raw input files; either demux results or per sample fastqs.

    Parameters
    ----------
    prefix : str
        Filepath for prefix directory
    sample : str
        fastq-prefix sample name.
    direction : str
        R1 or R2
    samplesheets : pd.DataFrame
        Global demultiplexing sample sheet.
    config : dict
        Global configuration.
    _type : str
        "files" or "dirs"

    Returns
    -------
    Unique set of file paths for input for function rejoin_samples.
    If _type == "dirs" the input for Snakemake is returned, i.e. dirs for demux,
    if _type == "files" the return list contains multiple fastq.gz file paths.
    """
    res = []
    for _, row in samplesheets[(samplesheets['fastq-prefix'] == sample) & (samplesheets['is_alias'] != True)].iterrows():
        if pd.isnull(row['Lane']):
            res.append('%s%s%s%s/%s_%s.fastq.gz' % (prefix, config['dirs']['inputs'], config['dirs']['persamplefastq'], row['run'], row['fastq-prefix'], direction))
        else:
            if _type == 'files':
                res.append('%s%s%s/%s/%s_L%03i_%s_001.fastq.gz' % (prefix, config['dirs']['intermediate'], config['stepnames']['join_demultiplex'], row['run'], row['fastq-prefix'], row['Lane'], direction))
            elif _type == 'dirs':
                res.append('%s%s%s/%s' % (prefix, config['dirs']['intermediate'], config['stepnames']['join_demultiplex'], row['run']))
            else:
                raise ValueError('get_rejoin_input: Unknown function type')
    res = sorted(list(set(res)))
    return res


def get_xenograft_hybridreference(sample, samplesheets, config):
    project = samplesheets[samplesheets['fastq-prefix'] == sample]['Sample_Project'].dropna().unique()
    if len(project) != 1:
        raise ValueError("_get_reference: Sample '%s' has ambiguous or missing project!" % sample)
    return config['projects'][project[0]]['xenograft']


def get_xenograft_stepname(sample, samplesheets, config):
    project = samplesheets[
        (samplesheets['fastq-prefix'] == sample) &
        (samplesheets['is_alias'] != True)]['Sample_Project'].dropna().unique()
    if len(project) != 1:
        sys.stderr.write('%s\n' % (project))
        raise ValueError("get_xenograft_stepname: Sample '%s' has ambiguous or missing project!" % sample)
    if 'xenograft' in config['projects'][project[0]] and config['projects'][project[0]]['xenograft'] != "":
        return config['stepnames']['xenograft_bwa_sampe']
    else:
        return config['stepnames']['trim']


def get_min_coverage(project, config):
    if project in config['projects']:
        if 'min_coverage' in config['projects'][project]:
            return int(config['projects'][project]['min_coverage'])
    return 30


def get_reverse_file(fp_reverse, wildcards, SAMPLESHEETS, config):
    single_paired_end = get_kind_of_run(wildcards, SAMPLESHEETS, config)
    if single_paired_end == 'Unpaired':
        return []
    elif single_paired_end == 'Paired':
        return fp_reverse
    else:
        raise ValueError('get_reverse_file: that should not have happened')


def get_kind_of_run(wildcards, SAMPLESHEETS, config):
    """Depending on the kind of run (single end or paired end) a sample might produce files in sub-directory Paired or Unpaird.
       This function returns the directory name.

    Parameters
    ----------

    Returns
    -------
    str : 'Paired' or 'Unpaired'"""
    kor = SAMPLESHEETS[SAMPLESHEETS['fastq-prefix'] == wildcards.sample]['header_kind_of_run']
    if (kor.dropna().shape[0] <= 0) or ('x' not in kor.dropna().unique()[0]):
        raise ValueError('Abort, since "kind_of_run" is not defined for sample %s. Please identify the matching samplesheet and add the header field "kind_of_run" which should hold the information if the run is either paired end (2xYYY) or single end (1xYYY), were YYY is the number of sequenced nucleotids.' % wildcards.sample)
    single_paired_end = kor.unique()[0].split('x')[0]
    if single_paired_end == '1':
        return 'Unpaired'
    elif single_paired_end == '2':
        return 'Paired'
    else:
        raise ValueError('The information in header field "kind_of_run" in the sample sheet for sample "%s" is not properly formatted. It should look like "2xYYY" were YYY is the length of the sequences.' % wildcards.sample)
