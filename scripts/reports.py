from os.path import join, exists, dirname, basename
from os import makedirs
import sys
import pandas as pd
from glob import glob
import seaborn as sns
import numpy as np
from scipy import stats
import xlsxwriter
import matplotlib.pyplot as plt
from scripts.parse_samplesheet import get_min_coverage, get_role, add_aliassamples, get_species
from scripts.snupy import check_snupy_status
import json
import datetime
import getpass
import socket
import requests
from requests.auth import HTTPBasicAuth
import urllib3
import yaml
import pickle
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
plt.switch_backend('Agg')


RESULT_NOT_PRESENT = -5


def report_undertermined_filesizes(fp_filesizes, fp_output, fp_error,
                                   zscorethreshold=1):
    # read all data
    fps_sizes = glob(join(dirname(fp_filesizes), '*.txt'))
    pds_sizes = []
    for fp_size in fps_sizes:
        data = pd.read_csv(
            fp_size, sep="\t", names=["filesize", "filename", "status"],
            index_col=1)
        # mark given read as isme=True while all other data in the dir
        # are isme=False
        data['isme'] = fp_filesizes in fp_size
        data['filesize'] /= 1024**3
        pds_sizes.append(data)
    pd_sizes = pd.concat(pds_sizes)

    # compute z-score against non-bad known runs
    pd_sizes['z-score'] = np.nan
    idx_nonbad = pd_sizes[pd_sizes['status'] != 'bad'].index
    pd_sizes.loc[idx_nonbad, 'z-score'] = stats.zscore(
        pd_sizes.loc[idx_nonbad, 'filesize'])

    # plot figure
    fig = plt.figure()
    ax = sns.distplot(
        pd_sizes[(pd_sizes['isme'] == np.False_) &
                 (pd_sizes['status'] != 'bad')]['filesize'],
        kde=False, rug=False, color="black", label='known runs')
    ax = sns.distplot(
        pd_sizes[(pd_sizes['isme'] == np.False_) &
                 (pd_sizes['status'] == 'bad')]['filesize'],
        kde=False, rug=False, color="red", label='bad runs')
    ax = sns.distplot(
        pd_sizes[pd_sizes['isme'] == np.True_]['filesize'],
        kde=False, rug=True, color="green", label='this run')
    _ = ax.set_ylabel('number of files')
    _ = ax.set_xlabel('file-size in GB')
    ax.set_title('run %s' % basename(fp_filesizes)[:-4])
    ax.legend()

    # raise error if current run contains surprisingly large undetermined
    # filesize
    if pd_sizes[(pd_sizes['isme'] == np.True_) &
                (pd_sizes['status'] == 'unknown')]['z-score'].max() > zscorethreshold:
        ax.set_title('ERROR: %s' % ax.get_title())
        fig.savefig(fp_error, bbox_inches='tight')
        raise ValueError(
            ("Compared to known historic runs, your run contains surprisingly "
             "(z-score > %f) large file(s) of undetermined reads. You will find"
             " an supporting image at '%s'. Please do the following things:\n"
             "1. discuss with lab personal about the quality of the run.\n"
             "2. should you decide to keep going with this run, mark file "
             "status (3rd column) in file '%s' as 'good'.\n"
             "3. for future automatic considerations, mark file status (3rd "
             "column) as 'bad' if you have decided to abort processing due to"
             " too low quality (z-score kind of averages about known values)."
             ) % (zscorethreshold, fp_error, fp_filesizes))
    else:
        fig.savefig(fp_output, bbox_inches='tight')


def report_exome_coverage(
    fps_sample, fp_plot,
    min_coverage=30, min_targets=80, coverage_cutoff=200):
    """Creates an exome coverage plot for multiple samples.

    Parameters
    ----------
    fps_sample : [str]
        A list of file-paths with coverage data in csv format.
    fp_plot : str
        Filepath of output graph.
    min_coverage : int
        Default: 30.
        An arbitraty threshold of minimal coverage that we expect.
        A vertical dashed line is drawn at this value.
    min_targets : float
        Default: 80.
        An arbitraty threshold of minimal targets that we expect to be covered.
        A horizontal dashed line is drawn at this value.
    coverage_cutoff : float
        Default: 200.
        Rightmost coverage cut-off value where X-axis is limited.

    Raises
    ------
    ValueError : If one of the sample's coverage falls below expected
    thresholds.
    """
    # Usually we aim for a 30X coverage on 80% of the sites.

    fig, ax = plt.subplots()
    ax.axhline(y=min_targets, xmin=0, xmax=coverage_cutoff, color='gray',
               linestyle='--')
    ax.axvline(x=min_coverage, ymin=0, ymax=100, color='gray', linestyle='--')

    samples_below_coverage_threshold = []
    for fp_sample in fps_sample:
        coverage = pd.read_csv(fp_sample, sep="\t")
        samplename = fp_sample.split('/')[-1].split('.')[0]
        linewidth = 1
        if coverage[coverage['#coverage'] == min_coverage]['percent_cumulative'].min() < min_targets:
            linewidth = 4
            samples_below_coverage_threshold.append(samplename)
        ax.plot(coverage['#coverage'],
                coverage['percent_cumulative'],
                label=samplename,
                linewidth=linewidth)

    ax.set_xlim((0, coverage_cutoff))
    ax.set_xlabel('Read Coverage')
    ax.set_ylabel('Targeted Exome Bases')
    ax.legend()

    if len(samples_below_coverage_threshold) > 0:
        fp_plot = fp_plot.replace('.pdf', '.error.pdf')

    fig.savefig(fp_plot, bbox_inches='tight')

    if len(samples_below_coverage_threshold) > 0:
        raise ValueError(
            "The following %i sample(s) have coverage below expected "
            "thresholds. Please discuss with project PIs on how to proceed. "
            "Maybe, samples need to be re-sequenced.\n\t%s\nYou will find more"
            " information in the generated coverage plot '%s'." % (
                len(samples_below_coverage_threshold),
                '\n\t'.join(samples_below_coverage_threshold),
                fp_plot))



ACTION_PROGRAMS = [
    {'action': 'background',
     'program': 'GATK',
     'fileending_snupy_extract': '.snp_indel.gatk',
     'fileending_spike_calls': '.gatk.snp_indel.vcf',
     'stepname_spike_calls': 'gatk_CombineVariants',
    },
    {'action': 'background',
     'program': 'Platypus',
     'fileending_snupy_extract': '.indel.ptp',
     'fileending_spike_calls': '.ptp.annotated.filtered.indels.vcf',
     'stepname_spike_calls': 'platypus_filtered',
    },
    {'action': 'tumornormal',
     'program': 'Varscan',
     'fileending_snupy_extract': '.somatic.varscan',
     'fileending_spike_calls':
        {'homo sapiens': '.snp.somatic_germline.vcf',
         'mus musculus': '.indel_snp.vcf'},
     'stepname_spike_calls': 'merge_somatic',
    },
    {'action': 'tumornormal',
     'program': 'Mutect',
     'fileending_snupy_extract': '.somatic.mutect',
     'fileending_spike_calls': '.all_calls.vcf',
     'stepname_spike_calls': 'mutect',
    },
    {'action': 'tumornormal',
     'program': 'Excavator2',
     'fileending_snupy_extract': '.somatic.cnv.excavator2',
     'fileending_spike_calls': '.vcf',
     'stepname_spike_calls': 'excavator_somatic',
    },
    {'action': 'trio',
     'program': 'Varscan\ndenovo',
     'fileending_snupy_extract': '.denovo.varscan',
     'fileending_spike_calls': '.var2denovo.vcf',
     'stepname_spike_calls': 'writing_headers',
    },
    {'action': 'trio',
     'program': 'Excavator2',
     'fileending_snupy_extract': '.trio.cnv.excavator2',
     'fileending_spike_calls': '.vcf',
     'stepname_spike_calls': 'excavator_trio',
    },
]


def _get_statusdata_demultiplex(samplesheets, prefix, config):
    demux_yields = []
    for flowcell in samplesheets['run'].unique():
        fp_yielddata = '%s%s%s/Data/%s.yield_data.csv' % (prefix, config['dirs']['intermediate'], config['stepnames']['yield_report'], flowcell)
        if exists(fp_yielddata):
            demux_yields.append(
                pd.read_csv(fp_yielddata, sep="\t").rename(columns={'Project': 'Sample_Project', 'Sample': 'Sample_ID', 'Yield': 'yield'})) #.set_index(['Project', 'Lane', 'Sample', 'Barcode sequence'])
    if len(demux_yields) <= 0:
        return pd.DataFrame()
    demux_yields = add_aliassamples(pd.concat(demux_yields, axis=0), config)

    # map yields of original sampels to aliases
    for idx, row in demux_yields[demux_yields['is_alias'] == True].iterrows():
        orig = demux_yields[(demux_yields['Sample_Project'] == row['fastq-prefix'].split('/')[0]) & (demux_yields['Sample_ID'] == row['fastq-prefix'].split('/')[1])]['yield']
        if orig.shape[0] > 0:
            demux_yields.loc[idx, 'yield'] = orig.sum()
    demux_yields = demux_yields.dropna(subset=['yield'])

    return pd.DataFrame(demux_yields).groupby(['Sample_Project', 'Sample_ID'])['yield'].sum()


def _get_statusdata_coverage(samplesheets, prefix, config, min_targets=80):
    coverages = []
    for (sample_project, sample_id), meta in samplesheets.groupby(['Sample_Project', 'Sample_ID']):
        role_sample_project, role_sample_id = sample_project, sample_id
        if (meta['is_alias'] == True).any():
            role_sample_project, role_sample_id = get_role(sample_project, meta['spike_entity_id'].unique()[0], meta['spike_entity_role'].unique()[0], samplesheets).split('/')
        fp_coverage = join(prefix, config['dirs']['intermediate'], config['stepnames']['exome_coverage'], role_sample_project, '%s.exome_coverage.csv' % role_sample_id)
        if exists(fp_coverage):
            coverage = pd.read_csv(fp_coverage, sep="\t")
            if coverage.shape[0] > 0:
                coverages.append({
                    'Sample_Project': sample_project,
                    'Sample_ID': sample_id,
                    'coverage': coverage.loc[coverage['percent_cumulative'].apply(lambda x: abs(x-min_targets)).idxmin(), '#coverage']})
    if len(coverages) <= 0:
        return pd.DataFrame()
    return pd.DataFrame(coverages).set_index(['Sample_Project', 'Sample_ID'])['coverage']


def _isKnownDuo(sample_project, spike_entity_id, config):
    """Checks if trio is a known duo, i.e. missing samples won't be available in the future.

    Parameters
    ----------
    sample_project : str
    spike_entity_id : str
    config : dict()
        Snakemake configuration.

    Returns
    -------
    Boolean: True, if spike_entity_id is in config list of known duos for given project.
    False, otherwise.
    """
    if 'projects' in config:
        if sample_project in config['projects']:
            if 'known_duos' in config['projects'][sample_project]:
                if spike_entity_id in config['projects'][sample_project]['known_duos']:
                    return True
    return False


def _get_statusdata_snupyextracted(samplesheets, prefix, snupy_instance, config):
    results = []
    for sample_project, meta in samplesheets.groupby('Sample_Project'):
        # project in config file is not properly configure for snupy!
        if config['projects'].get(sample_project, None) is None:
            continue
        if config['projects'][sample_project].get('snupy', None) is None:
            continue
        if config['projects'][sample_project]['snupy'][snupy_instance].get('project_id', None) is None:
            continue

        r = requests.get('%s/experiments/%s.json' % (config['credentials']['snupy'][snupy_instance]['host'], config['projects'][sample_project]['snupy'][snupy_instance]['project_id']),
            auth=HTTPBasicAuth(config['credentials']['snupy'][snupy_instance]['username'], config['credentials']['snupy'][snupy_instance]['password']),
            verify=False)
        check_snupy_status(r)
        samples = [sample['name'] for sample in r.json()['samples']]

        for sample_id, meta_sample in meta.groupby('Sample_ID'):
            for file_ending, action, program in [(ap['fileending_snupy_extract'], ap['action'], ap['program']) for ap in ACTION_PROGRAMS]:
                # in some cases "sample name" hold spike_entity_id, in others Sample_ID
                entity = sample_id
                runs = '+'.join(sorted(meta_sample['run'].unique()))
                if (action == 'trio'):
                    if meta_sample['spike_entity_role'].unique()[0] == 'patient':
                        entity = meta_sample['spike_entity_id'].iloc[0]
                        runs = '+'.join(sorted(samplesheets[samplesheets['spike_entity_id'] == meta_sample['spike_entity_id'].iloc[0]]['run'].unique()))
                if (action == 'tumornormal'):
                    if meta_sample['spike_entity_role'].unique()[0] == 'tumor':
                        entity = meta_sample['spike_entity_id'].iloc[0]
                        runs = '+'.join(sorted(samplesheets[samplesheets['spike_entity_id'] == meta_sample['spike_entity_id'].iloc[0]]['run'].unique()))
                name = '%s_%s/%s%s' % (runs, sample_project, entity, file_ending)

                if (sample_project in config['projects']) and (pd.notnull(meta_sample['spike_entity_role'].iloc[0])):
                    if ((action == 'trio') and (meta_sample['spike_entity_role'].iloc[0] in ['patient', 'sibling']) and (not _isKnownDuo(sample_project, meta_sample['spike_entity_id'].iloc[0], config))) or\
                       ((action == 'background')) or\
                       ((action == 'tumornormal') and (meta_sample['spike_entity_role'].iloc[0].startswith('tumor'))):
                        results.append({
                            'Sample_Project': sample_project,
                            'Sample_ID': sample_id,
                            'action': action,
                            'program': program,
                            'status': name in samples,
                            'snupy_sample_name': name
                        })
    if len(results) <= 0:
        return pd.DataFrame()
    return pd.DataFrame(results).set_index(['Sample_Project', 'Sample_ID', 'action', 'program'])


def _get_statusdata_numberpassingcalls(samplesheets, prefix, config, RESULT_NOT_PRESENT, verbose=sys.stderr):
    results = []

    # leave out samples aliases
    for (sample_project, spike_entity_id, spike_entity_role, fastq_prefix), meta in samplesheets[samplesheets['is_alias'] != True].fillna('not defined').groupby(['Sample_Project', 'spike_entity_id', 'spike_entity_role', 'fastq-prefix']):
        def _get_fileending(file_ending, fastq_prefix, samplesheets, config):
            if isinstance(file_ending, dict):
                return file_ending[get_species(fastq_prefix, samplesheets, config)]
            else:
                return file_ending
        for ap in ACTION_PROGRAMS:
            fp_vcf = None
            if (ap['action'] == 'background') and pd.notnull(spike_entity_role):
                if (ap['program'] == 'GATK'):
                    fp_vcf = '%s%s%s/%s%s' % (prefix, config['dirs']['intermediate'], config['stepnames'][ap['stepname_spike_calls']], fastq_prefix, _get_fileending(ap['fileending_spike_calls'], fastq_prefix, meta, config))
                elif (ap['program'] == 'Platypus'):
                    fp_vcf = '%s%s%s/%s%s' % (prefix, config['dirs']['intermediate'], config['stepnames'][ap['stepname_spike_calls']], fastq_prefix, _get_fileending(ap['fileending_spike_calls'], fastq_prefix, meta, config))
            elif (ap['action'] == 'tumornormal'):
                for (alias_sample_project, alias_spike_entity_role, alias_sample_id), alias_meta in samplesheets[(samplesheets['fastq-prefix'] == fastq_prefix) & (samplesheets['spike_entity_role'].apply(lambda x: x.split('_')[0] if pd.notnull(x) else x).isin(['tumor']))].groupby(['Sample_Project', 'spike_entity_role', 'Sample_ID']):
                    # for Keimbahn, the tumor sample needs to include the name of the original sample ID
                    instance_id = '%s/%s' % (alias_sample_project, alias_sample_id)
                    if alias_spike_entity_role == 'tumor':
                        # for Maus_Hauer, the filename holds the entity name, but not the Sample ID
                        instance_id = '%s/%s' % (sample_project, spike_entity_id)
                    if (alias_spike_entity_role.split('_')[0] in set(['tumor'])):
                        if (ap['program'] == 'Varscan'):
                            fp_vcf = '%s%s%s/%s%s' % (prefix, config['dirs']['intermediate'], config['stepnames'][ap['stepname_spike_calls']], instance_id, _get_fileending(ap['fileending_spike_calls'], fastq_prefix, meta, config))
                        elif (ap['program'] == 'Mutect'):
                            fp_vcf = '%s%s%s/%s%s' % (prefix, config['dirs']['intermediate'], config['stepnames'][ap['stepname_spike_calls']], instance_id, _get_fileending(ap['fileending_spike_calls'], fastq_prefix, meta, config))
                        elif (ap['program'] == 'Excavator2'):
                            fp_vcf = '%s%s%s/%s/Results/%s/EXCAVATORRegionCall_%s%s' % (prefix, config['dirs']['intermediate'], config['stepnames'][ap['stepname_spike_calls']], instance_id, fastq_prefix.split('/')[-1], fastq_prefix.split('/')[-1], _get_fileending(ap['fileending_spike_calls'], fastq_prefix, meta, config))
            elif (ap['action'] == 'trio'):
                for (alias_sample_project, alias_spike_entity_role, alias_sample_id, alias_spike_entity_id), alias_meta in samplesheets[(samplesheets['fastq-prefix'] == fastq_prefix) & (samplesheets['spike_entity_role'].isin(['patient', 'sibling']))].groupby(['Sample_Project', 'spike_entity_role', 'Sample_ID', 'spike_entity_id']):
                    # Trios are a more complicated case, since by default the result name is given by the
                    # spike_entity_id, but if computed for siblings, the name is given by the fastq-prefix
                    if (ap['program'] == 'Varscan\ndenovo'):
                        if (alias_spike_entity_role in set(['patient'])):
                            fp_vcf = '%s%s%s/%s/%s%s' % (prefix, config['dirs']['intermediate'], config['stepnames'][ap['stepname_spike_calls']], alias_sample_project, alias_spike_entity_id, _get_fileending(ap['fileending_spike_calls'], fastq_prefix, meta, config))
                        elif (alias_spike_entity_role in set(['sibling'])):
                            fp_vcf = '%s%s%s/%s%s' % (prefix, config['dirs']['intermediate'], config['stepnames'][ap['stepname_spike_calls']], fastq_prefix, _get_fileending(ap['fileending_spike_calls'], fastq_prefix, meta, config))
                    elif (ap['program'] == 'Excavator2'):
                        if (alias_spike_entity_role in set(['patient'])):
                            fp_vcf = '%s%s%s/%s/%s/Results/%s/EXCAVATORRegionCall_%s%s' % (prefix, config['dirs']['intermediate'], config['stepnames'][ap['stepname_spike_calls']], alias_sample_project, alias_spike_entity_id, fastq_prefix.split('/')[-1], fastq_prefix.split('/')[-1], _get_fileending(ap['fileending_spike_calls'], fastq_prefix, meta, config))
                        elif (alias_spike_entity_role in set(['sibling'])):
                            fp_vcf = '%s%s%s/%s/Results/%s/EXCAVATORRegionCall_%s%s' % (prefix, config['dirs']['intermediate'], config['stepnames'][ap['stepname_spike_calls']], fastq_prefix, fastq_prefix.split('/')[-1], fastq_prefix.split('/')[-1], _get_fileending(ap['fileending_spike_calls'], fastq_prefix, meta, config))
                    # remove entry, if it is known (config.yaml) that this trio is incomplete
                    if (spike_entity_role == 'patient') and (spike_entity_id in config.get('projects', []).get(sample_project, []).get('known_duos', [])):
                        fp_vcf = None

            results.append({
                    'Sample_Project': sample_project,
                    'Sample_ID': fastq_prefix.split('/')[-1],
                    'action': ap['action'],
                    'program': ap['program'],
                    'fp_calls': fp_vcf,
            })

    status = 0
    num_status = 20
    if verbose is not None:
        print('of %i: ' % num_status, file=verbose, end="")
    for i, res in enumerate(results):
        if (verbose is not None) and int(i % (len(results) / num_status)) == 0:
            status+=1
            print('%i ' % status, file=verbose, end="")
        nr_calls = RESULT_NOT_PRESENT
        if (res['fp_calls'] is not None) and exists(res['fp_calls']):
            try:
                if res['program'] == 'Varscan':
                    nr_calls = pd.read_csv(res['fp_calls'], comment='#', sep="\t", dtype=str, header=None, usecols=[7], squeeze=True).apply(lambda x: ';SS=2;' in x).sum()
                else:
                    nr_calls = pd.read_csv(res['fp_calls'], comment='#', sep="\t", dtype=str, header=None, usecols=[6], squeeze=True).value_counts()['PASS']
            except pd.io.common.EmptyDataError:
                nr_calls = 0
        res['number_calls'] = nr_calls
    if verbose is not None:
        print('done.', file=verbose)

    if len(results) <= 0:
        return pd.DataFrame()
    results = pd.DataFrame(results)
    results = results[pd.notnull(results['fp_calls'])].set_index(['Sample_Project', 'Sample_ID', 'action', 'program'])['number_calls']

    # add alias sample results
    for (sample_project, spike_entity_id, spike_entity_role, fastq_prefix), meta in samplesheets[samplesheets['is_alias'] == True].groupby(['Sample_Project', 'spike_entity_id', 'spike_entity_role', 'fastq-prefix']):
        for (_, _, action, program), row in results.loc[fastq_prefix.split('/')[0], fastq_prefix.split('/')[-1], :].iteritems():
            results.loc[sample_project, meta['Sample_ID'].unique()[0], action, program] = row

    # remove samples, that don't have their own role, but were used for aliases
    for (sample_project, sample_id), _ in samplesheets[pd.isnull(samplesheets['spike_entity_role'])].groupby(['Sample_Project', 'Sample_ID']):
        idx_to_drop = results.loc[sample_project, sample_id, ['tumornormal', 'trio'], :].index
        if len(idx_to_drop) > 0:
            results.drop(index=idx_to_drop, inplace=True)

    return results


def _get_genepanel_data(samplesheets, prefix, config):
    results = []
    columns = ['Sample_Project', 'Sample_ID', 'genepanel', 'gene']

    # leave out samples aliases
    for (sample_project, spike_entity_id, spike_entity_role, fastq_prefix), meta in samplesheets[samplesheets['is_alias'] != True].fillna('not defined').groupby(['Sample_Project', 'spike_entity_id', 'spike_entity_role', 'fastq-prefix']):
        #print(sample_project, spike_entity_id, spike_entity_role, fastq_prefix)
        for file in glob('%s%s%s/*/%s.tsv' % (prefix, config['dirs']['intermediate'], config['stepnames']['genepanel_coverage'], fastq_prefix)):
            #print("\t", file)
            coverage = pd.read_csv(file, sep="\t")

            parts = file.split('/')
            # determine genepanel name, project and sample_id from filename
            coverage['Sample_Project'] = sample_project
            coverage['Sample_ID'] = meta['Sample_ID'].unique()[0]
            coverage['genepanel'] = parts[-3][:-5]
            coverage = coverage.set_index(columns)

            results.append(coverage)
    if len(results) > 0:
        results = pd.concat(results).sort_values(by=columns)
    else:
        results = pd.DataFrame(columns=columns)

    # add alias sample results
    for (sample_project, spike_entity_id, spike_entity_role, fastq_prefix), meta in samplesheets[samplesheets['is_alias'] == True].groupby(['Sample_Project', 'spike_entity_id', 'spike_entity_role', 'fastq-prefix']):
        for (_, _, action, program), row in results.loc[fastq_prefix.split('/')[0], fastq_prefix.split('/')[-1], :].iterrows():
            results.loc[sample_project, meta['Sample_ID'].unique()[0], action, program] = row

    return results


def get_status_data(samplesheets, config, snupy_instance, prefix=None, verbose=sys.stderr):
    """
    Parameters
    ----------
    samplesheets : pd.DataFrame
        The global samplesheets.
    config : dict()
        Snakemake configuration object.
    prefix : str
        Default: None, i.e. config['dirs']['prefix'] is used.
        Filepath to spike main directory.
    verbose : StringIO
        Default: sys.stderr
        If not None: print verbose information.

    Returns
    -------
    4-tuple: (data_yields, data_coverage, data_snupy, data_calls)
    """
    global RESULT_NOT_PRESENT
    NUMSTEPS = 6

    if prefix is None:
        prefix = config['dirs']['prefix']
    if verbose is not None:
        print("Creating report", file=verbose)
    # obtain data
    if verbose is not None:
        print("1/%i) gathering demuliplexing yields: ..." % NUMSTEPS, file=verbose, end="")
    data_yields = _get_statusdata_demultiplex(samplesheets, prefix, config)
    if verbose is not None:
        print(" done.\n2/%i) gathering coverage: ..." % NUMSTEPS, file=verbose, end="")
    data_coverage = _get_statusdata_coverage(samplesheets, prefix, config)
    if verbose is not None:
        print(" done.\n3/%i) gathering snupy extraction status: ..." % NUMSTEPS, file=verbose, end="")
    data_snupy = _get_statusdata_snupyextracted(samplesheets, prefix, snupy_instance, config)
    if verbose is not None:
        print(" done.\n4/%i) gathering number of PASSing calls: ..." % NUMSTEPS, file=verbose, end="")
    data_calls = _get_statusdata_numberpassingcalls(samplesheets, prefix, config, RESULT_NOT_PRESENT, verbose=verbose)
    if verbose is not None:
        print(" done.\n5/%i) gathering gene coverage: ..." % NUMSTEPS, file=verbose, end="")
    data_genepanels = _get_genepanel_data(samplesheets, prefix, config)
    if verbose is not None:
        print("done.\n6/%i) generating Excel output: ..." % NUMSTEPS, file=verbose, end="")

    return (data_yields, data_coverage, data_snupy, data_calls, data_genepanels)


def write_status_update(data, filename, samplesheets, config, offset_rows=0, offset_cols=0, min_yield=5.0, verbose=sys.stderr):
    """
    Parameters
    ----------
    data : (pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame)
        yields, coverage, snupy, calls. Result of function get_status_data.
    filename : str
        Filepath to output Excel file.
    samplesheets : pd.DataFrame
        The global samplesheets.
    config : dict()
        Snakemake configuration object.
    offset_rows : int
        Default: 0
        Number if rows to leave blank on top.
    offset_cols : int
        Default: 0
        Number if columns to leave blank on the left.
    min_yield : float
        Default: 5.0
        Threshold when to color yield falling below this value in red.
        Note: I don't know what a good default looks like :-/
    verbose : StringIO
        Default: sys.stderr
        If not None: print verbose information.
    """
    global RESULT_NOT_PRESENT

    # for debugging purposes
    pickle.dump(data, open('%s.datadump' % filename, 'wb'))

    data_yields, data_coverage, data_snupy, data_calls, data_genepanels = data

    # start creating the Excel result
    workbook = xlsxwriter.Workbook(filename)
    worksheet = workbook.add_worksheet()

    format_good = workbook.add_format({'bg_color': '#ccffcc'})
    format_bad = workbook.add_format({'bg_color': '#ffcccc'})

    # date information
    format_info = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center',
        'font_size': 9})
    info_username = getpass.getuser()
    info_now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    info_machine = socket.gethostname()
    worksheet.merge_range(offset_rows, offset_cols, offset_rows+1, offset_cols+3, ('status report created\nat %s\nby %s\non %s' % (info_now, info_username, info_machine)),format_info)

    gene_order = []
    if data_genepanels.shape[0] > 0:
        for panel in sorted(data_genepanels.index.get_level_values('genepanel').unique()):
            for gene in sorted(data_genepanels.loc(axis=0)[:, :, panel, :].index.get_level_values('gene').unique()):
                gene_order.append((panel, gene))

    # header action
    format_action = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center',
        'bold': True})
    aps = pd.Series([ap['action'] for ap in ACTION_PROGRAMS]).to_frame()
    for caption, g in aps.groupby(0):
        left = offset_cols+6+g.index[0]
        right = offset_cols+6+g.index[-1]
        if left == right:
             worksheet.write(offset_rows, left, caption, format_action)
        else:
             worksheet.merge_range(offset_rows, left, offset_rows, right, caption, format_action)

    # header
    format_header = workbook.add_format({
        'rotation': 90,
        'bold': True,
        'valign': 'vcenter',
        'align': 'center'})
    worksheet.set_row(offset_rows+1, 80)
    for i, caption in enumerate(['yield (MB)', 'coverage'] + [ap['program'] for ap in ACTION_PROGRAMS]):
        worksheet.write(offset_rows+1, offset_cols+4+i, caption, format_header)
    format_spike_seqdate = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'font_size': 8})
    worksheet.write(offset_rows+1, offset_cols+6+len(ACTION_PROGRAMS), 'sequenced at', format_spike_seqdate)

    # header for gene panels
    format_header_genes = workbook.add_format({
        'rotation': 90,
        'bold': False,
        'valign': 'vcenter',
        'align': 'center',
        'font_size': 8})
    if len(gene_order) > 0:
        for caption, g in pd.DataFrame(gene_order).groupby(0):
            left = offset_cols+6+len(ACTION_PROGRAMS)+1+g.index[0]
            right = offset_cols+6+len(ACTION_PROGRAMS)+1+g.index[-1]
            if left == right:
                worksheet.write(offset_rows, left, caption, format_action)
            else:
                worksheet.merge_range(offset_rows, left, offset_rows, right, caption, format_action)
        for i, (panel, gene) in enumerate(gene_order):
            worksheet.write(offset_rows+1, offset_cols+6+len(ACTION_PROGRAMS)+1+i, gene, format_header_genes)
        worksheet.set_column(offset_cols+6+len(ACTION_PROGRAMS)+1, offset_cols+6+len(ACTION_PROGRAMS)+1+len(gene_order), 3)

    worksheet.freeze_panes(offset_rows+2, offset_cols+4)

    # body
    format_project = workbook.add_format({
        'rotation': 90,
        'bold': True,
        'valign': 'vcenter',
        'align': 'center'})
    format_spike_entity_id = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center'})
    format_spike_sampleID = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center'})
    format_spike_entity_role_missing = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center',
        'font_color': '#ff0000'})
    format_gene_coverage_good = workbook.add_format({
        'valign': 'vcenter',
        'align': 'right',
        'font_size': 6,
        'bg_color': '#ccffcc'})
    format_gene_coverage_bad = workbook.add_format({
        'valign': 'vcenter',
        'align': 'right',
        'font_size': 6,
        'bg_color': 'ffcccc'})
    row = offset_rows+2
    for sample_project, grp_project in samplesheets.groupby('Sample_Project'):
        # add in lines to indicate missing samples, e.g. for trios that are incomplete
        missing_samples = []
        for spike_entity_id, grp_spike_entity_group in grp_project.groupby('spike_entity_id'):
            if len(set(grp_spike_entity_group['spike_entity_role'].unique()) & set(['patient', 'father', 'mother', 'sibling'])) > 0:
                for role in ['patient', 'mother', 'father']:
                    if grp_spike_entity_group[grp_spike_entity_group['spike_entity_role'] == role].shape[0] <= 0:
                        missing_samples.append({
                            'spike_entity_id': spike_entity_id,
                            'Sample_ID': role,
                            'spike_entity_role': role,
                            'missing': True,
                        })

        # combine samples from samplesheets AND those that are expected but missing
        samples_and_missing = pd.concat([grp_project, pd.DataFrame(missing_samples)], sort=False).fillna(value={'spike_entity_id': ''})

        worksheet.merge_range(row, offset_cols, row+len(samples_and_missing.groupby(['spike_entity_id', 'Sample_ID']))-1, offset_cols, sample_project.replace('_', '\n'), format_project)
        worksheet.set_column(offset_cols, offset_cols, 4)

        # groupby excludes NaNs, thus I have to hack: replace NaN by "" here and
        # reset to np.nan within the loop
        for spike_entity_group, grp_spike_entity_group in samples_and_missing.groupby('spike_entity_id'):
            if spike_entity_group != "":
                label = spike_entity_group
                if _isKnownDuo(sample_project, spike_entity_group, config):
                    label += "\n(known duo)"
                if len(grp_spike_entity_group.groupby('Sample_ID')) > 1:
                    worksheet.merge_range(row, offset_cols+1, row+len(grp_spike_entity_group.groupby('Sample_ID'))-1, offset_cols+1, label, format_spike_entity_id)
                else:
                    worksheet.write(row, offset_cols+1, label, format_spike_entity_id)
            else:
                spike_entity_group = np.nan
            worksheet.set_column(offset_cols+1, offset_cols+1, 10)

            for nr_sample_id, (sample_id, grp_sample_id) in enumerate(grp_spike_entity_group.sort_values(by='spike_entity_role').groupby('Sample_ID')):
                worksheet.set_column(offset_cols+2, offset_cols+2, 4)
                role = grp_sample_id['spike_entity_role'].iloc[0]
                is_missing = ('missing' in grp_sample_id.columns) and (grp_sample_id[grp_sample_id['missing'] == np.True_].shape[0] > 0)

                # sample_ID, extend field if no spike_entity_group or spike_entity_role is given
                col_start = offset_cols+2
                col_end = offset_cols+2
                if pd.isnull(spike_entity_group):
                    col_start -= 1
                if pd.isnull(role):
                    col_end += 1
                sample_id_value = sample_id
                # if sample_id starts with name of the entity group, we are using "..." to make it visually more pleasing
                if pd.notnull(spike_entity_group) and sample_id_value.startswith(spike_entity_group):
                    sample_id_value = '%s' % sample_id[len(spike_entity_group):]
                frmt = format_spike_sampleID
                if is_missing:
                    sample_id_value = '?'
                    if not _isKnownDuo(sample_project, spike_entity_group, config):
                        frmt = format_spike_entity_role_missing
                if col_start != col_end:
                    worksheet.merge_range(row, col_start, row, col_end, sample_id_value, frmt)
                else:
                    worksheet.write(row, col_start, sample_id_value, frmt)

                # spike_entity_role
                if pd.notnull(role):
                    worksheet.write(row, offset_cols+3, str(role), frmt)

                if is_missing:
                    fmt = format_spike_entity_role_missing
                    if _isKnownDuo(sample_project, spike_entity_group, config):
                        fmt = format_spike_sampleID
                    worksheet.merge_range(row, offset_cols+4, row, offset_cols+3+2+len(ACTION_PROGRAMS)+1, "missing sample", fmt)
                else:
                    # demultiplexing yield
                    frmt = format_bad
                    value_yield = "missing"
                    if (sample_project, sample_id) in data_yields.index:
                        value_yield = float('%.1f' % (int(data_yields.loc[sample_project, sample_id]) / (1000**3)))
                        if value_yield >= 5.0:
                            frmt = format_good
                    if grp_sample_id['Lane'].dropna().shape[0] <= 0:
                        value_yield = 'per sample fastq'
                        frmt = format_good
                    worksheet.write(row, offset_cols+4, value_yield, frmt)
                    worksheet.set_column(offset_cols+4, offset_cols+4, 4)

                    # coverage
                    if ((sample_project, sample_id) in data_coverage.index) and (pd.notnull(data_coverage.loc[sample_project, sample_id])):
                        frmt = format_bad
                        value_coverage = "missing"
                        if (sample_project, sample_id) in data_coverage.index:
                            value_coverage = data_coverage.loc[sample_project, sample_id]
                            if value_coverage >= get_min_coverage(sample_project, config):
                                frmt = format_good
                        worksheet.write(row, offset_cols+5, value_coverage, frmt)
                        worksheet.set_column(offset_cols+5, offset_cols+5, 4)

                    for i, (action, program) in enumerate([(ap['action'], ap['program']) for ap in ACTION_PROGRAMS]):
                        value_numcalls = ""
                        frmt = None
                        if (sample_project, sample_id, action, program) in data_calls:
                            value_numcalls = data_calls.loc[sample_project, sample_id, action, program]
                            if not isinstance(value_numcalls, np.int64):
                                value_numcalls = value_numcalls.iloc[0]
                        if (sample_project, sample_id, action, program) in data_snupy.index:
                            if data_snupy.loc[sample_project, sample_id, action, program]['status']:
                                frmt = format_good
                            else:
                                frmt = format_bad
                        if value_numcalls == RESULT_NOT_PRESENT:
                            value_numcalls = 'vcf missing'
                            frmt = format_bad
                        elif value_numcalls != "" and frmt is None:
                            frmt = format_bad
                        if frmt is not None:
                            worksheet.write(row, offset_cols+6+i, value_numcalls, frmt)

                    # sequencing date
                    worksheet.write(row, offset_cols+6+len(ACTION_PROGRAMS), ' / '.join(sorted(map(
                        lambda x: datetime.datetime.strptime('20%s' % x.split('_')[0], '%Y%m%d').strftime("%Y-%m-%d"), grp_sample_id['run'].unique()))), format_spike_seqdate)
                    worksheet.set_column(offset_cols+6+len(ACTION_PROGRAMS), offset_cols+6+len(ACTION_PROGRAMS), 16)

                    # gene panel coverage
                    if pd.notnull(role):
                        for gene_index, (panel, gene) in enumerate(gene_order):
                            if (sample_project, sample_id, panel, gene) in data_genepanels.index:
                                cov = data_genepanels.loc[sample_project, sample_id, panel, gene]
                                #cov_text = '%i | %.1f | %i' % (cov['mincov'], cov['avgcov_0'], cov['maxcov'])
                                cov_text = '%.1f' % cov['avgcov_0']
                                frmt = format_gene_coverage_bad
                                if cov['avgcov_0'] >= get_min_coverage(sample_project, config):
                                    frmt = format_gene_coverage_good
                                worksheet.write(row, offset_cols+6+len(ACTION_PROGRAMS)+1+gene_index, cov_text, frmt)

                row += 1
    print("done.\n", file=verbose, end="")
    workbook.close()


def _divide_non_zero(numerator, denumerator):
    if denumerator <= 0:
        return 0
    return numerator /  denumerator
def collect_yield_data(dir_flowcell, verbose=None):
    """Composes yield report for (potentially) split demulitplexing.

    Notes
    -----
    People used different length barcodes in the same lane / flowcell.
    This should in general be avoided, due to potential clashes of barcodes.
    Therefore Illumina's bcl2fastq fails to process sample sheets formatted
    like this. I overcome this issue by splitting the samplesheet and running
    bcl2fast multiple times independently. However, this requires complicated
    logic especially for stats for the undetermined reads.

    Parameters
    ----------
    dir_flowcell : str
        Filepath to demultiplexing results in split fashion.

    Returns
    -------
    3-tuple of pd.DataFrame : (Flowcell Summary, Lane Summary, Top Unknown Barcodes)
    Formatting is required before data will resemble Illumina's original yield report.
    """
    lane_summary = []
    lane_meta = []
    cluster_meta = []
    unknown_barcodes = []
    if verbose:
        verbose.write('collect_yield_data(%s):\n' % dir_flowcell)
    parts = glob(join(dir_flowcell, 'part_*'))
    for num_part, dir_part in enumerate(parts):
        if verbose:
            verbose.write('  part %i of %i\n' % (num_part+1, len(parts)))
        clusters = []
        for fp_fastqsummary in glob(join(dir_part, 'Stats/FastqSummaryF*L*.txt')):
            cluster = pd.read_csv(fp_fastqsummary, sep="\t")
            perc_pf_clusters = pd.concat([cluster.groupby(['SampleNumber'])['NumberOfReadsPF'].sum(),
                                          cluster.groupby(['SampleNumber'])['NumberOfReadsRaw'].sum()], axis=1)
            perc_pf_clusters['Lane'] = fp_fastqsummary.split('/')[-1].split('.')[0].split('L')[-1]
            clusters.append(perc_pf_clusters)
        clusters = pd.concat(clusters).reset_index().set_index(['Lane', 'SampleNumber'])
        cluster_meta.append(clusters)

        meta_samples = pd.read_csv(join(dir_part, 'Stats/AdapterTrimming.txt'), sep="\t", usecols=[0,2,3,4])
        meta_samples = meta_samples.iloc[:meta_samples[meta_samples['Lane'].apply(lambda x: x.startswith('Lane:'))].index.max()-2,:].drop_duplicates()
        meta_samples.set_index(['Lane', 'Sample Id'], inplace=True)

        fp_json = join(dir_part, 'Stats/Stats.json')
        part_stats = json.load(open(fp_json, 'r'))
        # sample numbers in FastqSummary files match S-idx numbers, which are only increased if sample is not seen before, independent on lane
        sample_numbers = dict()
        for res_conv in part_stats['ConversionResults']:
            numq30bases = 0
            sumQualityScore = 0
            for res_demux in res_conv['DemuxResults']:
                q30bases = sum([res_metrics['YieldQ30'] for res_metrics in res_demux['ReadMetrics']])
                qualityScore = sum([res_metrics['QualityScoreSum'] for res_metrics in res_demux['ReadMetrics']])
                if res_demux['SampleId'] not in sample_numbers:
                    sample_numbers[res_demux['SampleId']] = len(sample_numbers)+1
                sample_number = sample_numbers[res_demux['SampleId']]
                sample_result = {
                    'Lane': res_conv['LaneNumber'],
                    'Project': meta_samples.loc[str(res_conv['LaneNumber']), res_demux['SampleId']]['Project'],
                    'Sample': res_demux['SampleId'],
                    'PF Clusters': res_demux['NumberReads'],
                    '% of the lane': _divide_non_zero(res_demux['NumberReads'], res_conv['TotalClustersPF']),
                    'Yield': res_demux['Yield'],
                    '% PF Clusters': _divide_non_zero(clusters.loc[str(res_conv['LaneNumber']), sample_number]['NumberOfReadsPF'], clusters.loc[str(res_conv['LaneNumber']), sample_number]['NumberOfReadsRaw']),
                    '% >= Q30 bases': _divide_non_zero(q30bases, res_demux['Yield']),
                    'Q30 bases': q30bases,
                    'QualityScoreSum': qualityScore,
                    'Mean Quality Score': _divide_non_zero(sum([res_metrics['QualityScoreSum'] for res_metrics in res_demux['ReadMetrics']]), res_demux['Yield']),
                    'Sample_Number': sample_number,
                    # default values
                    'Barcode sequence': 'unknown',
                    '% Perfect barcode': 1,
                    '% One mismatch barcode': np.nan,
                }
                if 'IndexMetrics' in res_demux:
                    sample_result['Barcode sequence'] = res_demux['IndexMetrics'][0]['IndexSequence']
                    sample_result['Barcode length'] = len(res_demux['IndexMetrics'][0]['IndexSequence'])
                    sample_result['% Perfect barcode'] = _divide_non_zero(sum([res_idx['MismatchCounts']['0'] for res_idx in res_demux['IndexMetrics']]), res_demux['NumberReads'])
                    sample_result['% One mismatch barcode'] = _divide_non_zero(sum([res_idx['MismatchCounts']['1'] for res_idx in res_demux['IndexMetrics']]), res_demux['NumberReads'])
                lane_summary.append(sample_result)
                numq30bases += q30bases
                sumQualityScore += qualityScore
            if 'Undetermined' in res_conv:
                numq30bases += sum([res_metrics['YieldQ30'] for res_metrics in res_conv['Undetermined']['ReadMetrics']])
                sumQualityScore += sum([res_metrics['QualityScoreSum'] for res_metrics in res_conv['Undetermined']['ReadMetrics']])

            lane_meta.append({
                'Lane': res_conv['LaneNumber'],
                "TotalClustersRaw" : res_conv['TotalClustersRaw'],
                "TotalClustersPF" : res_conv['TotalClustersPF'],
                "Yield": res_conv['Yield'],
                "YieldQ30": numq30bases,
                "QualityScoreSum": sumQualityScore,
                "Flowcell": part_stats['Flowcell'],
            })
        for res_unknown in part_stats['UnknownBarcodes']:
            for barcode in res_unknown['Barcodes'].keys():
                res_barcode = {
                    'Lane': res_unknown['Lane'],
                    'Count': res_unknown['Barcodes'][barcode],
                    'Barcode length': 0,
                    'Barcode sequence': 'unknown'}
                if 'Barcodes' in res_unknown:
                    res_barcode['Barcode sequence'] = barcode
                    res_barcode['Barcode length'] = len(barcode)
                unknown_barcodes.append(res_barcode)

    lane_meta = pd.DataFrame(lane_meta).drop_duplicates().set_index('Lane')
    lane_meta['run'] = dir_flowcell.split('/')[-1]
    lane_summary = pd.DataFrame(lane_summary)
    cluster_meta = pd.concat(cluster_meta)

    if len(set(lane_summary['Barcode sequence'].unique()) - set(['unknown'])) > 0:
        undetermined = []
        for lane, clst_lane in lane_summary.groupby('Lane'):
            undetermined.append({
                'Lane': lane,
                'Project': 'default',
                'Sample': 'Undetermined',
                'Barcode sequence': 'unknown',
                'PF Clusters': lane_meta.loc[lane, 'TotalClustersPF'] - clst_lane['PF Clusters'].sum(),
                '% of the lane': _divide_non_zero(lane_meta.loc[lane, 'TotalClustersPF'] - clst_lane['PF Clusters'].sum(), lane_meta.loc[lane, 'TotalClustersPF']),
                '% Perfect barcode': 1,
                '% One mismatch barcode': np.nan,
                'Yield': lane_meta.loc[lane, 'Yield'] - clst_lane['Yield'].sum(),
                '% PF Clusters': _divide_non_zero(lane_meta.loc[lane, 'TotalClustersPF'] - cluster_meta.loc[str(lane), range(1, cluster_meta.index.get_level_values('SampleNumber').max()), :]['NumberOfReadsPF'].sum(), (lane_meta.loc[lane, 'TotalClustersRaw'] - cluster_meta.loc[str(lane), range(1, cluster_meta.index.get_level_values('SampleNumber').max()), :]['NumberOfReadsRaw'].sum())),
                '% >= Q30 bases': _divide_non_zero(lane_meta.loc[lane, 'YieldQ30'] - lane_summary[lane_summary['Lane'] == lane]['Q30 bases'].sum(), lane_meta.loc[lane, 'Yield'] - lane_summary[lane_summary['Lane'] == lane]['Yield'].sum()),
                'Mean Quality Score': _divide_non_zero(lane_meta.loc[lane, 'QualityScoreSum'] - lane_summary[lane_summary['Lane'] == lane]['QualityScoreSum'].sum(), lane_meta.loc[lane, 'Yield'] - lane_summary[lane_summary['Lane'] == lane]['Yield'].sum()),
            })
        lane_summary = pd.concat([lane_summary, pd.DataFrame(undetermined)], sort=False)

    # traverse unknown barcodes and filter those that are actually used by samples
    # this is non-trivial, because due to splitting, used barcodes might have different sizes!
    unknown_barcodes = pd.DataFrame(unknown_barcodes)
    if len(set(unknown_barcodes['Barcode sequence'].unique()) - set(['unknown'])) > 0:
        idx_remove = []
        for lane in unknown_barcodes['Lane'].astype(int).unique():
            lane_known_bcs = lane_summary[(lane_summary['Lane'] == lane) & (lane_summary['Barcode sequence'] != 'unknown')][['Barcode sequence', 'Barcode length']]
            lane_unknown_bcs = unknown_barcodes[unknown_barcodes['Lane'] == lane]
            remove = set()
            for len_known, known in lane_known_bcs.groupby('Barcode length'):
                for len_unknown, unknown in lane_unknown_bcs.groupby('Barcode length'):
                    if len_known == len_unknown:
                        remove |= set(unknown['Barcode sequence']) & set(known['Barcode sequence'])
                    elif len_known < len_unknown:
                        for _, bc_unknown in unknown['Barcode sequence'].iteritems():
                            if bc_unknown[:int(len_known)] in set(known['Barcode sequence']):
                                remove |= set([bc_unknown])
                    elif len_known > len_unknown:
                        remove |= set(unknown['Barcode sequence']) & set(known['Barcode sequence'].apply(lambda x: x[:int(len_unknown)]))
            idx_remove.extend(lane_unknown_bcs[lane_unknown_bcs['Barcode sequence'].isin(remove)].index)
        top_unknown_barcodes = unknown_barcodes.loc[set(unknown_barcodes.index) - set(idx_remove),:]
    else:
        top_unknown_barcodes = unknown_barcodes

    return lane_meta, lane_summary, top_unknown_barcodes


def create_html_yield_report(fp_yield_report, lane_meta, lane_summary, top_unknown_barcodes, config):
    """Creates a HTML yield report.

    Parameters
    ----------
    fp_yield_report : str
        Filepath to resulting HTML file.
    lane_meta : pd.DataFrame
        Result of collect_yield_data, first component:
        Flowcell summary statistics, basically # clusters.
    lane_summary : pd.DataFrame
        Result of collect_yield_data, second component:
        Sample demultiplexing information.
    top_unknown_barcodes : pd.DataFrame
        Result of collect_yield_data, third component:
        Infos about high abundant unused barcodes.
    config : dict
        Snakemakes configuration dictionary.
    """
    #dir_flowcell : str
    #    Filepath to demultiplexing results in split fashion.

    out = '<html>\n<head>\n'
    #out += '<link rel="stylesheet" href="Report.css" type="text/css">\n'
    out += '<style>\n'
    out += 'body {font-size: 100%; font-family:monospace;}\n'
    out += 'table#ReportTable {border-width: 1px 1px 1px 1px; border-collapse: collapse;}\n'
    out += 'table#ReportTable td {text-align:right; padding: 0.3em;}\n'
    out += 'thead { display: table-header-group }\n'
    out += 'tfoot { display: table-row-group }\n'
    out += 'tr { page-break-inside: avoid }\n'
    out += '</style>\n'
    out += '<title>%s</title>\n' % lane_meta['run'].unique()[0]
    out += '</head>\n<body>'
    out += "%s / [all projects] / [all samples] / [all barcodes]" % lane_meta['Flowcell'].unique()[0]


    out += "<h2>Flowcell Summary</h2>"
    fc_summary = lane_meta.sum()[['TotalClustersRaw', 'TotalClustersPF', 'Yield']].to_frame().T
    fc_summary.rename(columns={'TotalClustersRaw': 'Clusters (Raw)',
                               'TotalClustersPF': 'Clusters(PF)',
                               'Yield': 'Yield (MBases)'}, inplace=True)
    for col in ['Clusters (Raw)', 'Clusters(PF)']:
        fc_summary[col] = fc_summary[col].apply(lambda x: '{:,}'.format(x))
    fc_summary['Yield (MBases)'] = fc_summary['Yield (MBases)'].apply(lambda x: '{:,}'.format(int(x/1000000)))
    out += fc_summary.to_html(index=False, table_id='ReportTable', justify='center')


    out += "<h2>Lane Summary</h2>"
    x = lane_summary.sort_values(['Lane', 'Sample_Number'])[['Lane', 'Project', 'Sample', 'Barcode sequence', 'PF Clusters', '% of the lane', '% Perfect barcode', '% One mismatch barcode', 'Yield', '% PF Clusters', '% >= Q30 bases', 'Mean Quality Score']].rename(columns={'Yield': 'Yield (Mbases)'})
    x['PF Clusters'] = x['PF Clusters'].apply(lambda x: '{:,}'.format(x))
    for col in ['% of the lane', '% Perfect barcode', '% PF Clusters', '% >= Q30 bases', '% One mismatch barcode']:
        x[col] = x[col].apply(lambda x: '' if x == 0 else '%.2f' % (x*100))
    x['Yield (Mbases)'] = x['Yield (Mbases)'].apply(lambda x: '{:,}'.format(int(x/1000000)))
    x['Mean Quality Score'] = x['Mean Quality Score'].apply(lambda x: '' if x == 0 else '%.2f' % x)
    x['% One mismatch barcode'] = x['% One mismatch barcode'].replace('nan', 'NaN')
    out += x.to_html(index=False, table_id='ReportTable', justify='center')


    out += "<h2>Top Unknown Barcodes</h2>"
    topX = 10
    lanes = []
    for lane, lane_barcodes in top_unknown_barcodes.groupby('Lane'):
        x = lane_barcodes.sort_values('Count', ascending=False).iloc[:topX][['Lane', 'Count', 'Barcode sequence']].rename(columns={'Barcode sequence': 'Sequence'})#.reset_index().set_index(['Lane', 'Count'])
        x.index = range(1,topX+1)[:x.shape[0]]
        x['Count'] = x['Count'].apply(lambda x: '{:,}'.format(x))
        lanes.append(x)
    topunknown = pd.concat(lanes, axis=1)
    out += '<table border="1" class="dataframe" id="ReportTable">\n<thead>\n<tr style="text-align: center;">\n'
    for col in topunknown.columns:
        out += '<th>%s</th>\n' % col
    out += '</tr>\n</thead>\n<tbody>\n'
    for i, row in topunknown.iterrows():
        out += '<tr>\n'
        for col, value in row.iteritems():
            if col == 'Lane':
                if i == 1:
                    out += '<th rowspan=%s>%s</th>\n' % (min(topX, topunknown.shape[0]), value)
            else:
                out += '<td>%s</td>\n' % value
        out += '</tr>\n'
    out += '</tbody>\n</table>\n'


    out += "Report generated by %s.</body>\n</html>" % config['name_program']

    with open(fp_yield_report, 'w') as f:
        f.write(out)


def _agilent_annotation_to_genenames(annotation, field):
    """Splits Agilent capture kit annotations into key-value pairs for different database sources.

    Parameters
    ----------
    annotation : str
        Annotation line of Agilent coverage.bed file, column 4.
    field : str
        Name of database entry to be returned.

    Returns
    -------
    str : Name of entry.

    Notes
    -----
    Should a reference database provide multiple names, only the first is used!"""

    gene_names = dict()
    if ',' not in annotation:
        return np.nan
    for entry in annotation.split(','):
        db_name, _id = entry.split('|')
        if db_name not in gene_names:
            gene_names[db_name] = _id
    return gene_names.get(field, np.nan)


def get_gene_panel_coverage(fp_genepanel, fp_bamstat, fp_agilent_coverage, fp_output):
    """Looks up gene coverage for given panel in given sample, based on bamstat.

    Parameters
    ----------
    fp_genepanel : str
        Filepath to yaml gene panel configuration file.
    fp_bamstat : str
        Filepath to bamstat output.
    fp_agilent_coverage: str
        Filepath to original Agilent coverage bed file, i.e. with gene names.
    fp_output : str
        Filepath for output filename.
    """
    # load gene panel definition
    if not exists(fp_genepanel):
        raise ValueError("Gene panel file '%s' does not exist." % fp_genepanel)
    panel = yaml.load(open(fp_genepanel, 'r'))

    # read capture kit probe positions, including gene names
    probes = pd.read_csv(fp_agilent_coverage, sep="\t", header=None, skiprows=2)
    probes[3] = probes[3].apply(lambda x: _agilent_annotation_to_genenames(x, panel['reference_name']))
    probes.columns = ['chromosome', 'start', 'end', 'gene']

    # subset probes to those covering genes of the panel
    probes_of_interest = probes[probes['gene'].isin(panel['genes'])]

    # load coverage information
    coverage = pd.read_csv(fp_bamstat, sep="\s+", dtype=str).rename(columns={'#chrom': 'chromosome'})
    coverage['chromosome'] = coverage['chromosome'].apply(lambda x: 'chr%s' % x)
    for field in ['start', 'end', 'mincov', 'maxcov']:
        coverage[field] = coverage[field].astype(int)
    coverage['avgcov_0'] = coverage['avgcov_0'].astype(float)

    # subset coverage for genes of interest
    # we might experience +/- 1 errors. Currently only observed for -1 if "start"
    for i in range(2):
        coverage_per_probe = probes_of_interest.merge(coverage, left_on=['chromosome', 'start', 'end'], right_on=['chromosome', 'start', 'end'])
        if coverage_per_probe.shape[0] == 0:
            coverage['start'] -= 1
        else:
            break

    # determine min, max, avg coverage across all probes of the chip per gene
    result = pd.concat([
        coverage_per_probe.groupby('gene')['mincov'].min(),
        coverage_per_probe.groupby('gene')['avgcov_0'].mean(),
        coverage_per_probe.groupby('gene')['maxcov'].max()], axis=1)

    result.to_csv(fp_output, sep="\t", index=True, index_label='gene')
