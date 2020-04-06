import sys
sys.path.append("..")

from unittest import TestCase, main
from scripts.parse_samplesheet import *
from scripts.parse_samplesheet import _run2date
import yaml
from io import StringIO


class ParseSamplesheetTests(TestCase):
    samplesheets = None
    config = None

    def setUp(self):
        self.samplesheets = get_global_samplesheets(
            'scripts/tests/data/', self.config)
        with open('config.yaml', 'r') as f:
            self.config = yaml.load(f, Loader=yaml.SafeLoader)

    def test_get_role(self):
        self.config = {'dirs': {'prefix': './',
                                'inputs': 'tests/',
                                'samplesheets': 'data'}}
        spike_project = 'Alpss'
        with self.assertRaisesRegex(
            ValueError,
            (('Could not find a spike project with name "%s". Available '
              'projects are:\n\t%s\n') %
             (spike_project,
              '\n\t'.join(sorted([
                'AG_Remke', 'Alps', 'Fischer_Geron', 'Maus_Hauer',
                'ALL_family_LB']))))):
            get_role(spike_project, None, None, self.samplesheets)

        spike_project = 'Alps'
        spike_entity_id = 'ALPS_6x'
        with self.assertRaisesRegex(
            ValueError,
            (('Could not find a spike entity group with name "%s". '
              'Available entities for projects "%s" are:\n\t%s\n') %
             (spike_entity_id, spike_project,
              '\n\t'.join(sorted(['ALPS_66']))))):
            get_role(spike_project, spike_entity_id, None, self.samplesheets)

        spike_entity_id = 'ALPS_66'
        spike_entity_role = 'grandma'
        with self.assertRaisesRegex(
            ValueError,
            (('Could not find a role "%s" for spike entity group with '
              'name "%s". Available roles are:\n\t%s\n') %
             (spike_entity_role, spike_entity_id,
              '\n\t'.join(sorted(['father', 'mother', 'patient']))))):
            get_role(spike_project, spike_entity_id,
                     spike_entity_role, self.samplesheets)

        spike_entity_role = 'father'
        exp = 'Alps/ALPS_66_a'
        obs = get_role(spike_project, spike_entity_id,
                       spike_entity_role, self.samplesheets)
        self.assertEqual(exp, obs)

    def test_get_reference_exometrack(self):
        exp = {'Agilent_SureSelect_V5plusUTR.bed': [
            'Keimbahn/KB0198_c',
            'Keimbahn/KB0198_f',
            'Keimbahn/KB0198_m',
            'Keimbahn/KB0005_m',
            'Keimbahn/KB0043_c',
            'Keimbahn/KB0043_m',
            'Keimbahn/KB0077_c',
            'Keimbahn/KB0077_m',
            'Keimbahn/KB0122_c',
            'Keimbahn/KB0122_f',
            'Keimbahn/KB0122_m',
            'Alps/ALPS_64',
            'Alps/ALPS_64_a',
            'Alps/ALPS_64_b',
            'Alps/ALPS_65',
            'Alps/ALPS_65_a',
            'Alps/ALPS_65_b',
            'Keimbahn/KB0135_c',
            'Keimbahn/KB0135_f',
            'Keimbahn/KB0135_m',
            'Keimbahn/KB0136_c',
            'Keimbahn/KB0136_m',
            'Keimbahn/KB0199_c',
            'Keimbahn/KB0199_f',
            'Keimbahn/KB0199_m',
            'Keimbahn/KB0190_c',
            'Keimbahn/KB0190_m',
            'Keimbahn/KB0174_c',
            'Keimbahn/KB0174_m',
            'Keimbahn/KB0206_c',
            'Keimbahn/KB0206_f',
            'Keimbahn/KB0145_c',
            'Keimbahn/KB0145_m',
            'Alps/ALPS_66',
            'Alps/ALPS_66_a',
            'Alps/ALPS_66_b',
            'Fischer_Geron/F49',
            'Fischer_Geron/42'],
        'Agilent_mouseexome_S0276129_Regions_mm10_modified_merged.bed': [
            'Maus_Hauer/mAID_292_N',
            'Maus_Hauer/mAID_292_T',
            'Maus_Hauer/mAID_293_N',
            'Maus_Hauer/mAID_293_T',
            'Maus_Hauer/mMicro295_N',
            'Maus_Hauer/mMicro295_T',
            'Maus_Hauer/mMicro296_N',
            'Maus_Hauer/mMicro296_T',
            'Maus_Hauer/mMicro297_N',
            'Maus_Hauer/mMicro297_T',
            'Maus_Hauer/mMicro298_N',
            'Maus_Hauer/mMicro298_T',
            'Maus_Hauer/286C',
            'Maus_Hauer/286T',
            'Maus_Hauer/287C',
            'Maus_Hauer/287T',
            'Maus_Hauer/288C',
            'Maus_Hauer/288T',
            'Maus_Hauer/291C',
            'Maus_Hauer/291T'],
        "Agilent_SureSelect_Human_AllExon_V6_r2.hGRC37.Covered.bed": [
            'ALL_family_LB/LB_4-2',
            'ALL_family_LB/LB_4-3',
            'ALL_family_LB/LB_5-1',
            'ALL_family_LB/LB_5-2',
            'ALL_family_LB/LB_5-3']}

        for entity in self.samplesheets[~self.samplesheets[
                'Sample_Project'].isin(['AG_Remke'])]['fastq-prefix'].unique():
            obs = get_reference_exometrack(entity, self.samplesheets,
                                           self.config)
            self.assertIn(entity, exp[obs])

    def test_get_bwa_mem_header(self):
        self.assertIn('SureSelectXTV5plusUTRautomated', get_bwa_mem_header(
            'Alps/ALPS_66', self.samplesheets, self.config))
        self.assertIn('SureSelectXTmouse', get_bwa_mem_header(
            'Maus_Hauer/286C', self.samplesheets, self.config))
        self.assertIn('SureSelectHumanAllExonV6r2', get_bwa_mem_header(
            'ALL_family_LB/LB_5-3', self.samplesheets, self.config))

    def test_validate_samplesheet(self):
        # WARNINGS
        data = [{'Lane': None,
                 'Sample_ID': "KB0005_c",
                 'Sample_Name': None,
                 'I7_Index_ID': 'A02',
                 'index': 'AGCAGGAA',
                 'Sample_Project': "Keimbahn",
                 'spike_entity_id': "KB0005",
                 'spike_entity_role': 'patient'}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn(('Sample_ID "KB0005_c" for Sample_Project "Keimbahn" '
                       'in line 22 uses unexpected demultiplexing barcode A02:'
                       ' "AGCAGGAA"'), obs.getvalue())

        data = [{'Lane': None,
                 'Sample_ID': "KB0005_c",
                 'Sample_Name': None,
                 'I7_Index_ID': 'A01',
                 'index': 'AAAACCCC',
                 'Sample_Project': "Keimbahn",
                 'spike_entity_id': "KB0005",
                 'spike_entity_role': 'patient'}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn(('Sample_ID "KB0005_c" for Sample_Project "Keimbahn" in'
                       ' line 22 uses unexpected combination of index and '
                       'I7_index_ID A01: "AAAACCCC"'), obs.getvalue())

        data = [{'Lane': None,
                 'Sample_ID': "KB0005_x",
                 'Sample_Name': None,
                 'I7_Index_ID': 'A01',
                 'index': 'ATGCCTAA',
                 'Sample_Project': "Keimbahn",
                 'spike_entity_id': "KB0005",
                 'spike_entity_role': 'patient'}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn(('Sample_ID "KB0005_x" does not match expected spike_'
                       'entity_role "patient" for Sample_Project "Keimbahn" '
                       'in line 22.'), obs.getvalue())

        data = [{'Lane': None,
                 'Sample_ID': "KB0004_c",
                 'Sample_Name': None,
                 'I7_Index_ID': 'A01',
                 'index': 'ATGCCTAA',
                 'Sample_Project': "Keimbahn",
                 'spike_entity_id': "KB0005",
                 'spike_entity_role': 'patient'}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn(('spike_entity_id "KB0005" is not part of the Sample_ID'
                       ' "KB0004_c" in line 22.'), obs.getvalue())

        data = [{'Lane': None,
                 'Sample_ID': "Kb0004_c",
                 'Sample_Name': None,
                 'I7_Index_ID': 'A01',
                 'index': 'ATGCCTAA',
                 'Sample_Project': "Keimbahn",
                 'spike_entity_id': "KB0005",
                 'spike_entity_role': 'patient'}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn(('Sample_ID "Kb0004_c" does not follow expected naming '
                       'schema "^KB\d{4}" in line 22.'), obs.getvalue())

        data = [{'Lane': None,
                 'Sample_ID': "paul_c",
                 'Sample_Name': None,
                 'I7_Index_ID': None,
                 'index': None,
                 'Sample_Project': "Alps",
                 'spike_entity_id': "paul",
                 'spike_entity_role': 'healthy'}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn(('Sample_ID "paul_c" does not follow expected naming '
                       'schema "^ALPS" in line 22.'), obs.getvalue())

        data = [{'Lane': None,
                 'Sample_ID': "ALPS100",
                 'Sample_Name': None,
                 'I7_Index_ID': None,
                 'index': None,
                 'Sample_Project': "Alpss",
                 'spike_entity_id': "ALPS100",
                 'spike_entity_role': 'patient'}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn('Sample_Project "Alpss" is not described in config.yaml.'
                      ' No processing other than demultiplexing will be'
                      ' applied.', obs.getvalue())

        data = [{'Lane': None,
                 'Sample_ID': "kurt_c",
                 'Sample_Name': None,
                 'I7_Index_ID': None,
                 'index': None,
                 'Sample_Project': "Maus_Hauer",
                 'spike_entity_id': "paul",
                 'spike_entity_role': 'healthy'}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn('spike_entity_id "paul" is not part of the Sample_ID '
                      '"kurt_c" in line 22.', obs.getvalue())

        data = [{'Lane': None,
                 'Sample_ID': "samplea",
                 'Sample_Name': None,
                 'I7_Index_ID': None,
                 'index': None,
                 'Sample_Project': "Maus_Hauer",
                 'spike_entity_id': None,
                 'spike_entity_role': 'brother'}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn('spike_entity_role "brother" in line 22 for Sample_Proj'
                      'ect "Maus_Hauer" is unknown!', obs.getvalue())

        data = [{'Lane': None,
                 'Sample_ID': "samplea",
                 'Sample_Name': None,
                 'I7_Index_ID': None,
                 'index': None,
                 'Sample_Project': "Keimbahn",
                 'spike_entity_id': None,
                 'spike_entity_role': 'brother'}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn('spike_entity_role "brother" in line 22 for Sample'
                      '_Project "Keimbahn" is unknown!', obs.getvalue())

        data = [{'Lane': None,
                 'Sample_ID': None,
                 'Sample_Name': "samplea",
                 'I7_Index_ID': None,
                 'index': None,
                 'Sample_Project': "NoProject",
                 'spike_entity_id': None,
                 'spike_entity_role': None}]
        obs = StringIO()
        validate_samplesheet(pd.DataFrame(data), self.config, err=obs)
        self.assertIn('is not described in config.yaml. No processing other'
                      ' than demultiplexing', obs.getvalue())


        # ERRORS
        with self.assertRaisesRegex(
                ValueError,
                'Samplesheet is missing column\(s\)'):
            data = pd.DataFrame(data=None, columns=['Sample_ID', 'Lane'])
            validate_samplesheet(pd.DataFrame(data), self.config)

        with self.assertRaisesRegex(
                ValueError,
                'contains a restricted character: "Wrong-Chars"'):
            data = [{'Lane': None,
                     'Sample_ID': None,
                     'Sample_Name': "Wrong-Chars",
                     'I7_Index_ID': None,
                     'index': None,
                     'Sample_Project': None,
                     'spike_entity_id': None,
                     'spike_entity_role': None}]
            validate_samplesheet(pd.DataFrame(data), self.config)

        with self.assertRaisesRegex(ValueError, 'has an empty Sample_Project'):
            data = [{'Lane': None,
                     'Sample_ID': None,
                     'Sample_Name': "samplea",
                     'I7_Index_ID': None,
                     'index': None,
                     'Sample_Project': "",
                     'spike_entity_id': None,
                     'spike_entity_role': None}]
            validate_samplesheet(pd.DataFrame(data), self.config)


class HelperFunctions(TestCase):
    def test__run2date(self):
        # Deya: I have tested this function for the following run
        self.assertEqual('2019/10/09', _run2date('191009_SN737_0479_ACDU95ACXX'))
        self.assertEqual('2019/08/21', _run2date('190821_SN737_0472_BCDPEFACXX'))
        # Stefan
        self.assertEqual('2012/03/23', _run2date('120323_SN737_0172_BC0G9NACXX'))
        self.assertEqual('2012/08/01', _run2date('120801_SN737_0199_AC0TNHACXX'))
        self.assertEqual('2014/11/23', _run2date('141123_M02092_0031_000000000-ABGCR'))
        self.assertEqual('2016/06/09', _run2date('160609_SN737_0385_AC8N1RACXX'))
        self.assertEqual('2017/07/25', _run2date('170725_SN737_0416_BCAD2RACXX'))
        self.assertEqual('2018/09/13', _run2date('180913_M04304_0134_000000000-C45VP'))
        self.assertEqual('2018/09/20', _run2date('180920_SN737_0447_ACC9UWACXX'))
        self.assertEqual('2019/08/02', _run2date('190802_SN737_0470_ACDT87ACXX'))
        self.assertEqual('2019/08/02', _run2date('190802_SN737_0471_BCDPE4ACXX'))
        self.assertEqual('2018/09/06', _run2date('180906_M04304_0133_000000000-C2VLP'))
        self.assertEqual('2018/09/20', _run2date('180920_SN737_0448_BCC9V5ACXX'))


if __name__ == '__main__':
    main()
