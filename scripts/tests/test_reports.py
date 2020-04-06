import sys
sys.path.append("..")

from unittest import TestCase, main
import tempfile
import pickle
import os
from scripts.reports import *
import pandas as pd
import sys

class ChecksTests(TestCase):
    fp_tmpfile = None

    def setUp(self):
        _, self.fp_tmpfile = tempfile.mkstemp()

    def tearDown(self):
        print("STEFAN", self.fp_tmpfile)
        #os.remove(self.fp_tmpfile)

    def test_write_status_update(self):
        with open('scripts/tests/data/statusupdate_keimbahn454.pickle',
                  'rb') as f:
            (data_yields, data_coverage, data_snupy,
             data_calls, ss, config) = pickle.load(f)
        write_status_update(
            (data_yields, data_coverage,
             data_snupy, data_calls, pd.DataFrame()),
            self.fp_tmpfile, ss, config)
        obs = pd.read_excel(self.fp_tmpfile, header=None)
        exp = pd.read_excel(
            'scripts/tests/data/statusupdate_keimbahn454.xlsx', header=None)
        for obj in [exp, obs]:
            # this cell hold timestamp and machine specific information
            obj.iloc[0, 0] = ""
        pd.testing.assert_frame_equal(obs, exp)


if __name__ == '__main__':
    main()
