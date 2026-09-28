"""The tab CSVs preserve distinct planned, observed and October target sets."""
import csv
import io
import json
import re
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


class ObserverTabTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        page=(ROOT/'docs/index.html').read_text()
        cls.data=json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>',page,re.S).group(1))

    def rows(self,tab):
        d=self.data
        return list(csv.DictReader(io.StringIO(d['files'][d['tab_downloads'][tab]])))

    def test_one_distinct_export_for_each_target_tab(self):
        downloads=self.data['tab_downloads']
        self.assertEqual(set(downloads),{'all','sep23','observed','oct26','oct27'})
        self.assertEqual(len(set(downloads.values())),5)
        for tab,name in downloads.items():
            self.assertTrue(all(c==c.upper() for c in self.rows(tab)[0]))
            self.assertEqual((ROOT/'docs/targets'/name).read_text(),self.data['files'][name])

    def test_all_contains_both_prepared_pool_and_added_october_targets(self):
        expected={t['name'] for t in self.data['targets']}|{r['name'] for r in self.data['october']['rows']}
        rows=self.rows('all')
        self.assertEqual({r['NAME'] for r in rows},expected)
        self.assertEqual(len(rows),len(expected))
        self.assertEqual(sum(bool(t['jwst_id']) for t in self.data['all_targets']),24)

    def test_planned_and_observed_preserve_their_own_exposures(self):
        planned=self.rows('sep23');observed=self.rows('observed')
        self.assertEqual(len(planned),14)
        self.assertEqual(len(observed),18)
        self.assertEqual([r['NAME'] for r in planned],[n for n,s in sorted(self.data['sep23_sequence'].items(),key=lambda v:v[1]['rank'])])
        self.assertEqual([r['NAME'] for r in observed],[r['name'] for r in self.data['observed']])
        self.assertEqual(sum(int(r['NEXP']) for r in observed),39)
        for r,o in zip(observed,self.data['observed']):
            self.assertEqual(float(r['TOTAL_EXPTIME_S']),o['exposure_seconds'])
        for night in ['oct26','oct27']:
            rows=self.rows(night)
            self.assertEqual(len(rows),40)
            self.assertTrue(all(r['NEXP']=='2' and r['EXPTIME']=='SET 300' for r in rows))


if __name__=='__main__':unittest.main()
