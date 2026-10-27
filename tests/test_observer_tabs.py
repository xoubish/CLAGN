"""The tab CSVs preserve distinct planned, observed and October target sets."""
import csv
import io
import json
import re
import shutil
import subprocess
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


class ObserverTabTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        page=(ROOT/'docs/index.html').read_text()
        cls.data=json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>',page,re.S).group(1))

    def displayed(self,name):
        return self.data['display_names'].get(name,name)

    def rows(self,tab):
        d=self.data
        return list(csv.DictReader(io.StringIO(d['files'][d['tab_downloads'][tab]])))

    @unittest.skipUnless(shutil.which('node'), 'Node is needed to execute page navigation')
    def test_about_navigation_renders_details_on_demand(self):
        result=subprocess.run(['node',str(ROOT/'tests/observer_navigation.cjs')],
            cwd=ROOT,capture_output=True,text=True,timeout=30)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_one_distinct_export_for_each_target_tab(self):
        downloads=self.data['tab_downloads']
        self.assertEqual(set(downloads),{'all','sep23','observed','oct26','oct27'})
        self.assertEqual(len(set(downloads.values())),5)
        for tab,name in downloads.items():
            self.assertTrue(all(c==c.upper() for c in self.rows(tab)[0]))
            self.assertEqual((ROOT/'docs/targets'/name).read_text(),self.data['files'][name])

    def test_all_contains_both_prepared_pool_and_added_october_targets(self):
        expected={t['name'] for t in self.data.get('field_cards',[])}|{self.displayed(t['name']) for t in self.data['targets']}|{self.displayed(r['name']) for r in self.data['october']['rows']}
        rows=self.rows('all')
        self.assertEqual({r['NAME'] for r in rows},expected)
        self.assertEqual(len(rows),len(expected))
        if self.data.get('field_cards'):
            star=next(r for r in rows if r['NAME']=='Z16CygABC')
            self.assertEqual(star['EXPTIME'],'SET 1')
            self.assertEqual(star['SLITANGLE'],'133.04')
        with (ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv').open() as f:
            protected={r['internal_id'] for r in csv.DictReader(f)}
        self.assertTrue({self.displayed(n) for n in protected}<={r['NAME'] for r in rows})
        for t in self.data['all_targets']:
            self.assertEqual(t['display_name'], ('J'+t['name']) if t['name'] in protected and not t['name'].startswith('J') else t['name'])

    def test_public_artifacts_omit_proposal_metadata_without_dropping_targets(self):
        page=(ROOT/'docs/index.html').read_text()
        # Random base64 image bytes are not text or metadata.
        page=re.sub(r'data:[^;]+;base64,[A-Za-z0-9+/=]+','data:omitted',page)
        self.assertIsNone(re.search(r'jwst|miri|proposal',page.lower()),'Proposal metadata leaked into public page')
        self.assertFalse((ROOT/'docs/october/jwst_coverage.csv').exists())
        self.assertNotRegex((ROOT/'docs/README.md').read_text().lower(),r'jwst|miri|proposal')
        for directory in ['october','targets']:
            for p in (ROOT/'docs'/directory).glob('*.csv'):
                self.assertNotRegex(p.read_text().lower(),r'jwst|miri|proposal',str(p))
        for night in ['oct26','oct27']:
            source=ROOT/f'observing/oct26_27_2026/{night}_ngps.csv'
            if not source.exists():continue
            with source.open() as f:
                local=list(csv.DictReader(f))
            public=self.rows(night)
            for r,p in zip(local,public):
                for key in ['NAME','RA','DECL','EXPTIME','NEXP','SLITWIDTH','BINSPAT','BINSPECT','SLITANGLE','AIRMASS_MAX']:
                    self.assertEqual(self.displayed(r[key]) if key=='NAME' else r[key],p[key],(night,r['NAME'],key))
            self.assertEqual(len(public),len(local))

    def test_planned_and_observed_preserve_their_own_exposures(self):
        planned=self.rows('sep23');observed=self.rows('observed')
        self.assertEqual(len(planned),14)
        self.assertEqual(len(observed),18)
        self.assertEqual([r['NAME'] for r in planned],[self.displayed(n) for n,s in sorted(self.data['sep23_sequence'].items(),key=lambda v:v[1]['rank'])])
        self.assertEqual([r['NAME'] for r in observed],[self.displayed(r['name']) for r in self.data['observed']])
        self.assertEqual(sum(int(r['NEXP']) for r in observed),39)
        for r,o in zip(observed,self.data['observed']):
            self.assertEqual(float(r['TOTAL_EXPTIME_S']),o['exposure_seconds'])
        for night in ['oct26','oct27']:
            rows=self.rows(night)
            self.assertEqual(len(rows),43 if night=='oct26' else 42)
            rows=rows[1:-1]
            rows=[r for r in rows if r['NAME']!='Z16CygABC']
            self.assertTrue(all(r['NEXP']==({'P8544':'3','HIRAS-F23290-1252':'4'}.get(r['NAME'],'2')) and r['EXPTIME']=='SET 300' for r in rows))

    def test_october_uploads_include_standards_at_both_ends(self):
        for night in ['oct26','oct27']:
            rows=self.rows(night)
            standards=list(csv.DictReader(io.StringIO(self.data['files'][night+'_standards_ngps.csv'])))
            self.assertEqual(rows[0],standards[0])
            self.assertEqual(rows[-1],standards[-1])
            self.assertEqual(rows[0]['NAME'],'BD+28 4211')
            self.assertEqual(rows[-1]['NAME'],'Feige 34')
            self.assertEqual([r['EXPTIME'] for r in [rows[0],rows[-1]]],['SET 10','SET 20'])
            for standard in [rows[0],rows[-1]]:
                self.assertEqual((standard['NEXP'],standard['SLITWIDTH'],standard['BINSPAT'],standard['BINSPECT']),('2','SET 1.5','2','3'))
                self.assertEqual(sum(r['NAME']==standard['NAME'] for r in rows),1)
            coordinates=list(csv.reader(io.StringIO(self.data['files'][night+'_coordinates.csv'])))
            self.assertEqual(coordinates,[[r[k] for k in ['NAME','RA','DECL']] for r in rows])

    def test_review_and_public_downloads_match_packet(self):
        october=self.data['october']
        self.assertEqual(october['review_date'],'2026-10-22')
        self.assertAlmostEqual(october['september_performance']['median_overhead_minutes'],4.4083664,places=5)
        rows=october['rows']
        for night in ['oct26','oct27']:
            name=night+'_ngps.csv'
            self.assertEqual((ROOT/'docs/october'/name).read_text(),self.data['files'][name])
            exported=self.rows(night)[1:-1]
            selected=[r for r in rows if r['night']==night]+[dict(r,name=r['upload_name']) for r in october.get('manual_visits',[]) if r['night']==night]
            selected.sort(key=lambda r:r['start_pdt'])
            self.assertEqual([r['NAME'] for r in exported],[self.displayed(r['name']) for r in selected])
            for upload,r in zip(exported,selected):
                self.assertEqual(int(upload['NEXP']),int(r['exposures']))
        protected=next(r for r in rows if r['name']=='P8544')
        self.assertEqual(protected['night'],'oct27')
        self.assertEqual(int(protected['visit_minutes']),20)


if __name__=='__main__':unittest.main()
