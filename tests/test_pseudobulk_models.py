import csv
import json
import os
from pathlib import Path
import random
import tempfile
import unittest

from extremarank.models import run_refits
from extremarank.pseudobulk import aggregate_pseudobulk
from extremarank.preparation import read_table


def write(path, rows):
    with path.open('w',newline='') as f:csv.writer(f,delimiter='\t').writerows(rows)


class PseudobulkTests(unittest.TestCase):
    def setUp(self):self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
    def tearDown(self):self.temp.cleanup()
    def metadata(self):
        path=self.root/'cells.tsv'
        write(path,[['cell_id','sample_id','donor_id','group','cell_type','batch'],
                    ['c1','s1','d1','T','B','a'],['c2','s2','d1','T','B','a'],
                    ['c3','s3','d2','R','B','b'],['c4','s4','d2','R','Mono','b']])
        return path
    def test_exact_sums_and_technical_library_aggregation(self):
        matrix=self.root/'matrix.tsv';write(matrix,[['feature_id','c1','c2','c3','c4'],['g1',1,2,3,4],['g2',0,5,0,7]])
        report=aggregate_pseudobulk(matrix,self.metadata(),self.root/'out',min_cells=1,merge_technical=True)
        self.assertEqual(report['n_cells'],4);self.assertEqual(report['n_pseudobulks'],3)
        h,r,_=read_table(self.root/'out/celltype_1/matrix.tsv');self.assertEqual(r,[['g1','3','3'],['g2','5','0']])
        h,r,_=read_table(self.root/'out/celltype_1/metadata.tsv');meta=[dict(zip(h,x)) for x in r]
        self.assertEqual(meta[0]['n_cells'],'2');self.assertEqual(json.loads(meta[0]['source_samples']),['s1','s2'])
    def test_10x_duplicate_coordinates_are_additive_and_non_gene_excluded(self):
        folder=self.root/'10x';folder.mkdir();(folder/'barcodes.tsv').write_text('c1\nc2\nc3\nc4\n')
        (folder/'features.tsv').write_text('g1\tG1\tGene Expression\nab\tAB\tAntibody Capture\n')
        (folder/'matrix.mtx').write_text('%%MatrixMarket matrix coordinate integer general\n% comment\n2 4 4\n1 1 1\n1 1 2\n1 2 4\n2 3 100\n')
        report=aggregate_pseudobulk(folder,self.metadata(),self.root/'out','10x',1,merge_technical=True)
        self.assertEqual(report['excluded_non_gene_features'],1)
        self.assertEqual(read_table(self.root/'out/celltype_1/matrix.tsv')[1],[['g1','7','0']])
    def test_missing_cells_fractional_counts_and_inconsistent_metadata_rejected(self):
        matrix=self.root/'matrix.tsv';write(matrix,[['feature_id','c1'],['g1',1]])
        with self.assertRaises(ValueError):aggregate_pseudobulk(matrix,self.metadata(),self.root/'out',min_cells=1)
        write(matrix,[['feature_id','c1','c2','c3','c4'],['g1',1,2.5,3,4]])
        with self.assertRaises(ValueError):aggregate_pseudobulk(matrix,self.metadata(),self.root/'out',min_cells=1)
    def test_low_cell_exclusions_are_explicit(self):
        matrix=self.root/'matrix.tsv';write(matrix,[['feature_id','c1','c2','c3','c4'],['g1',1,2,3,4]])
        report=aggregate_pseudobulk(matrix,self.metadata(),self.root/'out',min_cells=2,merge_technical=True)
        self.assertEqual(report['n_pseudobulks'],1);self.assertEqual(len(report['excluded_low_cell_units']),2)

    def test_multiple_biological_samples_are_not_merged_implicitly(self):
        matrix=self.root/'matrix.tsv';write(matrix,[['feature_id','c1','c2','c3','c4'],['g1',1,2,3,4]])
        with self.assertRaisesRegex(ValueError,'merge-technical'):
            aggregate_pseudobulk(matrix,self.metadata(),self.root/'out',min_cells=1)


class ModelInputTests(unittest.TestCase):
    def test_tpm_and_repeated_donors_rejected_before_fitting(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);matrix=root/'x.tsv';meta=root/'m.tsv'
            write(matrix,[['feature_id','a','b','c','d','e','f'],['g',1,2,3.5,4,5,6]])
            write(meta,[['sample_id','group','donor_id'],*[[s,'T' if i<3 else 'R',s] for i,s in enumerate('abcdef')]])
            with self.assertRaises(ValueError):run_refits(matrix,meta,'T','R','DESeq2',root/'out',k=1)
            write(meta,[['sample_id','group','donor_id'],*[[s,'T' if i<3 else 'R','same'] for i,s in enumerate('abcdef')]])
            with self.assertRaises(ValueError):run_refits(matrix,meta,'T','R','limma',root/'out',k=1)


@unittest.skipUnless(os.environ.get('EXTREMARANK_R_TESTS')=='1','optional installed R/Bioconductor integration')
class RModelTests(unittest.TestCase):
    def test_all_four_models_and_R_topk_export(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);rng=random.Random(77331);matrix=root/'x.tsv';meta=root/'m.tsv';ids=tuple('abcdefgh')
            write(matrix,[['feature_id',*ids],*[[f'g{g:03}',*[rng.randint(20,180)+(80 if g<10 and i<4 else 0) for i in range(8)]] for g in range(100)]])
            write(meta,[['sample_id','group','donor_id','batch'],*[[s,'T' if i<4 else 'R',s,'A' if i%2 else 'B'] for i,s in enumerate(ids)]])
            for model in ('limma','limma-voom','edgeR','DESeq2'):
                out=root/model;report=run_refits(matrix,meta,'T','R',model,out,k=5,categorical=['batch'])
                self.assertIn(report['status'],('OBSERVED_STABLE','OBSERVED_CHANGED'));self.assertEqual(report['valid_refits'],8)
                h,r,_=read_table(out/'expected_topk.tsv');expected=[dict(zip(h,row)) for row in r]
                self.assertEqual([row['feature_id'] for row in expected if row['scenario_id']=='baseline' and row['direction']=='up'],report['baseline_topk'])
    def test_confounded_covariate_keeps_non_evaluable_state(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);matrix=root/'x.tsv';meta=root/'m.tsv'
            write(matrix,[['feature_id',*'abcdef'],['g1',1,2,3,4,5,6],['g2',8,7,6,1,2,3]])
            write(meta,[['sample_id','group','batch'],*[[s,'T' if i<3 else 'R','T' if i<3 else 'R'] for i,s in enumerate('abcdef')]])
            report=run_refits(matrix,meta,'T','R','limma',root/'out',k=1,categorical=['batch'])
            self.assertEqual(report['status'],'NOT_EVALUABLE');self.assertIn('rank deficient',report['fit_status'][0]['detail'])


if __name__=='__main__':unittest.main()
