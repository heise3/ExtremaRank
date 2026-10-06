"""Optional dense/CSR/CSC H5AD counts; explicit layer and no whole-matrix densification."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

from extremarank.pseudobulk import aggregate_pseudobulk
from extremarank.preparation import read_table


@unittest.skipUnless(importlib.util.find_spec('anndata'),'optional AnnData dependency')
class H5ADTests(unittest.TestCase):
    def test_dense_and_sparse_layers_preserve_exact_sums(self):
        import anndata as ad
        import numpy as np
        import pandas as pd
        from scipy import sparse
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);meta=root/'m.tsv'
            meta.write_text('cell_id\tsample_id\tdonor_id\tgroup\tcell_type\n'+'c1\ts1\td1\tT\tB\n'+'c2\ts1\td1\tT\tB\n'+'c3\ts2\td2\tR\tB\n')
            counts=np.array([[1,2],[3,0],[5,6]],dtype=np.int64)
            for kind,raw in [('dense',counts),('csr',sparse.csr_matrix(counts)),('csc',sparse.csc_matrix(counts))]:
                data=ad.AnnData(X=np.log1p(counts),obs=pd.DataFrame(index=['c1','c2','c3']),var=pd.DataFrame(index=['g1','g2']))
                data.layers['counts']=raw;path=root/(kind+'.h5ad');data.write_h5ad(path)
                report=aggregate_pseudobulk(path,meta,root/kind,'h5ad',1,counts_layer='counts')
                self.assertEqual(report['counts_layer'],'counts')
                self.assertEqual(read_table(root/kind/'celltype_1/matrix.tsv')[1],[['g1','4','5'],['g2','2','6']])
                with self.assertRaises(ValueError):aggregate_pseudobulk(path,meta,root/'bad','h5ad',1,counts_layer='X')
                with self.assertRaises(ValueError):aggregate_pseudobulk(path,meta,root/'bad','h5ad',1)
    def test_metadata_cell_metrics_are_ignored_unless_declared_covariates(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);matrix=root/'x.tsv';meta=root/'m.tsv'
            matrix.write_text('feature_id\tc1\tc2\ng1\t1\t2\n')
            meta.write_text('cell_id\tsample_id\tdonor_id\tgroup\tcell_type\tpercent_mito\n'+'c1\ts1\td1\tT\tB\t1\n'+'c2\ts1\td1\tT\tB\t20\n')
            report=aggregate_pseudobulk(matrix,meta,root/'out',min_cells=1)
            self.assertEqual(report['ignored_cell_metadata_columns'],['percent_mito'])
            with self.assertRaises(ValueError):aggregate_pseudobulk(matrix,meta,root/'bad',min_cells=1,covariate_columns=['percent_mito'])


if __name__=='__main__':unittest.main()
