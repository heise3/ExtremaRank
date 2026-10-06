"""Streaming 10X or dense raw-count aggregation by donor/group/cell type."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
from pathlib import Path

from .cli import _write_json, _write_tsv
from .preparation import read_table


def _integer(text):
    try:
        number=Decimal(text)
    except InvalidOperation as e:
        raise ValueError("raw counts must be nonnegative integers") from e
    if not number.is_finite() or number < 0 or number != number.to_integral_value() or number > 2**53-1:
        raise ValueError("raw counts must be nonnegative exact integers <=2^53-1; no normalized/log counts")
    return int(number)


def _open(path):
    return gzip.open(path,"rt",encoding="utf-8-sig",newline="") if path.suffix==".gz" else path.open(encoding="utf-8-sig",newline="")


def _hash(path):
    digest=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def _locate(folder, names):
    for name in names:
        for suffix in ("",".gz"):
            candidate=folder/(name+suffix)
            if candidate.is_file():
                return candidate
    raise ValueError("10X directory lacks "+"/".join(names))


def aggregate_pseudobulk(matrix, metadata, output, input_format="dense", min_cells=10,
                        covariate_columns=(), merge_technical=False, counts_layer=None):
    matrix,metadata,output=Path(matrix),Path(metadata),Path(output)
    if input_format not in {"dense","10x","h5ad"} or isinstance(min_cells,bool) or not isinstance(min_cells,int) or min_cells<1:
        raise ValueError("format must be dense/10x/h5ad; min-cells must be a positive integer")
    h, rows, meta_sha=read_table(metadata)
    required={"cell_id","sample_id","donor_id","group","cell_type"}
    if not required<=set(h):
        raise ValueError("cell metadata requires cell_id/sample_id/donor_id/group/cell_type")
    cells={}; unit_samples={}; samples={}; counts=Counter()
    extra=tuple(covariate_columns)
    if len(set(extra))!=len(extra) or not set(extra)<=set(h)-required:
        raise ValueError("covariates must be distinct nonreserved columns in cell metadata")
    if set(extra)&{"n_cells","source_samples"}:
        raise ValueError("cell metadata uses reserved pseudobulk output fields n_cells/source_samples")
    covariates={}
    for row in rows:
        rec=dict(zip(h,(x.strip() for x in row))); cell=rec["cell_id"]
        if cell in cells or any(not rec[x] for x in required):
            raise ValueError("cell IDs must be unique; required cell metadata cannot be empty")
        sid=rec["sample_id"]; sample=(rec["donor_id"],rec["group"])
        if sid in samples and samples[sid]!=sample:
            raise ValueError("one sample cannot refer to different donors/conditions")
        samples[sid]=sample
        key=(rec["cell_type"],rec["donor_id"],rec["group"])
        values=tuple(rec[x] for x in extra)
        if key in covariates and covariates[key]!=values:
            raise ValueError("covariates must be constant within each donor/condition/cell type")
        covariates[key]=values
        cells[cell]=key; counts[key]+=1; unit_samples.setdefault(key,set()).add(sid)
    if not merge_technical and any(len(samples)>1 for samples in unit_samples.values()):
        raise ValueError("multiple samples share a donor/condition/cell type; declare --merge-technical only for genuine technical libraries, not repeated biological visits")
    keys=tuple(sorted(key for key in counts if counts[key]>=min_cells))
    if not keys:
        raise ValueError("no donor/condition/cell-type units pass min-cells")
    key_index={key:i for i,key in enumerate(keys)}
    retained_cells={cell:key_index[key] for cell,key in cells.items() if key in key_index}
    features=[]; aggregate=[]; sources=[metadata]; excluded_non_gene=0; n_entries=0
    def add(j,cell,value):
        i=retained_cells.get(cell)
        if i is not None:
            aggregate[j][i]+=value
            if aggregate[j][i]>2**53-1:
                raise ValueError("aggregated count exceeds exact downstream binary64 integer range")
    if input_format=="dense":
        mh,mr,msha=read_table(matrix)
        if mh[0]!="feature_id" or set(mh[1:])!=set(cells):
            raise ValueError("dense matrix requires feature_id x cells and exact cell metadata ID sets")
        seen=set()
        for row in mr:
            fid=row[0].strip()
            if not fid or fid in seen:
                raise ValueError("gene IDs must be unique/nonempty")
            seen.add(fid);features.append(fid);aggregate.append([0]*len(keys))
            for cell,value in zip(mh[1:],row[1:]):
                count=_integer(value);n_entries+=count!=0;add(len(features)-1,cell,count)
        sources.append(matrix)
    elif input_format=="h5ad":
        if not counts_layer:
            raise ValueError("H5AD requires an explicit raw-count --counts-layer; use X only when X contains raw counts")
        try:
            import h5py
            from anndata.io import read_elem, sparse_dataset
        except ImportError as e:
            raise ValueError("H5AD input requires the optional extremarank[single-cell] dependencies") from e
        with h5py.File(matrix,'r') as store:
            def names(group):
                index=group.attrs.get('_index')
                if isinstance(index,bytes):index=index.decode('utf-8')
                if not isinstance(index,str) or index not in group:
                    raise ValueError("H5AD obs/var must have an encoded unique index")
                return [str(x) for x in read_elem(group[index])]
            barcodes_list=tuple(names(store['obs']))
            features=names(store['var'])
            if len(set(barcodes_list))!=len(barcodes_list) or set(barcodes_list)!=set(cells):
                raise ValueError("H5AD cell IDs and cell metadata must match exactly")
            if not features or len(set(features))!=len(features) or any(not x for x in features):
                raise ValueError("H5AD gene IDs must be unique/nonempty")
            if counts_layer!='X' and ('layers' not in store or counts_layer not in store['layers']):
                raise ValueError("declared raw-count H5AD layer does not exist; no fallback to X/raw")
            element=store['X'] if counts_layer=='X' else store['layers'][counts_layer]
            data=sparse_dataset(element) if isinstance(element,h5py.Group) else element
            if data.shape!=(len(barcodes_list),len(features)):
                raise ValueError("H5AD count layer dimensions disagree with cell/gene indices")
            aggregate=[[0]*len(keys) for fid in features]
            for start in range(0,len(barcodes_list),256):
                block=data[start:start+256]
                if hasattr(block,'tocoo'):
                    block=block.tocoo()
                    for i,j,value in zip(block.row,block.col,block.data):
                        count=_integer(str(value));n_entries+=count!=0;add(int(j),barcodes_list[start+int(i)],count)
                else:
                    for i,row in enumerate(block):
                        for j,value in enumerate(row):
                            if value!=0:
                                count=_integer(str(value));n_entries+=1;add(j,barcodes_list[start+i],count)
            sources.append(matrix)
    else:
        barcodes=_locate(matrix,("barcodes.tsv",))
        gene_file=_locate(matrix,("features.tsv","genes.tsv"))
        mtx=_locate(matrix,("matrix.mtx",))
        sources.extend((barcodes,gene_file,mtx))
        with _open(barcodes) as f:
            barcodes_list=[line.rstrip("\r\n").split("\t")[0] for line in f]
        if len(set(barcodes_list))!=len(barcodes_list) or set(barcodes_list)!=set(cells):
            raise ValueError("10X barcodes and cell metadata must match exactly without duplicate IDs")
        feature_map=[];seen=set()
        with _open(gene_file) as f:
            for row in csv.reader(f,delimiter="\t"):
                if not row or not row[0] or row[0] in seen:
                    raise ValueError("10X feature IDs must be unique/nonempty")
                seen.add(row[0])
                if len(row)>=3 and row[2]!="Gene Expression":
                    feature_map.append(None);excluded_non_gene+=1
                else:
                    feature_map.append(len(features));features.append(row[0]);aggregate.append([0]*len(keys))
        with _open(mtx) as f:
            header=f.readline().strip().lower()
            if header not in {"%%matrixmarket matrix coordinate integer general","%%matrixmarket matrix coordinate real general"}:
                raise ValueError("only general coordinate Matrix Market raw counts are supported")
            def content():
                for line in f:
                    if line.strip() and not line.lstrip().startswith("%"):
                        yield line.split()
            lines=content()
            dims=next(lines,None)
            if dims is None or len(dims)!=3:
                raise ValueError("invalid Matrix Market dimensions")
            nf,nc,nnz=map(int,dims)
            if nf!=len(feature_map) or nc!=len(barcodes_list) or nnz<0:
                raise ValueError("Matrix Market dimensions disagree with feature/barcode IDs")
            for entry in lines:
                if len(entry)!=3:
                    raise ValueError("invalid Matrix Market coordinate")
                j,i=int(entry[0])-1,int(entry[1])-1; value=_integer(entry[2]);n_entries+=1
                if not 0<=j<nf or not 0<=i<nc:
                    raise ValueError("Matrix Market coordinate out of range")
                dest=feature_map[j]
                if dest is not None:
                    add(dest,barcodes_list[i],value)
            if n_entries!=nnz:
                raise ValueError("Matrix Market nonzero entry count does not match its header")
    if not features:
        raise ValueError("no gene-expression features")
    cell_types=tuple(sorted({key[0] for key in keys})); outputs=[]
    reserved={output/"pseudobulk.json"}
    for i in range(len(cell_types)):
        reserved.update({output/("celltype_"+str(i+1))/name for name in ("matrix.tsv","metadata.tsv")})
    if {p.resolve() for p in sources}&{p.resolve() for p in reserved}:
        raise ValueError("output filenames would overwrite pseudobulk source inputs")
    output.mkdir(parents=True,exist_ok=True)
    for i,cell_type in enumerate(cell_types):
        indices=[j for j,key in enumerate(keys) if key[0]==cell_type]
        sids=["PB"+str(j+1).zfill(6) for j in indices]
        folder=output/("celltype_"+str(i+1));folder.mkdir(exist_ok=True)
        _write_tsv(folder/"matrix.tsv",("feature_id",*sids),[{"feature_id":fid,**{sid:aggregate[g][j] for sid,j in zip(sids,indices)}} for g,fid in enumerate(features)])
        metadata_rows=[{"sample_id":sid,"donor_id":keys[j][1],"group":keys[j][2],"cell_type":cell_type,
                        "n_cells":counts[keys[j]],"source_samples":json.dumps(sorted(unit_samples[keys[j]])),
                        **dict(zip(extra,covariates[keys[j]]))} for sid,j in zip(sids,indices)]
        _write_tsv(folder/"metadata.tsv",tuple(metadata_rows[0]),metadata_rows)
        outputs.append({"cell_type":cell_type,"directory":folder.name,"n_pseudobulks":len(indices),
                        "donors":sorted({keys[j][1] for j in indices}),"groups":dict(Counter(keys[j][2] for j in indices))})
    report={"status":"PREPARED","input_format":input_format,"n_cells":len(cells),"retained_cells":len(retained_cells),
            "n_genes":len(features),"min_cells":min_cells,"n_pseudobulks":len(keys),"cell_types":outputs,
            "excluded_low_cell_units":[{"cell_type":key[0],"donor_id":key[1],"group":key[2],"n_cells":count} for key,count in sorted(counts.items()) if count<min_cells],
            "excluded_non_gene_features":excluded_non_gene,"source_entries":n_entries,
            "input_sha256":{str(p):_hash(p) for p in sources},
            "counts_layer":counts_layer if input_format=='h5ad' else None,
            "merge_technical_declared":merge_technical,"carried_covariates":list(extra),
            "ignored_cell_metadata_columns":[x for x in h if x not in required and x not in extra],
            "aggregation":"Exact integer sum by donor + condition + cell type; multiple sample libraries require an explicit technical-merge declaration; duplicate MTX coordinates are additive.",
            "scope":"Cells are never independent deletion units. Counts/annotations/QC supplied by caller; cell identity and genotype checks NOT_EVALUATED."}
    _write_json(output/"pseudobulk.json",report)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(prog="extremarank pseudobulk",description=__doc__)
    p.add_argument("matrix",type=Path);p.add_argument("--metadata",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    p.add_argument("--format",choices=("dense","10x","h5ad"),default="dense");p.add_argument("--min-cells",type=int,default=10)
    p.add_argument("--counts-layer",help="explicit raw-count layer for H5AD; X is an explicit option")
    p.add_argument("--merge-technical",action="store_true")
    p.add_argument("--covariate",action="append",default=[]);p.add_argument("--numeric-covariate",action="append",default=[])
    p.add_argument("--audit-design",choices=("paired","welch"));p.add_argument("--audit-model",choices=("limma-voom","edgeR","DESeq2"))
    p.add_argument("--target");p.add_argument("--reference");p.add_argument("--budget",type=int,default=1);p.add_argument("--top-k",type=int,default=20)
    p.add_argument("--direction",choices=("up","down","absolute"),default="up")
    p.add_argument("--min-total-count",type=int,default=10)
    p.add_argument("--min-count-samples",type=int,default=0)
    a=p.parse_args(argv)
    if (a.audit_design or a.audit_model) and (not a.target or not a.reference):
        p.error("auditing requires explicit target/reference")
    try:
        report=aggregate_pseudobulk(a.matrix,a.metadata,a.output,a.format,a.min_cells,
                                  a.covariate+a.numeric_covariate,a.merge_technical,a.counts_layer)
        audits=[]
        for cell in report["cell_types"]:
            folder=a.output/cell["directory"]
            if a.audit_design or a.audit_model:
                from .preparation import prepare_study
                from .models import run_refits
                try:
                    study=prepare_study(folder/"matrix.tsv",folder/"metadata.tsv",a.target,a.reference,a.audit_design or "welch",transform="cpm-log2")
                    maximum=len(study.unit_ids)-(2 if (a.audit_design or "welch")=="paired" else 4)
                    if a.budget>maximum:
                        raise ValueError("requested budget exceeds available independent donor replication")
                    if a.audit_model:
                        result=run_refits(folder/"matrix.tsv",folder/"metadata.tsv",a.target,a.reference,a.audit_model,folder/"refit",
                                          a.top_k,a.budget,a.direction,a.covariate,a.numeric_covariate,paired=a.audit_design=="paired",
                                          min_total_count=a.min_total_count,min_count_samples=a.min_count_samples)
                        audits.append({"cell_type":cell["cell_type"],"status":result["status"],"directory":str(cell["directory"])+"/refit"})
                    else:
                        from .study import main as study_main
                        study_main([str(folder/"matrix.tsv"),"--metadata",str(folder/"metadata.tsv"),"--target",a.target,"--reference",a.reference,
                                    "--design",a.audit_design,"--transform","cpm-log2","--budget",str(a.budget),"--top-k",str(a.top_k),
                                    "--direction",a.direction,"--feature-audit","--output",str(folder/"audit")])
                        result=json.loads((folder/"audit"/"audit.json").read_text())
                        audits.append({"cell_type":cell["cell_type"],"status":result["status"],"directory":str(cell["directory"])+"/audit"})
                except ValueError as e:
                    audits.append({"cell_type":cell["cell_type"],"status":"NOT_EVALUABLE","reason":str(e)})
        report["audits"]=audits
        _write_json(a.output/"pseudobulk.json",report)
    except (OSError,ValueError) as e:
        p.error(str(e))
    print("PREPARED: "+str(report["n_pseudobulks"])+" donor/condition/cell-type units; "+str(a.output.resolve()))
    return 0
