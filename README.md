# Single Cell TE Quantification Snakemake Pipeline with SoloTE

**Prior to running:** 

- If using on Minerva, please create a [Snakemake LSF profile with cookiecutter](https://github.com/Christian-Heyer/snakemake-lsf).
- Download [SoloTE](https://github.com/bvaldebenitom/SoloTE) and modify the path to SoloTE in the relevant rules (or add it to your .bashrc).
- Initialize conda and/or mamba to automatically create the relevant environment per rule.

**To get DAG:** 

```
snakemake --profile lsf  --forceall --dag | dot -Tpdf > dag.pdf
```

**To get rule graph:** 

```
snakemake --profile lsf  --forceall --rulegraph | dot -Tpdf > rulegraph.pdf
```

**To get file graph:** 

```
snakemake --profile lsf  --forceall --filegraph | dot -Tpdf > filegraph.pdf
```

**To run pipeline:**

```
snakemake --profile lsf 
```

**To modify pipeline:**

Modify the sample table in the config directory, and change the appropriate column names (with the sample/run IDs) in the Snakefile. 
