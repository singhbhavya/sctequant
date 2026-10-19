######################################################################
## Final outputs: one AnnData (.h5ad) per sample per TE tool, a merged
## AnnData per tool, and a QC table.
######################################################################


def te_output(wc):
    if wc.tool == "solote":
        return rules.solote.output[0].format(sample=wc.sample)
    if wc.tool == "scte":
        return rules.scte.output[0].format(sample=wc.sample)
    return rules.irescue.output[0].format(sample=wc.sample)


rule sample_metadata:
    output: f"{OUT}/{{sample}}/metadata.csv"
    localrule: True
    run:
        SAMPLES.loc[[wildcards.sample]].to_csv(output[0], index=False)


rule make_anndata:
    input:
        genes=rules.starsolo.output.genes_filtered,
        te=te_output,
        families=rules.te_annotation.output.families,
        meta=rules.sample_metadata.output,
    output: f"{OUT}/anndata/{{tool}}/{{sample}}.h5ad"
    params:
        level=SOLOTE_LEVELS[config["solote_level"]],
        chemistry=lambda wc: SAMPLES.loc[wc.sample, "chemistry"],
    resources: mem_mb=16000, runtime=60
    log: f"{OUT}/logs/{{sample}}/anndata_{{tool}}.log"
    conda: "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/make_anndata.py --sample {wildcards.sample} --tool {wildcards.tool} \
            --genes {input.genes} --te {input.te} --families {input.families} \
            --sample-sheet-row {input.meta} --solote-level {params.level} \
            --chemistry {params.chemistry} --genome {GENOME} --out {output} > {log} 2>&1
        """


rule merge_anndata:
    input: expand(f"{OUT}/anndata/{{{{tool}}}}/{{sample}}.h5ad", sample=SAMPLES.index)
    output: f"{OUT}/anndata/{{tool}}_all_samples.h5ad"
    resources: mem_mb=64000, runtime=120
    log: f"{OUT}/logs/merge_anndata_{{tool}}.log"
    conda: "../envs/python.yaml"
    shell: "python workflow/scripts/merge_anndata.py {output} {input} > {log} 2>&1"


rule qc_summary:
    input: expand(rules.starsolo.output.summary, sample=SAMPLES.index)
    output: f"{OUT}/qc/starsolo_summary.csv"
    localrule: True
    conda: "../envs/python.yaml"
    shell: "python workflow/scripts/qc_summary.py {output} {input}"
