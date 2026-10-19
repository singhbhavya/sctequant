######################################################################
## TE quantification. Each tool reads the same STARsolo BAM.
##   SoloTE  - locus-aware; counts TE reads not assigned to genes
##   scTE    - TE subfamilies; reads in gene exons go to genes first
##   IRescue - TE subfamilies; shares multi-mapping UMIs between TEs (EM)
######################################################################

BAM = f"{OUT}/{{sample}}/starsolo/Aligned.sortedByCoord.out.bam"
BAI = f"{BAM}.bai"

SOLOTE_LEVELS = {"subfamily": "subfamilytes", "family": "familytes", "class": "classtes", "locus": "locustes"}


rule solote:
    input:
        bam=BAM,
        bai=BAI,
        bed=rules.te_annotation.output.solote_bed,
        code=rules.get_solote.output,
    output:
        directory(f"{OUT}/{{sample}}/solote/{{sample}}_SoloTE_output"),
    params:
        outdir=f"{OUT}/{{sample}}/solote",
    threads: config["te_threads"]
    resources: mem_mb=48000, runtime=12 * 60
    log: f"{OUT}/logs/{{sample}}/solote.log"
    benchmark: f"{OUT}/benchmarks/{{sample}}/solote.tsv"
    conda: "../envs/solote.yaml"
    shell:
        """
        bam=$(realpath {input.bam}); bed=$(realpath {input.bed})
        code=$(realpath {input.code}); log=$(realpath -m {log})
        mkdir -p {params.outdir} && cd {params.outdir}
        python $code --threads {threads} --bam $bam --teannotation $bed \
            --outputprefix {wildcards.sample} --outputdir . > $log 2>&1
        rm -rf {wildcards.sample}_SoloTE_temp
        """


rule scte:
    input:
        bam=BAM,
        bai=BAI,
        index=rules.scte_index.output,
    output: f"{OUT}/{{sample}}/scte/{{sample}}.h5ad"
    params:
        outdir=f"{OUT}/{{sample}}/scte",
    threads: 4  # scTE uses ~10 GB of memory per thread
    resources: mem_mb=48000, runtime=12 * 60
    log: f"{OUT}/logs/{{sample}}/scte.log"
    benchmark: f"{OUT}/benchmarks/{{sample}}/scte.tsv"
    conda: "../envs/scte.yaml"
    shell:
        """
        bam=$(realpath {input.bam}); index=$(realpath {input.index}); log=$(realpath -m {log})
        mkdir -p {params.outdir} && cd {params.outdir}
        scTE -i $bam -o {wildcards.sample} -x $index -p {threads} \
            -CB CB -UMI UB --hdf5 True --min_genes 100 --expect-cells 50000 > $log 2>&1
        """


rule irescue:
    input:
        bam=BAM,
        bai=BAI,
        bed=rules.te_annotation.output.subfamily_bed,
        cells=rules.starsolo.output.genes_filtered,
    output: directory(f"{OUT}/{{sample}}/irescue/counts")
    params:
        outdir=f"{OUT}/{{sample}}/irescue",
    threads: config["te_threads"]
    resources: mem_mb=32000, runtime=12 * 60
    log: f"{OUT}/logs/{{sample}}/irescue.log"
    benchmark: f"{OUT}/benchmarks/{{sample}}/irescue.tsv"
    conda: "../envs/irescue.yaml"
    shell:
        """
        irescue -b {input.bam} -r {input.bed} -w {input.cells}/barcodes.tsv \
            -c CB -u UB -p {threads} -o {params.outdir} > {log} 2>&1
        """
