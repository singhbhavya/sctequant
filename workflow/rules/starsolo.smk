######################################################################
## STARsolo: align reads, call cells, count genes.
## The BAM keeps multi-mapping reads (needed for TE quantification) and
## carries the CB (cell barcode), UB (UMI) and GN (gene) tags the TE
## tools read.
######################################################################


rule starsolo:
    input:
        r1=lambda wc: FASTQS[wc.sample][0],
        r2=lambda wc: FASTQS[wc.sample][1],
        index=rules.star_index.output,
        whitelist=lambda wc: whitelist_path(wc.sample),
    output:
        bam=f"{OUT}/{{sample}}/starsolo/Aligned.sortedByCoord.out.bam",
        genes_filtered=directory(f"{OUT}/{{sample}}/starsolo/Solo.out/Gene/filtered"),
        summary=f"{OUT}/{{sample}}/starsolo/Solo.out/Gene/Summary.csv",
    params:
        prefix=f"{OUT}/{{sample}}/starsolo/",
        index_dir=f"{REF}/star_index",
        # STAR takes the cDNA read first, then the barcode read.
        cdna=lambda wc: ",".join(FASTQS[wc.sample][1]),
        barcode=lambda wc: ",".join(FASTQS[wc.sample][0]),
        cb_len=lambda wc: chem(wc.sample)["cb_len"],
        umi_start=lambda wc: chem(wc.sample)["cb_len"] + 1,
        umi_len=lambda wc: chem(wc.sample)["umi_len"],
        strand=lambda wc: chem(wc.sample)["strand"],
        sort_ram=int(config["star"]["mem_mb"] * 0.4) * 1000000,
        extra=config["star"]["extra"],
    threads: config["star"]["threads"]
    resources:
        mem_mb=config["star"]["mem_mb"],
        runtime=24 * 60,
    log: f"{OUT}/logs/{{sample}}/starsolo.log"
    benchmark: f"{OUT}/benchmarks/{{sample}}/starsolo.tsv"
    conda: "../envs/star.yaml"
    shell:
        """
        rm -rf {params.prefix}_STARtmp
        STAR --runThreadN {threads} \
            --genomeDir {params.index_dir} \
            --readFilesIn {params.cdna} {params.barcode} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist {input.whitelist} \
            --soloCBstart 1 --soloCBlen {params.cb_len} \
            --soloUMIstart {params.umi_start} --soloUMIlen {params.umi_len} \
            --soloBarcodeReadLength 0 \
            --soloStrand {params.strand} \
            --soloFeatures Gene GeneFull \
            --soloCellFilter EmptyDrops_CR \
            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloUMIfiltering MultiGeneUMI_CR \
            --soloUMIdedup 1MM_CR \
            --clipAdapterType CellRanger4 \
            --outFilterScoreMin 30 \
            --outFilterMultimapNmax 100 \
            --winAnchorMultimapNmax 100 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --limitBAMsortRAM {params.sort_ram} \
            {params.extra} > {log} 2>&1
        """


rule index_bam:
    input: rules.starsolo.output.bam
    output: f"{OUT}/{{sample}}/starsolo/Aligned.sortedByCoord.out.bam.bai"
    threads: 4
    resources: mem_mb=8000, runtime=120
    conda: "../envs/star.yaml"
    shell: "samtools index -@ {threads} {input}"
