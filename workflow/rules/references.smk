######################################################################
## Reference files: genome, gene annotation, TE annotation, STAR index,
## scTE index, 10x whitelists and the SoloTE code.
## These are built once in reference_dir and reused by every run.
######################################################################

# Downloads run on the machine that runs Snakemake (Minerva compute
# nodes may not have internet access).
localrules: download_genome, download_gtf, download_rmsk, download_whitelist, get_solote, prepare


rule download_genome:
    output: f"{REF}/genome/genome.fa"
    params: url=config["genome"]["fasta_url"]
    shell: "curl -fsSL --retry 3 {params.url} | gunzip -c > {output}.tmp && mv {output}.tmp {output}"


rule download_gtf:
    output: f"{REF}/genome/genes.gtf"
    params: url=config["genome"]["gtf_url"]
    shell: "curl -fsSL --retry 3 {params.url} | gunzip -c > {output}.tmp && mv {output}.tmp {output}"


rule download_rmsk:
    output: f"{REF}/te/{GENOME}.fa.out.gz"
    params: url=config["genome"]["rmsk_url"]
    shell: "curl -fsSL --retry 3 -o {output}.tmp {params.url} && mv {output}.tmp {output}"


rule download_whitelist:
    output: f"{REF}/whitelists/{{whitelist}}"
    params:
        url=lambda wc: next(c["whitelist_url"] for c in CHEMISTRIES.values() if c["whitelist"] == wc.whitelist)
    shell: "curl -fsSL --retry 3 {params.url} | gunzip -c > {output}.tmp && mv {output}.tmp {output}"


rule get_solote:
    output: f"{REF}/tools/SoloTE/SoloTE_pipeline.py"
    params:
        repo=config["solote_repo"],
        commit=config["solote_commit"],
        dir=f"{REF}/tools/SoloTE",
    shell:
        """
        rm -rf {params.dir}
        git clone -q {params.repo} {params.dir}
        git -C {params.dir} checkout -q {params.commit}
        """


rule te_annotation:
    input: rules.download_rmsk.output
    output:
        subfamily_bed=f"{REF}/te/{GENOME}_te_subfamily.bed",
        solote_bed=f"{REF}/te/{GENOME}_te_solote.bed",
        families=f"{REF}/te/{GENOME}_te_families.tsv",
    conda: "../envs/python.yaml"
    resources: mem_mb=8000, runtime=60
    shell:
        """
        python workflow/scripts/make_te_annotation.py --rmsk {input} \
            --subfamily-bed {output.subfamily_bed} \
            --solote-bed {output.solote_bed} \
            --families {output.families}
        """


rule star_index:
    input:
        fasta=rules.download_genome.output,
        gtf=rules.download_gtf.output,
    output: f"{REF}/star_index/SAindex"
    params: dir=f"{REF}/star_index"
    threads: 12
    resources: mem_mb=64000, runtime=360
    log: f"{REF}/logs/star_index.log"
    conda: "../envs/star.yaml"
    shell:
        """
        rm -rf {params.dir}/_tmp
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.dir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --outTmpDir {params.dir}/_tmp \
            --outFileNamePrefix {params.dir}/ > {log} 2>&1
        """


rule scte_index:
    input:
        bed=rules.te_annotation.output.subfamily_bed,
        gtf=rules.download_gtf.output,
    output: f"{REF}/scte/{GENOME}.exclusive.idx"
    params: prefix=f"{REF}/scte/{GENOME}"
    resources: mem_mb=32000, runtime=180
    log: f"{REF}/logs/scte_index.log"
    conda: "../envs/scte.yaml"
    shell: "scTE_build -te {input.bed} -gene {input.gtf} -m exclusive -o {params.prefix} > {log} 2>&1"



rule prepare:
    """Download everything that needs the internet (run on a login node)."""
    input:
        rules.download_genome.output,
        rules.download_gtf.output,
        rules.download_rmsk.output,
        lambda wc: sorted({whitelist_path(s) for s in SAMPLES.index}),
        lambda wc: rules.get_solote.output if config["tools"].get("solote") else [],
