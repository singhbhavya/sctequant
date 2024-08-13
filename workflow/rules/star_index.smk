# -*- coding utf-8 -*-

rule star_index_gdc38_gencode38:
    conda:
        "../envs/star.yaml"
    input:
        genome = "databases/sequences/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        annotation = "databases/annotations/Homo_sapiens.GRCh38.93.gtf"
    output:
       directory("databases/indexes/STAR_GRCh38"),
       "databases/indexes/STAR_GRCh38/SAindex"
    params:
        sjdbOverhang = 74
    threads: 8
    shell:
        """
        STAR\
            --runThreadN {threads}\
            --runMode genomeGenerate\
            --genomeDir {output}\
            --outFileNamePrefix {output}\
            --genomeFastaFiles {input.genome}\
            --sjdbGTFfile {input.annotation}\
            --sjdbOverhang {params.sjdbOverhang}\
            --genomeSAsparseD 3
        """
