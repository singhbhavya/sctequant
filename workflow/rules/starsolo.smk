######################################################
## Analyze droplet single-cell RNA sequencing data
## with STARsolo

# Setup
import csv
from collections import defaultdict

# Function to get files
def create_comma_separated_list_for_biospecimen(all_dict, biospecimen, read_type):
    if biospecimen not in all_dict:
        return ""

    filenames = all_dict[biospecimen]

    if read_type not in ['R1', 'R2']:
        raise ValueError("read_type must be 'R1' or 'R2'")

    # Create an empty list to store unique filenames
    unique_file_list = []
    
    # Add filenames to the list, avoiding duplicates
    for filename in filenames:
        full_filename = f"/sc/arion/projects/Perturb-map/scRNAseq/scRNAseq_PDAC/Zhou_NatGenetics/Zhou_NatGen_scRNAseq_raw/{filename}_{read_type}_001.fastq.gz"
        if full_filename not in unique_file_list:
            unique_file_list.append(full_filename)

    return ",".join(unique_file_list)


# Get whitelists
rule get_whitelist:
    """
    Get appropriate 10x whitelist of used barcodes
    """
    input:
        "databases/remotefiles/3M-february-2018.txt.gz"
    output:
        config['whitelist']['v3']
    shell:
        '''
        gunzip -c {input} > {output}
        '''

# Run STARSolo 
rule starsolo_alignment:
    """
    Align sequencing reads from a 10x V3 single-cell RNA-seq experiment using STARsolo
    """
    conda:
        "../envs/star.yaml"
    params:
        cDNA = lambda wildcards: create_comma_separated_list_for_biospecimen(all_dict, wildcards.samid, "R2"),
        barcodes = lambda wildcards: create_comma_separated_list_for_biospecimen(all_dict, wildcards.samid, "R1"),
        out_prefix="results/starsolo_algn/{samid}/{samid}_GDC38.",
        tmp="/sc/arion/scratch/singhb22/{samid}_tmp"
    input:
        genome = config['indexes']['star'],
        whitelist = "databases/remotefiles/3M-february-2018.txt"
    output: "results/starsolo_algn/{samid}/{samid}_GDC38.Aligned.sortedByCoord.out.bam","results/starsolo_algn/{samid}/{samid}_GDC38.Solo.out/Gene/filtered/barcodes.tsv"
    benchmark: "benchmarks/star_alignment/{samid}_star_alignment.tsv"
    resources: mem_mb=50000, disk_mb=100000, mem_mib=50000, time_min=1840
    threads: 8
    shell:
        '''
        #--- STARsolo (turned on by --soloType CB_UMI_Simple)
        STAR\
             --outSAMattributes NH HI nM AS CR UR CY UY CB UB GX GN sS sQ sM\
             --readFilesCommand zcat\
             --outSAMtype BAM SortedByCoordinate\
             --outFilterMultimapNmax 100\
             --winAnchorMultimapNmax 100\
             --outMultimapperOrder Random\
             --outSAMmultNmax 1\
             --runRNGseed 777\
             --runThreadN 8\
             --soloUMIlen 12\
             --soloType CB_UMI_Simple\
             --soloCBwhitelist {input.whitelist}\
             --soloCellFilter EmptyDrops_CR\
             --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts\
             --soloUMIfiltering MultiGeneUMI_CR\
             --soloUMIdedup 1MM_CR\
             --soloMultiMappers EM\
             --readFilesIn {params.cDNA} {params.barcodes}\
             --genomeDir {input.genome}\
             --outFileNamePrefix {params.out_prefix}\
             --outTmpDir {params.tmp}\
             --limitBAMsortRAM 31000000000
        '''



