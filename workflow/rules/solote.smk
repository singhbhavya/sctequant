
# Run SoloTE 

localrules: soloTE_annotation
rule soloTE_annotation:
    conda:
        "/sc/arion/work/singhb22/test-env/envs/SoloTE"
    output: "databases/annotations/hg38_rmsk_edited.bed"
    shell:
        '''
        python SoloTE_RepeatMasker_to_BED.py -g hg38
        sed 's/^chr//' hg38_rmsk.bed > hg38_rmsk_edited.bed
        '''


rule soloTE:
    conda:
        "/sc/arion/work/singhb22/test-env/envs/SoloTE"
    input: 
        bam = "results/starsolo_algn/{samid}/{samid}_GDC38.Aligned.sortedByCoord.out.bam",
        annotation = "databases/annotations/hg38_rmsk_edited.bed"
    output:
        features = "results/SoloTE/{samid}/{samid}_SoloTE_output/{samid}_locustes_MATRIX/features.tsv",
        barcodes = "results/SoloTE/{samid}/{samid}_SoloTE_output/{samid}_locustes_MATRIX/barcodes.tsv",
        matrix = "results/SoloTE/{samid}/{samid}_SoloTE_output/{samid}_locustes_MATRIX/matrix.mtx"
    resources: time_min=1000, solote_jobs=1
    params:
        prefix = "{samid}",
        out_dir = "results/SoloTE/{samid}"
    threads: 8
    shell:
        '''
        #--- SoloTE
        python /sc/arion/projects/Perturb-map/scRNAseq/scRNAseq_PDAC/Zhou_NatGenetics/sctequant/packages/SoloTE/SoloTE_pipeline.py\
        --threads {threads}\
        --bam {input.bam}\
        --teannotation {input.annotation}\
        --outputprefix {params.prefix}\
        --outputdir {params.out_dir}
        '''

rule soloTE_complete:
    input: rules.soloTE.output
    output: touch("results/completed/{samid}_completed.txt")

rule soloTE_features_processing:
    input: 
       	features = "results/SoloTE/{samid}/{samid}_SoloTE_output/{samid}_locustes_MATRIX/features.tsv",
       	gtf = "databases/annotations/Homo_sapiens.GRCh38.93.gtf"
    output: "results/SoloTE/{samid}/{samid}_SoloTE_output/{samid}_locustes_MATRIX/features.tsv.gz"
    shell:
        """
        python workflow/scripts/solote_features_to_standard.py {input.gtf} {input.features}
        """

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

rule SoloTE_gzip:
    input: 
        barcodes = "results/SoloTE/{samid}/{samid}_SoloTE_output/{samid}_locustes_MATRIX/barcodes.tsv",
        matrix = "results/SoloTE/{samid}/{samid}_SoloTE_output/{samid}_locustes_MATRIX/matrix.mtx"
    output:
        barcodes = "results/SoloTE/{samid}/{samid}_SoloTE_output/{samid}_locustes_MATRIX/barcodes.tsv.gz",
        matrix = "results/SoloTE/{samid}/{samid}_SoloTE_output/{samid}_locustes_MATRIX/matrix.mtx.gz"
    shell:
        '''
        gzip {input.barcodes}
        gzip {input.matrix}
        '''

rule cleanup:
    input: rules.SoloTE_gzip.output, rules.soloTE_features_processing.output
    output: "results/cleanup/{samid}.cleanup.done.txt"
    run:
        # Get the cDNA and barcode file lists as comma-separated strings
        cDNA_files_str = create_comma_separated_list_for_biospecimen(all_dict, wildcards.samid, 'R1')
        barcode_files_str = create_comma_separated_list_for_biospecimen(all_dict, wildcards.samid, 'R2')

        # Split the comma-separated strings into lists
        cDNA_files = cDNA_files_str.split(",")
        barcode_files = barcode_files_str.split(",")

        # Combine the lists of files to delete
        files_to_delete = cDNA_files + barcode_files

        # Debugging: Print the files to be deleted
        print("cDNA files:")
        for f in cDNA_files:
            print(f)
        print("Barcode files:")
        for f in barcode_files:
            print(f)

        # Validate and delete the files
        for f in files_to_delete:
            if os.path.exists(f):
                print(f"Deleting file: {f}")
                shell(f"rm {f}")
            else:
                print(f"File not found, skipping: {f}")

        # Create the output file indicating cleanup is done
        with open(output[0], 'w') as f:
            f.write("Cleanup completed for samid {samid}\n".format(samid=wildcards.samid))

