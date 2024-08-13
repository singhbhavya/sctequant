#! /usr/bin/env python
# -*- coding: utf-8 -*-

################################# DOWNLOAD REFS #################################
localrules: download_remote

rule download_remote:
    """ Downloads a remote file and checks the md5sum.
        Filenames, URLs and md5 checksums are configured in the `remotefiles` element of
        the configfile
    """
    output:
        'databases/remotefiles/{f}'
    params:
        url = lambda wildcards: config['downloads'][wildcards.f]['url']
    shell:
        '''
wget {params.url} > {output[0]}
        '''

################################# EXTRACT REFS #################################

rule extract_gtf_gz:
    input:
        'databases/remotefiles/{f}.gtf.gz'
    output:
        'databases/annotations/{f}.gtf'
    shell:
        '''
echo {resources.tmpdir}
gunzip -c {input} > {output}
        '''

rule extract_genome:
    input:
        'databases/remotefiles/{f}.fa.gz'
    output:
        'databases/sequences/{f}.fa'
    shell:
        '''
echo {resources.tmpdir}
gunzip -c {input} > {output}
        '''
