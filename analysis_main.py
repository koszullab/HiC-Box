#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import hic_exp
import time
import socket
import os
from Bio.Restriction import AllEnzymes, RestrictionBatch

toolbox_directory = os.path.dirname(os.path.abspath(__file__))
working_directory = toolbox_directory

def analyze(name, 
                 name_bank, 
                 enzyme,
                 genome_index, 
                 genome_fasta,
                 genome_fastq,
                 tag_length=6, 
                 looping=False, 
                 speed_looping=4, 
                 quality_min=30, 
                 tot_len_read=700,
                 len_paired_wise_fastq=3,
                 paired_wise_fastq=True,
                 bowtie2=os.path.join(working_directory, "bowtie2"), 
                 bank_folder=os.path.join(toolbox_directory, "results"), 
                 ncpu=4):
    
    start_total = time.clock()
    hostname = socket.gethostname()
    print "Host name:", hostname
    ordi = hostname.split('.')[0]
    
    # if ordi == 'renoir':
    #     folder_alignment_toolbox = '/Volumes/Data/HiC_project/alignment_toolbox'
    #     bowtie2 = '/Volumes/Data/HiC_project/alignment_toolbox/bowtie2-2.0.0-beta5/'
    #     out_foldR = '/Volumes/Data/hic_data/27_08_2013/results'
    #     base_folder= '/Volumes/Data/hic_data/27_08_2013'
    #     ncpu = 24
    # else:
    folder_alignment_toolbox = toolbox_directory
#     bowtie2 = os.path.join(working_directory, "bowtie2")
#     out_foldR = os.path.join(toolbox_directory, "results")
    base_folder= working_directory
#     ncpu = 4
    out_foldR = bank_folder
    if not(os.path.exists(out_foldR)):
        os.mkdir(out_foldR)
    
    
    
    
    ########################## PARAMETERS #################################################################################
#     name = 'Vibrio_WT'
#     name_bank = 'Vibrio_209'
    
    folder_a = os.path.join(base_folder,'')
    folder_b = os.path.join(base_folder,'')
    
    if enzyme in AllEnzymes:
        restriction_site = RestrictionBatch([enzyme]).get(enzyme).site # HpaII
    else:
        restriction_site = enzyme
    
    #folder_alignment_toolbox +
#     genome_index = os.path.join(working_directory, "index/"+name)
#     genome_fasta = os.path.join(working_directory, "fasta/"+name+".fa")
    
#     tag_length = 6
#     looping = False
#     speed_looping = 4
#     quality_min = 30
    print name_bank
    print out_foldR
    output_folder = os.path.join(out_foldR,bank_folder)
#     paired_wise_fastq = True
#     len_paired_wise_fastq = 3
#     tot_len_read = 700
    print genome_fastq
    print genome_fastq.split(',')
    if len(genome_fastq.split(','))<=1:
        motif_read_1 = genome_fastq
        print "Reading "+motif_read_1
        if motif_read_1[-1]=='2':
            motif_read_2 = motif_read_1[:-1]+'1'
            print "Reading "+motif_read_2
        elif motif_read_1[-1]=='1':
            motif_read_2 = motif_read_1[:-1]+'2'
            print "Reading "+motif_read_2
        else:
            motif_read_2 = ""
            print "Warning: no second fastq file found"
    else:
        motif_read_1 = genome_fastq.split(',')[0]
        motif_read_2 = genome_fastq.split(',')[1]
    
    print motif_read_1
    print motif_read_2
    #######################################################################################################################
    hic_bank = hic_exp.hic_exp(name_bank, tot_len_read, folder_a, folder_b,motif_read_1,motif_read_2,paired_wise_fastq,
        restriction_site,ncpu,tag_length,genome_index, genome_fasta,bowtie2,looping,quality_min,output_folder,speed_looping,
        len_paired_wise_fastq)
    start = time.clock()
    hic_bank.align()
    hic_bank.pcr_free()
    hic_bank.paired_reads_2_fragments()
    elapsed = (time.clock() - start)
    print 'Paired reads aligned in  ' + str(elapsed)+ ' s'
    print " start computing biases..."
    hic_bank.gc_size_bias()
    print " done."
    print "writing abs weighted contacts"
    hic_bank.fragments_contacts_2_weighted_contacts()
    elapsed = (time.clock() - start_total)
    print "all done in " + str(elapsed)+ " s"
    print "ready for computation"
    return True
