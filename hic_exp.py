# coding: utf-8
# -*- coding: utf-8 -*-_
_author__ = 'hervemn'
import os
import glob
import fastq_alignment
import shutil
import hic_analysis
import gephi_lib

class hic_exp():
    def __init__(self, name, tot_len_read, folder_a, folder_b,motif_read_1,motif_read_2,paired_wise_fastq,restriction_site,ncpu,len_tag,genome_index, genome_fasta, bowtie2,looping,
                 quality_min,output_folder,speed,len_paired_wise_fastq):
        self.len_paired_wise_fastq = len_paired_wise_fastq
        self.name = name
        self.tot_len_read = tot_len_read ## fragment size sent to the sequencer!
        self.seq = restriction_site
        self.paired_wise_fastq = paired_wise_fastq
        self.paired_files = zip(glob.glob(os.path.join(folder_a,motif_read_1)),glob.glob(os.path.join(folder_b,motif_read_2)))
        self.ncpu = ncpu
        self.output_folder = output_folder
        if not(os.path.exists(output_folder)):
            os.mkdir(output_folder)
        self.folder_analysis =os.path.join(output_folder,'analysis/')
        if not(os.path.exists(self.folder_analysis)):
            os.mkdir(self.folder_analysis)
        if os.path.exists(os.path.join(output_folder,name+'_all_paired_pcr_free.sam')):
            self.contacts_pcr_free = os.path.join(output_folder,name+'_all_paired_pcr_free.sam')
        self.genome_index = genome_index
        self.genome_fasta = genome_fasta
        self.bowtie2 = bowtie2
        self.looping = looping
        self.quality = quality_min
        self.len_tag = len_tag
        self.speed_looping = speed

        #### define restriction map ####
        self.file_contig_txt = os.path.join(output_folder,'list_contig_names.txt')
        os.system(os.path.join(bowtie2,'./bowtie2-inspect')+' '+ genome_index+' -n>'+self.file_contig_txt)
        
    
        self.restriction_site_list = os.path.join(output_folder,'restriction_site_list.txt')
        self.fragments_list =os.path.join(self.folder_analysis,'fragments_list.txt')
        self.dict_position_restriction_sites  = dict()
        self.dict_contigs = dict()
        print "making restriction map .."
        fastq_alignment.restriction_map(self.genome_fasta, self.seq, self.restriction_site_list,
            self.fragments_list, self.dict_position_restriction_sites, self.dict_contigs)
        
        self.fragments_contacts_file = os.path.join(self.output_folder,'fragments_contacts.txt')
        self.fragments_hetero_contacts_file = os.path.join(self.output_folder,'fragments_hetero_contacts.txt')
        self.fragments_contacts_file_abs = os.path.join(self.output_folder,'fragments_contacts_abs.txt')
        self.contacts_pcr_free = os.path.join(output_folder,name+'_all_paired_pcr_free.sam')
        
        file_list_contig = open(self.file_contig_txt,'r')
        list_contigs = []
        print "filling list of contigs .."
        for contig in file_list_contig:
            list_contigs.append(contig[:-1])
        self.list_contigs = list_contigs
        print list_contigs
        dict_fragments = dict()
        for ele in list_contigs:
            dict_fragments[ele] = []
        i = 0
        file_list_contig.close()
        frag_list = open(self.fragments_list,'r')
        frag_list.readline()
        print "filling dictionnary of fragments ..."
        while 1:
            line = frag_list.readline()
            if not line:
                frag_list.close()
                break
            a_tmp = line.split('\t')
            dict_fragments[a_tmp[1]].append(int(a_tmp[0]))

        dict_cumul_length = dict()
        dict_cumul_length[list_contigs[0]] = 0
        n_contigs = len(list_contigs)

        for i in xrange(1,n_contigs):
#            chunk = dict_cumul_length[list_contigs[i-1]]+len(dict_fragments[list_contigs[i-1]]) modif herve 4 06 2012
            chunk = dict_cumul_length[list_contigs[i-1]]+int(dict_fragments[list_contigs[i-1]][-1])
            print 'length '+ list_contigs[i] + ' : ' + str(chunk)
            dict_cumul_length[list_contigs[i]] = chunk
        self.dict_cumul_length = dict_cumul_length
        dict_fragments = dict()
        frag_list = open(self.fragments_list,'r')
        i = 0
        for line in frag_list:
            if i>0:
                data = line.split()
                id_frag = data[0]+'-'+data[1]
                dict_fragments[id_frag] = {'size':data[4],'gc_content':data[5],'start':data[2]}
            i = i+1
        frag_list.close()
        self.dict_fragments = dict_fragments
        self.file_cumul_length = os.path.join(self.folder_analysis,'info_contigs.txt')
        file_cumul_length = open(self.file_cumul_length,'w')
        file_cumul_length.write("%s\t%s\t%s\t%s\n" %('contig','length_kb','n_frags','cumul_length'))
        for ele in list_contigs:
            file_cumul_length.write("%s\t%s\t%s\t%s\n" %(ele,self.dict_contigs[ele]['length_kb'],
                                                         self.dict_contigs[ele]['n_frags'], self.dict_cumul_length[ele]))
        file_cumul_length.close()
        self.fragments_contacts_files_weighted = os.path.join(self.folder_analysis,'fragments_contacts_weighted.txt')
        self.fragments_abs_contacts_files_weighted = os.path.join(self.folder_analysis,'abs_fragments_contacts_weighted.txt')

    def align(self,):
        quality = self.quality
        gen_index = self.genome_index
        bowtie2 = self.bowtie2
        out_dir = self.output_folder
        experience_name = self.name
        i = 0
        ncpu_pp = 1
        speed = self.speed_looping
        ncpu_bowtie = self.ncpu/ncpu_pp
        looping = self.looping
        len_tag = self.len_tag
        paired_wise_fastq = self.paired_wise_fastq
        jobs = []
        ppservers = ()
        print self.bowtie2
        print self.paired_files
        for paired_file in self.paired_files:
            print paired_file[0] + '  ' + paired_file[1]
            id = str(i)
            in_a = paired_file[0]
            in_b = paired_file[1]
            fastq_alignment.bowtie_fastq(bowtie2, in_a, in_b, gen_index, out_dir,id,ncpu_bowtie, looping,
                                         quality, len_tag, paired_wise_fastq, speed, self.len_paired_wise_fastq)
#            jobs.append(job_server.submit(fastq_alignment.bowtie_fastq, (bowtie2, in_a, in_b, gen_index, out_dir,id,ncpu_bowtie, looping,quality,len_tag,paired_wise_fastq,speed)
#                , (fastq_alignment.sam_filter,), ("fastq_alignment",)))
            i = i +1
            print str(i) + " jobs launched over "+ str(ncpu_bowtie) + " cores"
#        for job in jobs:
#            print "la?..."
#            job()
        
    def pcr_free(self):
        pcr_free = fastq_alignment.pcr_amplification_extract(self.output_folder,self.name)
        self.contacts_pcr_free = pcr_free
        print 'pcr amplification detection done!'

    def paired_reads_2_fragments(self,):
        print 'building fragments contacts from data'
        print "fragment contact file = "
        print self.fragments_contacts_file
        fastq_alignment.paired_reads_2_frag_contacts(self.tot_len_read,self.contacts_pcr_free,self.len_tag,self.dict_position_restriction_sites, self.fragments_contacts_file)
        fastq_alignment.remove_self_fragments_contacts(self.fragments_contacts_file,self.fragments_hetero_contacts_file)
        ## creer index absolu des fragments!! ###

    def load_exp_2_db(self,):
        fastq_alignment.fragments_2_db(self.name, self.fragments_list)
        fastq_alignment.load_fragments_contacts_2_db(self.name,self.fragments_contacts_file)

    def rel_frag_2_abs_frag(self,):
        fastq_alignment.rel_frag_2_abs_frag(self.list_contigs,self.fragments_contacts_file,self.fragments_contacts_file_abs,self.fragments_list)

    def show_contact_matrix(self,):
        hic_analysis.draw_matrix(self.folder_analysis,self.fragments_abs_contacts_files_weighted,self.dict_fragments)

    def gc_size_bias(self):
        mat_gc, mat_size,steps_gc,steps_length =  hic_analysis.gc_size_bias(self.folder_analysis,self.dict_fragments, self.fragments_contacts_file)
        self.mat_gc_bias = mat_gc
        self.mat_length_bias = mat_size
        self.steps_length = steps_length
        self.steps_gc = steps_gc

    def data2gephi(self):
        gephi_lib.data2gexf(self.folder_analysis,self.fragments_list,self.fragments_contacts_file_abs,self.dict_cumul_length)

    def contact_vs_distance(self):
        hic_analysis.contact_vs_gen_distance(self.folder_analysis,self.dict_fragments,self.fragments_contacts_file,self.dict_contigs)

    def fragments_contacts_2_weighted_contacts(self,):
        self.fragments_contacts_files_weighted = os.path.join(self.folder_analysis,'fragments_contacts_weighted.txt')
        self.fragments_abs_contacts_files_weighted = os.path.join(self.folder_analysis,'abs_fragments_contacts_weighted.txt')
        hic_analysis.fragments_contacts_2_weighted_contacts(self.dict_cumul_length,
            self.dict_fragments,
            self.fragments_abs_contacts_files_weighted,
            self.fragments_contacts_files_weighted,
            self.fragments_contacts_file,
            self.mat_gc_bias,
            self.mat_length_bias,
            self.steps_gc,
            self.steps_length)
