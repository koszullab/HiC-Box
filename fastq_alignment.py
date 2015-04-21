# coding: utf-8
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 18:01:51 2012

@author: hervemn
"""

import os
import shutil
import time
import pymongo
import re
import numpy
import sys
from progressbar import ProgressBar

def restriction_map(genome_fasta, restriction_site_seq, restriction_site_list,
                    fragments_list,dict_Pos_restriction,dict_contigs):
    f = open(genome_fasta, 'r')
    handle_fragments_list = open(fragments_list, 'w')
    handle_site_list = open(restriction_site_list, 'w')
    D = dict()
    chrom_list = []
    all_lines = f.readlines()
    id_chrom = all_lines[0][1:-1].replace('\r', '').replace('\n', '')
    # chrom_list.append(id_chrom)
    D[id_chrom] = ''
    start = 1
    start_time = time.clock()

    # for i in xrange(1, len(all_lines)):
    #     if all_lines[i][0] == '>':
    #         print id_chrom
    #         chrom_list.append(id_chrom)
    #         D[id_chrom] = ''.join(all_lines[start:i])
    #         dict_Pos_restriction[id_chrom] = []
    #         start = i
    #         id_chrom = all_lines[i][1:-1]
    #         D[id_chrom] = ''
    for line in all_lines:
        if line[0] == '>':
            id_chrom = line[1:-1].replace('\r', '').replace('\n', '')
            D[id_chrom] = ''
            chrom_list.append(id_chrom)
            dict_Pos_restriction[id_chrom] = []
            print 'id_chrom =', id_chrom
        else:
            tmp_str = D[id_chrom] + line
            D[id_chrom] = tmp_str

    print id_chrom
    # D[id_chrom] = ''.join(all_lines[start:-1])
    # chrom_list.append(id_chrom)
#    chrom_list = D.keys()
#    chrom_list.sort()
    print chrom_list

    handle_fragments_list.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ('id','chrom',
                                                              'start_pos','end_pos','size',
                                                              'gc_content'))

    for chrom in chrom_list:
        tmp_str = D[chrom].replace('\n', '').replace('\r','')
        D[chrom] = tmp_str.upper()

    print 'writing files and fiding restriction sites'
    for chrom in chrom_list:
        dict_contigs[chrom] = dict()
        dict_contigs[chrom]['length_kb']  = len(D[chrom])
        len_cont = dict_contigs[chrom]['length_kb']
        print 'contig: ' + chrom
        print 'length = ', len_cont
        starts = [match.start() for match in re.finditer(re.escape(restriction_site_seq), D[chrom])]

        dict_Pos_restriction[chrom] = numpy.array(starts)
        pos_init = 0
        id_frag = 0

        for pos in starts:
            id_frag += 1
            handle_site_list.write("%s\t%s\n" % (chrom, str(pos)))
            size_frag = pos - pos_init
            if not((id_frag ==1) and (size_frag ==0)):
                gc_content = float(D[chrom][pos_init:pos].count('G')+D[chrom][pos_init:pos].count('C'))/float(size_frag)
                handle_fragments_list.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (str(id_frag),chrom,
                                                                          str(pos_init),str(pos),str(size_frag),
                                                                          str(gc_content)))
            else:
                id_frag = 0
            pos_init = pos

        pos = len(D[chrom])
        size_frag = pos - pos_init
        if size_frag >0:
            id_frag +=  1
            handle_site_list.write("%s\t%s\n" % (chrom,str(pos_init)))
            gc_content = float(D[chrom][pos_init:pos].count('G')+D[chrom][pos_init:pos].count('C'))/float(size_frag)
            handle_fragments_list.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (str(id_frag),chrom,
                                                                      str(pos_init),str(pos),str(size_frag),
                                                                      str(gc_content)))
        dict_contigs[chrom]['n_frags'] = id_frag
    handle_site_list.close()
    handle_fragments_list.close()
    elapsed = (time.clock() - start_time)
    print 'Restriction map generated in ' + str(elapsed)+ ' s'

def sam_filter(sam,sam_aligned,fastq_unaligned_1,fastq_unaligned_2,quality,len_tag,dict_tag,seq_len):
    import os
    input_a = open(sam,'r')
    print "sam file = ", sam
    handle_aligned = open(sam_aligned,'w')
    handle_unaligned_1 = open(fastq_unaligned_1,'w')
    handle_unaligned_2 = open(fastq_unaligned_2,'w')    
    i= 0
#    if len_tag>0:
#        handle_aligned.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('_id','contig_a','pos_a','or_read_a','tag_a','len_seq_a',
#                                                                       'contig_b','pos_b','or_read_b','tag_b','len_seq_b'))
#    else:
#        handle_aligned.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('_id','contig_a','pos_a','or_read_a','len_seq_a',
#                                                                       'contig_b','pos_b','or_read_b','len_seq_b'))
    while 1:
        line_a = input_a.readline()
        if not line_a:
            handle_aligned.close()
            handle_unaligned_1.close()
            handle_unaligned_2.close()
            break
        if i%2==0:
            read_a = line_a
        elif i%2==1:
            read_b = line_a
            tmp_a = read_a.split('\t')
            tmp_b = read_b.split('\t')
            if ((tmp_a[2] == '*' ) or (tmp_b[2] == '*') or (int(tmp_a[4])<quality)
                                                    or (int(tmp_b[4])<quality) or (int(tmp_a[1]) in (73,89)) or (int(tmp_b[1]) in (137,153)) ):

                handle_unaligned_1.write("@%s\n%s\n+%s\n%s\n" % (tmp_a[0],
                                                                 tmp_a[9],tmp_a[0],tmp_a[10]))
                handle_unaligned_2.write("@%s\n%s\n+%s\n%s\n" % (tmp_b[0],
                                                                 tmp_b[9],tmp_b[0],tmp_b[10]))       
            else:
                guess_or = int(tmp_b[1])-int(tmp_a[1]) -64
                if guess_or == 16:
                    or_read_a = 'c'
                    or_read_b = 'w'
                elif guess_or == -16:
                    or_read_a = 'w'
                    or_read_b = 'c'
                elif (guess_or == 0) and ( (int(tmp_a[4]) == 65) or (int(tmp_a[4]) == 67) ):
                    or_read_a = 'w'
                    or_read_b = 'w'
                elif (guess_or == 0) and ( (int(tmp_a[4]) == 113) or (int(tmp_a[4]) == 115) ):
                    or_read_a = 'c'
                    or_read_b = 'c'
                elif (guess_or == 0) and ( (int(tmp_a[1]) == 65) and (int(tmp_b[1]) == 129) ):
                    or_read_a = 'w'
                    or_read_b = 'w'
                elif (guess_or == 0) and ( (int(tmp_a[1]) == 113) and (int(tmp_b[1]) == 177) ):
                    or_read_a = 'c'
                    or_read_b = 'c'
                else:
                    print "!!!probleme!!!!"
                    print tmp_a[1]
                    print tmp_b[1]
                if len_tag>0:
                    data_paired_read = dict_tag[tmp_a[0]]
                    tag_a = data_paired_read[0]
                    tag_b = data_paired_read[1]
                    id_p = tmp_a[2] + '-' + tmp_a[3] + '-' + tag_a + '-' + tmp_b[2] + '-' + tmp_b[3] + '-' + tag_b
                    handle_aligned.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                         (id_p,tmp_a[2],tmp_a[3],or_read_a,tag_a,str(seq_len),
                                          tmp_b[2],tmp_b[3],or_read_b,tag_b,str(seq_len),
                                          tmp_a[1],tmp_a[4],tmp_a[5],tmp_b[1],tmp_b[4],tmp_b[5]))
                else:
                    id_p = tmp_a[2] + '-' + tmp_a[3] + '-' + tmp_b[2] + '-' + tmp_b[3]
                    handle_aligned.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                                         (id_p,tmp_a[2],tmp_a[3],or_read_a,str(seq_len),
                                          tmp_b[2],tmp_b[3],or_read_b,str(seq_len),
                                         tmp_a[1],tmp_a[4],tmp_a[5],tmp_b[1],tmp_b[4],tmp_b[5]))
        i = i+1

        
def bowtie_fastq(bowtie2, in_a, in_b, gen_index, out_dir,id,n_cpu, looping,quality,len_tag,paired_wise_fastq,speed,len_paired_wise_fastq):
    print "here we go ..." 
    import shutil
    import os
    trimming = ''
    ### extract tag ###################################################################################
    tmp_dir = os.path.join(out_dir,'tmp'+id)
    file_out_aligned = os.path.join(tmp_dir,id+'_aligned.sam')
    fastq_out_unaligned_1 = os.path.join(tmp_dir,id+'_A_fastq.txt')
    fastq_out_unaligned_2 = os.path.join(tmp_dir,id+'_B_fastq.txt')

    if not(os.path.exists(file_out_aligned) and os.path.exists(fastq_out_unaligned_1)
           and os.path.exists((fastq_out_unaligned_2))):
        dict_tag = dict()
        if len_tag>0:
            trimming = ' -5 '+ str(len_tag)
            input_a = open(in_a,'r')
            input_b = open(in_b,'r')
            i = 0
            while 1:
                line_a = input_a.readline()
                line_b = input_b.readline()
                if not line_a:
                    input_a.close()
                    input_b.close()
                    print "tag dictionnary created."
                    break
                if (i%4==0):
                    if paired_wise_fastq:
                        title_a = line_a[1:-len_paired_wise_fastq]
                        title_a = title_a.split(' ')[0]
                        title_b = line_b[1:-len_paired_wise_fastq]
                        title_b = title_b.split(' ')[0]
                    else:
                        title_a = line_a[1:-1]
                        title_b = line_b[1:-1]
                elif (i%4==1):
                    tag_a = line_a[:len_tag]
                    tag_b = line_b[:len_tag]
    #                print "tag read =  ",title_a
                    if (title_a == title_b):
                        dict_tag[title_a] = (tag_a,tag_b)
                    else:
                        print "problem: files unsychronized!!!!!!!!!!"
                        print "Title a is equal to "+str(title_a)
                        print "Title b is equal to "+str(title_b)
                        break
                i = i+1

        #################################################################################################
        #    parameters = ' --sam-no-hd --sam-no-sq --very-sensitive-local --minins 0 --maxins 5 -M 5 -p' + str(n_cpu)
        parameters = ' --sam-no-hd --sam-no-sq --quiet --very-sensitive --maxins 5 -p ' + str(n_cpu)
        print 'ncpu = ',n_cpu
        genome_index = ' -x '+ gen_index
        reads_a = ' -1 ' + in_a
        reads_b = ' -2 ' + in_b
        sam_out = os.path.join(out_dir,id+'.sam')
        ## check que les directories existent!!)
        if not(os.path.exists(os.path.join(out_dir,'tmp'+id))):
            os.mkdir(os.path.join(out_dir,'tmp'+id))
        tmp_dir = os.path.join(out_dir,'tmp'+id)
        file_out = os.path.join(tmp_dir,id+'.sam')
        handle = open(in_a,'r')
        handle.next()
        seq_len = len(handle.next()) - len_tag
    #    tr3 =
    #    trimming_3 = ' -3 '+str(tr3)
        if not(os.path.exists(file_out)):
            print 'start bowtie2'
            bowtie_align = os.path.join(bowtie2,'./bowtie2')
            print bowtie_align + parameters + genome_index + trimming+ reads_a + reads_b + ' -S '+ file_out+'>1'
            os.system(bowtie_align + parameters + genome_index + trimming+reads_a + reads_b + ' -S '+ file_out+'>1')
            print 'leaves bowtie2'
        else:
            print "data already aligned..."
        ############ extract the paired reads that did not align from the list #########
        file_out_aligned = os.path.join(tmp_dir,id+'_aligned.sam')
        fastq_out_unaligned_1 = os.path.join(tmp_dir,id+'_A_fastq.txt')
        fastq_out_unaligned_2 = os.path.join(tmp_dir,id+'_B_fastq.txt')

        if not(os.path.exists(file_out_aligned) and os.path.exists(fastq_out_unaligned_1)
               and os.path.exists(fastq_out_unaligned_2)):
            sam_filter(file_out,file_out_aligned,fastq_out_unaligned_1,fastq_out_unaligned_2,quality,len_tag,dict_tag,seq_len)
            ############ writing output file or startin looping ############################
            shutil.copy(file_out_aligned,sam_out)
        if not looping:
             print 'Alignment done...'
        else:
            print 'Start looping...'
            for trim_3 in range(1,seq_len - 20,speed):
                seq_len_loop = seq_len -trim_5
                print 'trimming '+str(trim_3) +' bp ...'
                reads_a = ' -1 ' + fastq_out_unaligned_1
                reads_b = ' -2 ' + fastq_out_unaligned_2
                file_out = os.path.join(tmp_dir,'tmp_out_loop.sam')
                trimming =' -3 '+ str(trim_3)
                os.system(bowtie2 + parameters + genome_index + trimming + reads_a + reads_b +' -S '+ file_out+'>1')
                file_out_aligned_loop = os.path.join(tmp_dir,'tmp_out_loop_aligned.sam')
                sam_filter(file_out,file_out_aligned_loop,fastq_out_unaligned_1,fastq_out_unaligned_2,quality,len_tag,dict_tag,seq_len_loop)
                os.system('cat '+sam_out + ' ' + file_out_aligned_loop + ' > '+os.path.join(tmp_dir,'total_sam_tmp'))
                shutil.copy(os.path.join(tmp_dir,'total_sam_tmp'),sam_out)
                os.system('rm ' + file_out_aligned_loop)
                os.system('rm ' + os.path.join(tmp_dir,'total_sam_tmp'))
        print 'looping finished'
        dict_tag.clear()
    else:
        print "alignement already performed"


def pcr_amplification_extract(output_folder, experience_name):
    import os
    import string
    from coll_counter_py2_6 import Counter
    print 'pcr amplification detection'
    # concatenate all results

    big_sam = os.path.join(output_folder, experience_name + '_all_paired.sam')
    search_sam_file = os.path.join(output_folder,'*.sam')
    if not(os.path.exists(big_sam)):
        os.system('cat '+search_sam_file+' >'+big_sam)
    # import list of paired reads
    f = open(big_sam,'r')
    d = Counter()
    all_lines = f.readlines()
    for line in all_lines:
        ind= string.joinfields(line.split()[1:],'\t')
        d[ind] +=1
    file_out = os.path.join(output_folder,experience_name+'_all_paired_pcr_free.sam')
    g = open(file_out,'w')
    for clef in d.keys():
        g.write('%s\t%s\n' % (str(d[clef]),clef))
    d.clear()
    return file_out
def map_read_on_fragment(dict_Pos_restriction, pos_read, chr):
    """      return corresponding fragment """
    id_tmp = ((dict_Pos_restriction[chr] - pos_read)<=0)
    tot_len = len(dict_Pos_restriction[chr])
    len_pre_pos = len(id_tmp.nonzero()[0])
    if len_pre_pos == 0:
        id = 1
    elif (len_pre_pos == tot_len):
        id = tot_len + 1
    else:
        id = id_tmp.nonzero()[0][-1] +2
    return id

def paired_reads_2_frag_contacts(tot_len_read,paired_reads_file,len_tag,dict_Pos_restriction,fragments_contacts_file):
    import numpy as np
    start = time.clock()
    i = 0
    contacts = open(paired_reads_file,'r')
    all_contacts = contacts.readlines()
    print "file fragment = ", fragments_contacts_file
    handle = open(fragments_contacts_file,'w')
    step = 0.0

    n_total_step = float(len(all_contacts))
    for ele in all_contacts:
        data = ele.split('\t')
        ####
        if len_tag>0:
            chr_a = data[1]
            tmp_a = int(data[2])
            or_a = data[3]
            seq_len_a = data[5]
            chr_b = data[6]
            tmp_b = int(data[7])
            or_b = data[8]
            seq_len_b = data[10]
        else:
            chr_a = data[1]
            tmp_a = int(data[2])
            or_a = data[3]
            seq_len_a = data[4]
            chr_b = data[5]
            tmp_b = int(data[6])
            or_b = data[7]
            seq_len_b = data[8]
        ####
        if or_a == 'w':
            pos_a = tmp_a
        else:
            pos_a = tmp_a - int(seq_len_a)
        if or_b == 'w':
            pos_b = tmp_b
        else:
            pos_b = tmp_b - int(seq_len_b)
        ####
        id_a = map_read_on_fragment(dict_Pos_restriction,pos_a,chr_a)
        id_b = map_read_on_fragment(dict_Pos_restriction,pos_b,chr_b)

#        ### closest restriction site read a###
#        if or_a == 'w':
#            list_dist_read_a = ((dict_Pos_restriction[chr_a] - pos_a))
#            dist_clos_rest_a =  list_dist_read_a[list_dist_read_a>=0][0]
#        else:
#            list_dist_read_a = ((dict_Pos_restriction[chr_a] - pos_a))
#            tmp_list_read_a = list_dist_read_a[list_dist_read_a>=0]
#            if tmp_list_read_a == []:
#                dist_clos_rest_a =  abs(list_dist_read_a[list_dist_read_a<=0][0])
#            else:
#                dist_clos_rest_a =  abs(list_dist_read_a[list_dist_read_a>=0][0])
#         ### closest restriction site read b ###
#        if or_b == 'w':
#            list_dist_read_b = ((dict_Pos_restriction[chr_b] - pos_b))
#            dist_clos_rest_b =  list_dist_read_b[list_dist_read_b>=0][0]
#        else:
#            list_dist_read_b = ((dict_Pos_restriction[chr_b] - pos_b))
#            tmp_list_read_b = list_dist_read_b[list_dist_read_b>=0]
#            if tmp_list_read_b == []:
#                dist_clos_rest_b =  abs(list_dist_read_b[list_dist_read_b<=0][0])
#            else:
#                dist_clos_rest_b =  abs(list_dist_read_b[list_dist_read_b>=0][0])

                
                
        list_dist_read_a = abs((dict_Pos_restriction[chr_a] - pos_a))
        list_dist_read_b = abs((dict_Pos_restriction[chr_b] - pos_b))
        if (len(list_dist_read_a) == 0) or (len(list_dist_read_b) == 0) :
            print " alignment on contigs where no restriction site has been detected!!"
        else:
            dist_clos_rest_a =  list_dist_read_a.min()
            dist_clos_rest_b =  list_dist_read_b.min()
            dist_frag_send_2_aligner = dist_clos_rest_a + dist_clos_rest_b

    #        if dist_frag_send_2_aligner<=tot_len_read:
            if 1<=2:
                handle.write("%s\t%s\t%s\t%s\t%s\n" % (str(id_a),chr_a,str(id_b),chr_b,str(dist_frag_send_2_aligner)))
                i = i+1
                # print i

        pt = step * 100 / n_total_step
        print pt, " percent complete\r",
        step += 1

    elapsed = (time.clock() - start)
    print "fragments contacts generated in" + str(elapsed) +' sec'


def load_fragments_contacts_2_db(experience_name,paired_reads):
    input_a = open(paired_reads,'r')
    start = time.clock()
    i= 0
    connection = pymongo.Connection('localhost',27017)
    db = connection[experience_name]
    db.fragments_contacts.drop()
    while 1:
        line_a = input_a.readline()
        frag_contact = line_a.split()
        if not line_a:
            input_a.close()
            elapsed = (time.clock() - start)
            print 'Fragments contacts loaded in data base in ' + str(elapsed)+ ' s'
            break

        i = i+1
        contact = {"_id": i,
                   "frag_a" :frag_contact[0],
                   "contig_a": frag_contact[1],
                   "frag_b":frag_contact[2],
                   "contig_b":frag_contact[3]}
        db.fragments_contacts.insert(contact)
    connection.close()

def fragments_2_db(experience_name, restriction_fragments):
    input_a = open(restriction_fragments,'r')
    start = time.clock()
    i= 0
    connection = pymongo.Connection('localhost',27017)
    db = connection[experience_name]
    db.fragments.drop()
    while 1:
        line_a = input_a.readline()
        data = line_a.split('\t')
        if not line_a:
            input_a.close()
            elapsed = (time.clock() - start)
            print 'fragments loaded in data base in ' + str(elapsed)+ ' s'
            break
        if not(data[0]=='id'):
            print data[0]
            fragment = {"_id":data[0]+','+data[1],
                        "contig":data[1],
                        "start":int(data[2]),
                        "end":int(data[3]),
                        "size":int(data[4]),
                        "gc_content":int(data[5])}
            db.fragments.save(fragment)
        i = i+1
    connection.close()

def remove_self_fragments_contacts(input_file,output_file):
    a = open(input_file,'r')
    b = open(output_file,'w')
    i = 0
    for line in a:
        dat = line.split()
        if not((dat[0] == dat[2]) and (dat[1] == dat[3])):
            i = i+1
            b.write(line)
    print str(i) + ' hetero fragments contacts written'


def rel_frag_2_abs_frag(list_contigs, fragments_contacts_file, fragments_contacts_file_absolute, fragments_list):
    frag_list = open(fragments_list,'r')
    dict_fragments = dict()
    input_contacts = open(fragments_contacts_file,'r')
    ouput_contacts = open(fragments_contacts_file_absolute,'w')
    i = 0
    for ele in list_contigs:
        print 'contig: ' + ele + ' index = ' + str(i)
        dict_fragments[ele] = []
        i = i+1
    i = 0
    for line in frag_list:
        if i >0:
            a = line.split()
            dict_fragments[a[1]].append(a[0])
        i = i+1
    dict_adapt = dict()
    dict_adapt[list_contigs[0]] = 0
    n_contigs = len(list_contigs)
    for i in range(1,n_contigs):
        chunk = dict_adapt[list_contigs[i-1]]+len(dict_fragments[list_contigs[i-1]])
        print 'length '+ list_contigs[i] +' : ' + str(chunk)
        dict_adapt[list_contigs[i]] = chunk
    print 'converting relative contacts...'
    for line in input_contacts:
        line_split = line.split()
        print line_split
        frag_a = dict_adapt[line_split[1]]+int(line_split[0])
        frag_b = dict_adapt[line_split[3]]+int(line_split[2])
        ouput_contacts.write( "%s\t%s\n" % ( str(frag_a),str(frag_b) ) )
    return dict_adapt
