__author__ = 'hervemarie-nelly'
import os
import numpy
import pymongo
import shutil
import re
import time

class hic:
    def __init__(self,experience_name, restriction_sequence,reads_a,reads_b,data_output):
        self.experience_name = experience_name
        self.restriction_sequence = restriction_sequence
        self.reads_a = reads_a
        self.reads_b = reads_b
        self.output_folder = data_output
        if not(os.path.exists(data_output + 'tmp_files')):
            os.mkdir(data_output + 'tmp_files')
        self.tmp = data_output + 'tmp_files/'
        self.sam_out = data_output + 'test.sam'
        self.out_unaligned = data_output + 'unaligned.fastq'
        self.dict_Pos_restriction = dict()
        self.contacts_fragments = self.tmp + 'contacts_fragments.txt'

    def align_fastq(self,bowtie2,gen_index,n_cpu, looping, quality):
        parameters = ' --sam-no-hd --sam-no-sq --very-sensitive --minins 0 --maxins 5 -M 5 -p' + str(n_cpu)
        genome_index = ' -x '+ gen_index
        reads_a = ' -1 ' + self.reads_a
        reads_b = ' -2 ' + self.reads_b
        file_out = self.tmp+'tmp_out.sam'
        handle = open(self.reads_a,'r')
        handle.next()
        seq_len = len(handle.next())
        handle.close()
        os.system(bowtie2 + parameters + genome_index + reads_a + reads_b + ' -S '+ file_out)
        ############ extract the paired reads that did not align from the list #########
        file_out_aligned = self.tmp+'tmp_out_aligned.sam'
        fastq_out_unaligned_1 =self.tmp+'tmp_1_fastq.txt'
        fastq_out_unaligned_2 =self.tmp+'tmp_2_fastq.txt'
        sam_filter(file_out,file_out_aligned,fastq_out_unaligned_1,fastq_out_unaligned_2,quality)
        ############ writing output file or startin looping ############################
        shutil.copy(file_out_aligned,self.sam_out)
        if not looping:
            print 'Alignment done...'
        else:
            print 'Start looping...'
            for trim_3 in range(1,seq_len - 20 ):
                reads_a = ' -1 ' + fastq_out_unaligned_1
                reads_b = ' -2 ' + fastq_out_unaligned_2
                file_out = self.tmp+'tmp_out_loop.sam'
                os.system(bowtie2 + parameters + genome_index + reads_a + reads_b +' -S '+ file_out)
                file_out_aligned_loop = self.tmp+'tmp_out_loop_aligned.sam'
                sam_filter(file_out,file_out_aligned_loop,fastq_out_unaligned_1,fastq_out_unaligned_2,quality)
                os.system('cat '+self.sam_out + ' ' + file_out_aligned_loop + ' > '+self.tmp+'total_sam_tmp')
                shutil.copy(self.tmp+'total_sam_tmp',self.sam_out)
                os.system('rm ' + file_out_aligned_loop)
                os.system('rm ' + self.tmp+'total_sam_tmp')
        self.seq_len = seq_len

    def pre_build_contact(self,use_tag,):
        input_a = open(self.sam_out,'r')
        start = time.clock()
        self.paired = self.output_folder+'contacts.paired'
        handle_aligned = open(self.paired,'w')
        i= 0
        if use_tag:
            connection = pymongo.Connection('localhost',27017)
            db = connection[self.experience_name]
        while 1:
            line_a = input_a.readline()
            if not line_a:
                input_a.close()
                handle_aligned.close()
                elapsed = (time.clock() - start)
                print 'Extract aligned paired reads in ' + str(elapsed)+ ' s'
                break
            if (i%2==0):
                read_a = line_a
            elif (i%2==1):
                read_b = line_a
                tmp_a = read_a.split('\t')
                tmp_b = read_b.split('\t')
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
                if use_tag:
                    data_paired_read = db.reads_tags.find_one({"_id":tmp_a[0]})
                    tag_a = data_paired_read['tag_a']
                    tag_b = data_paired_read['tag_b']
                    handle_aligned.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (tmp_a[2],tmp_a[3],or_read_a,tag_a,
                                                                               tmp_b[2],tmp_b[3],or_read_b,tag_b))
                else:
                    handle_aligned.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (tmp_a[2],tmp_a[3],or_read_a,
                                                                       tmp_b[2],tmp_b[3],or_read_b))
            i = i+1

    def load_data_2_mongo(self,):
        input_a = open(self.paired,'r')
        start = time.clock()
        i= 0
        connection = pymongo.Connection('localhost',27017)
        db = connection[self.experience_name]
        while 1:
            line_a = input_a.readline()
            if not line_a:
                input_a.close()
                elapsed = (time.clock() - start)
                print 'Paired reads loaded in data base in ' + str(elapsed)+ ' s'
                break
            if db.contacts.find_one({"_id":line_a[:-1]}) is None:
                db.contacts.save({"_id":line_a[:-1], "count":1})
            else:
                db.contacts.update({"_id":line_a[:-1]},{"$inc":{ "count":1} })
            i = i+1
    def restriction_map(self, genome_fasta):
        input_a = open(genome_fasta,'r')
        start = time.clock()
        self.restriction_site_list = self.output_folder+restriction_site_list
        self.fragments_list = self.output_folder+fragments_list
        handle_site_list = open(self.restriction_site_list,'w')
        handle_fragments_list = open(self.fragments_list,'w')
        i = 0
        D = dict()
        while 1:
            line_a = input_a.readline()
            if not line_a:
                input_a.close()
                break
            if line_a[0] == '>':
                tmp_key = line_a[1:-1]
                D[tmp_key] = ''
                self.dict_Pos_restriction[tmp_key] = []
            else:
                D[tmp_key] = D[tmp_key]+ line_a[:-1]
            i = i+1
        handle_fragments_list.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ('id','chrom',
                                                                  'start_pos','end_pos','size','gc_content'))
        for chrom in D.keys():
            print restriction_site
            starts = [match.start() for match in re.finditer(re.escape(self.restriction_sequence), D[chrom])]
            self.dict_Pos_restriction[chrom] = numpy.array(starts)
            pos_init = 0
            id_frag = 1
            for pos in starts:
                handle_site_list.write("%s\t%s\n" % (chrom,str(pos)))
                size_frag = pos - pos_init
                gc_content = D[chrom][pos_init:pos].count('G')+D[chrom][pos_init:pos].count('C')
                handle_fragments_list.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (str(id_frag),chrom,
                                                                          str(pos_init),str(pos),str(size_frag),
                                                                          str(gc_content)))
                pos_init = pos
                id_frag = id_frag +1

            handle_site_list.write("%s\t%s\n" % (chrom,str(pos)))
            pos = len(D[chrom])
            size_frag = pos - pos_init
            gc_content = D[chrom][pos_init:pos].count('G')+D[chrom][pos_init:pos].count('C')
            handle_fragments_list.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (str(id_frag),chrom,
                                                                      str(pos_init),str(pos),str(size_frag),
                                                                      str(gc_content)))
        handle_site_list.close()
        handle_fragments_list.close()
        elapsed = (time.clock() - start)
        print 'Restriction map generated in ' + str(elapsed)+ ' s'

    def fragments_2_db(self,):
        input_a = open(self.fragments_list,'r')
        start = time.clock()
        i= 0
        connection = pymongo.Connection('localhost',27017)
        db = connection[self.experience_name]
        db.fragments.drop()
        while 1:
            line_a = input_a.readline()
            data = line_a.split('\t')
            if not line_a:
                input_a.close()
                elapsed = (time.clock() - start)
                print 'Paired reads loaded in data base in ' + str(elapsed)+ ' s'
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

    def fragments_contacts_2_db(self,use_tag,contact_file):
        start = time.clock()
        i = 0
        handle = open(contact_file,'w')
        connection = pymongo.Connection('localhost',27017)
        db = connection[self.experience_name]
        db.fragments_contacts.drop()
        for ele in db.contacts.find():
            data = ele["_id"].split('\t')
            if use_tag:
                chr_a = data[0]
                tmp_a = int(data[1])
                or_a = data[2]
                chr_b = data[4]
                tmp_b = int(data[5])
                or_b = data[6]
            else:
                chr_a = data[0]
                tmp_a = int(data[1])
                or_a = data[2]
                chr_b = data[3]
                tmp_b = int(data[4])
                or_b = data[5]
            if or_a == 'w':
                pos_a = tmp_a
            else:
                pos_a = tmp_a - self.seq_len
            if or_b == 'w':
                pos_b = tmp_b
            else:
                pos_b = tmp_b - self.seq_len

            id_a_tmp = ((self.dict_Pos_restriction[chr_a] - pos_a)<=0)
            if len(id_a_tmp.nonzero()[0]) == 0:
                id_a = 1
            else:
                id_a = id_a_tmp.nonzero()[0][-1] +2

            id_b_tmp = ((self.dict_Pos_restriction[chr_b] - pos_b)<=0)
            if len(id_b_tmp.nonzero()[0]) == 0:
                id_b = 1
            else:
                id_b = id_b_tmp.nonzero()[0][-1] +2

            contact = {"_id":i,
                       "fragment_a":str(id_a)+','+chr_a,
                       "fragment_b":str(id_b)+','+chr_b}
            db.fragments_contacts.save(contact)
            handle.write("%s\t%s\t%s\t%s\n" % (str(id_a),chr_a,str(id_b),chr_b))
            i = i+1
        elapsed = (time.clock() - start)
        print "fragments contacts generated in" + str(elapsed) +' sec'




def sam_filter(sam,sam_aligned,fastq_unaligned_1,fastq_unaligned_2,quality):
    input_a = open(sam,'r')
    start = time.clock()
    handle_aligned = open(sam_aligned,'w')
    handle_unaligned_1 = open(fastq_unaligned_1,'w')
    handle_unaligned_2 = open(fastq_unaligned_2,'w')
    i= 0
    while 1:
        line_a = input_a.readline()
        if not line_a:
            handle_aligned.close()
            handle_unaligned_1.close()
            handle_unaligned_2.close()
            elapsed = (time.clock() - start)
            print 'Extract aligned paired reads in ' + str(elapsed)+ ' s'
            break
        if i%2==0:
            read_a = line_a
        elif i%2==1:
            read_b = line_a
            tmp_a = read_a.split('\t')
            tmp_b = read_b.split('\t')
            if len(tmp_a)<=11:
                print read_a
                print read_b
            if ((tmp_a[2] == '*' ) or (tmp_b[2] == '*') or (int(tmp_a[4])<quality)
                or (int(tmp_b[4])<quality)):
                handle_unaligned_1.write("@%s\n%s\n+%s\n%s\n" % (tmp_a[0],
                                                                 tmp_a[9],tmp_a[0],tmp_a[10]))
                handle_unaligned_2.write("@%s\n%s\n+%s\n%s\n" % (tmp_b[0],
                                                                 tmp_b[9],tmp_b[0],tmp_b[10]))
            else:
                handle_aligned.write("%s%s" % (read_a,read_b))
        i = i+1