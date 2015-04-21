# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 15:29:39 2012

@author: hervemn
"""
import os
import csv
import sys
import getopt
import numpy
import string
import time
import pymongo

def extract_tag(in_a,in_b,l_tag,experience_name):
    connection = pymongo.Connection('localhost',27017)
    db = connection[experience_name]
    collection = db.reads_tags
    collection.drop()
    collection = db.reads_tags
    input_a = open(in_a,'r')
    input_b = open(in_b,'r')
    start = time.clock()
    i = 0
    while 1:
        line_a = input_a.readline()
        line_b = input_b.readline()
        if (i%4==0):
            title_a = line_a[1:-1]
        elif (i%4==1):
            tag_a = line_a[:l_tag-1]
            tag_b = line_b[:l_tag-1]
            read = {"_id" : title_a,
                    "tag_a": tag_a,
                    "tag_b": tag_b}
            collection.insert(read,{"upsert":"true"})
        i = i+1
        if not line_a:
            elapsed = (time.clock() - start)
            print 'data base loaded in '+str(elapsed)+' s'
            break


def bank_split_by_tag(in_a,in_b,tag_bank_1,tag_bank_2,len_random_tag,name_garbage,dir_bank_1,dir_bank_2,dir_garbage):
    import string
    start = time.clock()
    len_tag_bank = len(tag_bank_1)
    rev_tag_bank_1 =tag_bank_1.translate(string.maketrans('TAGCtagc', 'ATCGATCG'))
    if not(os.path.exists(dir_bank_1)):
        os.mkdir(dir_bank_1)
    if not(os.path.exists(dir_garbage)):
        os.mkdir(dir_garbage)
    if len(tag_bank_2)>0:
        rev_tag_bank_2 =tag_bank_2.translate(string.maketrans('TAGCtagc', 'ATCGATCG'))
        if not(os.path.exists(dir_bank_2)):
            os.mkdir(dir_bank_2)
        bank_b_1= open(os.path.join(dir_bank_2,tag_bank_2+'_'+os.path.basename(in_a)), "w")
        bank_b_2= open(os.path.join(dir_bank_2,tag_bank_2+'_'+os.path.basename(in_b)), "w")
    bank_a_1= open(os.path.join(dir_bank_1,tag_bank_1+'_'+os.path.basename(in_a)), "w")
    bank_a_2= open(os.path.join(dir_bank_1,tag_bank_1+'_'+os.path.basename(in_b)), "w")
    bank_garbage_1= open(os.path.join(dir_garbage,name_garbage+'_'+os.path.basename(in_a)), "w")
    bank_garbage_2= open(os.path.join(dir_garbage,name_garbage+'_'+os.path.basename(in_b)), "w")
    input_a = open(in_a,'r')
    input_b = open(in_b,'r')
    i = 0
    while 1:
        line_a = input_a.readline()
        line_b = input_b.readline()
        if not line_a and not line_b:
            bank_a_1.close()
            bank_a_2.close()
            bank_garbage_1.close()
            bank_garbage_2.close()
            if len(tag_bank_2)>0:
                bank_b_1.close()
                bank_b_2.close()
            break
            elapsed = (time.clock() - start)
            print 'bank splited in ' + str(elapsed)+ ' s'
        if (i%4==0):
            title_a = line_a[1:]
            title_b = line_b[1:]
        elif (i%4==1):
            seq_a = line_a
            seq_b = line_b
        elif (i%4==3):
            qual_a = line_a
            qual_b = line_b
            random_tag_a = seq_a[0:len_random_tag]
            random_tag_b = seq_b[0:len_random_tag]
            tag_a = seq_a[len_random_tag:len_random_tag+len_tag_bank]
            tag_b = seq_b[len_random_tag:len_random_tag+len_tag_bank]
            block_a = seq_a[len_random_tag+len_tag_bank:]
            block_b = seq_b[len_random_tag+len_tag_bank:]
            block_a_qual = qual_a[:len_random_tag]+qual_a[len_random_tag+len_tag_bank:]
            block_b_qual = qual_b[:len_random_tag]+qual_b[len_random_tag+len_tag_bank:]
            garbage_block_a = seq_a[len_random_tag:]
            garbage_block_b = seq_b[len_random_tag:]
            qual_garbage_a = qual_a[len_random_tag:]
            qual_garbage_b = qual_b[len_random_tag:]
            #            print "tag read 1"+tag_a
            #            print "tag read 2"+tag_b
            if ((tag_a[0:2] == tag_bank_1[0:2]) or (tag_a[0:2] == rev_tag_bank_1[0:2])) and ((tag_b[0:2] == tag_bank_1[0:2]) or (tag_b[0:2] == rev_tag_bank_1[0:2])):
                bank_a_1.write("@%s%s+%s%s" % (title_a, random_tag_a+block_a,title_a,
                                               block_a_qual))
                bank_a_2.write("@%s%s+%s%s" % (title_b, random_tag_b+block_b,title_b,
                                               block_b_qual))
            elif (len(tag_bank_2)>0) and ((tag_a[0:2] == tag_bank_2[0:2]) or (tag_a[0:2] == rev_tag_bank_2[0:2])) and\
                 ((tag_b[0:2] == tag_bank_2[0:2]) or (tag_b[0:2] == rev_tag_bank_2[0:2])):
                bank_b_1.write("@%s%s+%s%s" % (title_a, random_tag_a+block_a,title_a,
                                               block_a_qual))
                bank_b_2.write("@%s%s+%s%s" % (title_b, random_tag_b+block_b,title_b,
                                               block_b_qual))
            else:
                bank_garbage_1.write("@%s%s+%s%s" % (title_a, garbage_block_a,title_a,
                                                     qual_garbage_a))
                bank_garbage_2.write("@%s%s+%s%s" % (title_b,garbage_block_b,title_b,
                                                     qual_garbage_b))
        i = i+1


def bank_split_by_tag_1(in_a,in_b,dir_bank_1,dir_bank_2,motif1):
    import string
    rev_motif1 = motif1.translate(string.maketrans('TAGCtagc', 'ATCGATCG'))
    if not(os.path.exists(dir_bank_1)):
        os.mkdir(dir_bank_1)
    if not(os.path.exists(dir_bank_2)):
        os.mkdir(dir_bank_2)
    trim =len(motif1)
    handle_a_1= open(os.path.join(dir_bank_1,motif1+'_'+os.path.basename(in_a)), "w")
    handle_a_2= open(os.path.join(dir_bank_2,os.path.basename(in_a)), "w")
    handle_b_1= open(os.path.join(dir_bank_1,motif1+'_'+os.path.basename(in_b)), "w")
    handle_b_2= open(os.path.join(dir_bank_2,os.path.basename(in_b)), "w")
    i = 0;
    input_a = open(in_a,'r')
    input_b = open(in_b,'r')
    start = time.clock()
    while 1:
        line_a = input_a.readline()
        line_b = input_b.readline()
        if (i%4==0):
            title_a = line_a[1:]
            title_b = line_b[1:]
        elif (i%4==1):
            seq_a = line_a
            seq_b = line_b
        elif (i%4==3):
            qual_a = line_a
            qual_b = line_b
            if ((seq_a[:trim] == motif1) or (seq_a[:trim] == rev_motif1)) and ((seq_b[:trim] == motif1) or (seq_b[:trim] == rev_motif1)):
                handle_a_1.write("@%s%s+%s%s" % (title_a, seq_a[trim:],title_a,
                                                 qual_a[trim:]))
                handle_b_1.write("@%s%s+%s%s" % (title_b, seq_b[trim:],title_b,
                                                 qual_b[trim:]))
            else:
                handle_a_2.write("@%s%s+%s%s" % (title_a, seq_a,title_a,
                                                 qual_a))
                handle_b_2.write("@%s%s+%s%s" % (title_b, seq_b,title_b,
                                                 qual_b))
        i = i+1
        if not line_a and not line_b:
            handle_a_1.close()
            handle_b_1.close()
            handle_a_2.close()
            handle_b_2.close()
            break
            elapsed = (time.clock() - start)
            print 'bank splited in ' + str(elapsed)+ ' s'


def remove_bank_tag(in_a,in_b,dir_bank_1,start_random_tag,end_random_tag,start_tag,end_tag):
    if not(os.path.exists(dir_bank_1)):
        os.mkdir(dir_bank_1)
    motif1 = 'trimmed'
    handle_a_1= open(os.path.join(dir_bank_1,motif1+'_'+os.path.basename(in_a)), "w")
    handle_b_1= open(os.path.join(dir_bank_1,motif1+'_'+os.path.basename(in_b)), "w")
    i = 0
    input_a = open(in_a,'r')
    input_b = open(in_b,'r')
    start = time.clock()
    while 1:
        line_a = input_a.readline()
        line_b = input_b.readline()
        if (i%4==0):
            title_a = line_a[1:]
            title_b = line_b[1:]
        elif (i%4==1):
            seq_a = line_a
            seq_b = line_b
        elif (i%4==3):
            qual_a = line_a
            qual_b = line_b
            random_tag_a = seq_a[start_random_tag-1:end_random_tag]
            random_tag_b = seq_b[start_random_tag-1:end_random_tag]
            tag_a = seq_a[start_tag-1:end_tag]
            tag_b = seq_b[start_tag-1:end_tag]
            block_a = seq_a[end_tag:]
            block_b = seq_b[end_tag:]
            block_a_qual = qual_a[start_random_tag-1:end_random_tag]+qual_a[end_tag:]
            block_b_qual = qual_b[start_random_tag-1:end_random_tag]+qual_b[end_tag:]
            handle_a_1.write("@%s%s+%s%s" % (title_a, random_tag_a+block_a,title_a,
                                             block_a_qual))
            handle_b_1.write("@%s%s+%s%s" % (title_b, random_tag_b+block_b,title_b,
                                             block_b_qual))
        i = i+1
        if not line_a and not line_b:
            handle_a_1.close()
            handle_b_1.close()
            break
            elapsed = (time.clock() - start)
            print 'bank trimmed in ' + str(elapsed)+ ' s'



def bank_split_by_tag_gene(in_a,in_b,info_bank,len_random_tag,info_garbage):
    import string
    import numpy as np
    import gzip
    list_bank = info_bank.keys()
    n_bank = len(list_bank)

    for ele in list_bank:
        tag = info_bank[ele]['tag']
        len_tag_bank = len(tag)
        rev_tag =tag.translate(string.maketrans('TAGCtagc', 'ATCGATCG'))
        info_bank[ele]['rev_tag'] = rev_tag
        dir_bank = info_bank[ele]['dir_bank']
        if not(os.path.exists(dir_bank)):
            os.mkdir(dir_bank)
        info_bank[ele]['handle_1'] = open(os.path.join(dir_bank,tag+'_'+os.path.basename(in_a)+'bet'), "w")
        info_bank[ele]['handle_2'] = open(os.path.join(dir_bank,tag+'_'+os.path.basename(in_b)+'bet'), "w")

    dir_garbage = info_garbage['dir_bank']
    if not(os.path.exists(dir_garbage)):
        os.mkdir(dir_garbage)
    info_garbage['handle_1'] = open(os.path.join(dir_garbage,'garbage_'+os.path.basename(in_a)), "w")
    info_garbage['handle_2'] = open(os.path.join(dir_garbage,'garbage_'+os.path.basename(in_b)), "w")

    ext_a = in_a.split('.')[-1]
    ext_b = in_b.split('.')[-1]
    if ((ext_a == 'gz') or (ext_b == 'gz')):
        input_a = gzip.open(in_a,'r')
        input_b = gzip.open(in_b,'r')
    else:
        input_a = open(in_a,'r')
        input_b = open(in_b,'r')
    i = 0
    start = time.clock()
    while 1:
        line_a = input_a.readline()
        line_b = input_b.readline()
        if not line_a and not line_b:
            for ele in list_bank:
                info_bank[ele]['handle_1'].close()
                info_bank[ele]['handle_2'].close()
            info_garbage['handle_1'].close()
            info_garbage['handle_2'].close()
            elapsed = (time.clock() - start)
            print 'bank splited in ' + str(elapsed)+ ' s'
            break
        if (i%4==0):
            title_a = line_a[1:]
            title_b = line_b[1:]
        elif (i%4==1):
            seq_a = line_a
            seq_b = line_b
        elif (i%4==3):
            qual_a = line_a
            qual_b = line_b
            random_tag_a = seq_a[0:len_random_tag]
            random_tag_b = seq_b[0:len_random_tag]
            tag_a = seq_a[len_random_tag:len_random_tag+len_tag_bank]
            tag_b = seq_b[len_random_tag:len_random_tag+len_tag_bank]
            block_a = seq_a[len_random_tag+len_tag_bank:]
            block_b = seq_b[len_random_tag+len_tag_bank:]
            block_a_qual = qual_a[:len_random_tag]+qual_a[len_random_tag+len_tag_bank:]
            block_b_qual = qual_b[:len_random_tag]+qual_b[len_random_tag+len_tag_bank:]
#            garbage_block_a = seq_a[len_random_tag:]
#            garbage_block_b = seq_b[len_random_tag:]
#            qual_garbage_a = qual_a[len_random_tag:]
#            qual_garbage_b = qual_b[len_random_tag:]

            res_test = np.zeros(n_bank)
            for k in range(0,n_bank):
                tmp_t1 = len({l for l, (left, right) in enumerate(zip(tag_a,info_bank[list_bank[k]]['tag'])) if left == right})
                tmp_t1p = len({l for l, (left, right) in enumerate(zip(tag_a,info_bank[list_bank[k]]['rev_tag'])) if left == right})

                tmp_t2 = len({l for l, (left, right) in enumerate(zip(tag_b,info_bank[list_bank[k]]['tag'])) if left == right})
                tmp_t2p = len({l for l, (left, right) in enumerate(zip(tag_b,info_bank[list_bank[k]]['rev_tag'])) if left == right})

                t1 = np.max([tmp_t1,tmp_t1p])
                t2 = np.max([tmp_t2,tmp_t2p])

#                score = t1 + t2
                score = np.max([t1,t2])
#                res_test[k] = score/(len_tag_bank *2.)
                res_test[k] = score/(len_tag_bank)
            if res_test.max()>=0.65:
                """ write corresponding bank"""
                index_bank = res_test.argmax()
                name_bank = list_bank[index_bank]
                info_bank[name_bank]['handle_1'].write("@%s%s+%s%s" % (title_a, random_tag_a+block_a,title_a,
                                                                       block_a_qual))
                info_bank[name_bank]['handle_2'].write("@%s%s+%s%s" % (title_b, random_tag_b+block_b,title_b,
                                                                       block_b_qual))
            else:
                """ write garbage"""
                info_garbage['handle_1'].write("@%s%s+%s%s" % (title_a, random_tag_a+block_a,title_a,
                                                               block_a_qual))
                info_garbage['handle_2'].write("@%s%s+%s%s" % (title_b, random_tag_b+block_b,title_b,
                                                               block_b_qual))
        i = i+1
