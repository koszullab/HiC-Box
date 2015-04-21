# coding: utf-8
# -*- coding: utf-8 -*-_
__author__ = 'hervemn'

import os
import numpy as np
class data_hic():
    def __init__(self,name,folder_analysis):
        print name
        self.name = name
        self.file_info_contigs = os.path.join(folder_analysis,'info_contigs.txt')
        self.file_abs_frag_contact_weighted = os.path.join(folder_analysis,'abs_fragments_contacts_weighted.txt')
        self.file_bias_gc_contact = os.path.join(folder_analysis,'bias_gc_contacts.txt')
        self.file_bias_size_contact = os.path.join(folder_analysis,'bias_size_contacts.txt')
        self.file_bins_gc = os.path.join(folder_analysis,'bins_gc.txt')
        self.file_bins_length = os.path.join(folder_analysis,'bins_length.txt')
        self.file_contacts_vs_genomic_distance = os.path.join(folder_analysis,'contacts_vs_genomic_distance.txt')
        self.file_info_contigs = os.path.join(folder_analysis,'info_contigs.txt')
        self.file_frag_contacts_weighted = os.path.join(folder_analysis,'fragments_contacts_weighted.txt')
        self.file_mat_bias_gc = os.path.join(folder_analysis,'gc_bias_mat_np.txt')
        self.file_length_bias_gc =  os.path.join(folder_analysis,'length_bias_mat_np.txt')
        self.file_weigth_distances =  os.path.join(folder_analysis,'weight_norm_step_2000.txt')

        ######################## Define contigs  ####################################
        info_contigs = open(self.file_info_contigs,'r')
        header_info = info_contigs.readline().split()
        self.dict_contigs = dict()
        while 1:
            line = info_contigs.readline()
            if not(line):
                info_contigs.close()
                break
            data = line.split()
            for i in range(0,len(header_info)):
                if i == 0:
                    self.dict_contigs[data[0]] = dict()
                else:
                    self.dict_contigs[data[0]][header_info[i]] = int(data[i])


        for contig in self.dict_contigs.keys():
            n_frags = self.dict_contigs[contig]['n_frags']
            if n_frags ==1:
                self.dict_contigs[contig]['ext_left'] = 1
                self.dict_contigs[contig]['ext_right'] = 1
            else:
                self.dict_contigs[contig]['ext_left'] = 1+(n_frags/2 -1)
                self.dict_contigs[contig]['ext_right'] = n_frags - (n_frags/2 -1)
        ###################### define numpy contacts ################################
        handle_abs_frag_contacts  = open(self.file_abs_frag_contact_weighted,'r')
        handle_abs_frag_contacts.readline()
        contacts = []
        while 1:
            line = handle_abs_frag_contacts.readline()
            if not(line):
                handle_abs_frag_contacts.close()
                break
            data = line.split()
            contacts.append([int(data[0]),int(data[1]),float(data[2]),float(data[3]) ] )
        self.contacts_np = np.array(contacts)

    def polarize(self):
        print 'start polarizing'

        for contig in self.dict_contigs.keys():

            out_votes_index = np.nonzero(self.contacts_np[:,0])

            if self.dict_contigs[contig]['n_frags'] == 1:

            else:

