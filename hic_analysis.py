# coding: utf-8
# -*- coding: utf-8 -*-_
__author__ = 'hervemn'

def binomialCoeff(n, k):
    result = 1
    for i in range(1, k+1):
        result = result * (n-i+1) / i
    return result


def draw_matrix(output_folder,fragments_contacts_file_absolute,dict_fragments):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    input_contacts = open(fragments_contacts_file_absolute,'r')
    frag_a_list = []
    frag_b_list = []
    weight_length_list = []
    weight_gc_list = []
    print 'load contacts fragments...'
    n_fragments = len(dict_fragments)
    line = input_contacts.readline()
    while 1:
        line = input_contacts.readline()
        if not line:
            input_contacts.close()
            break
        line_split = line.split()
        if not( line_split[0] == line_split[1]):
            frag_a_list.append(int(line_split[0]))
            frag_b_list.append(int(line_split[1]))
            weight_length_list.append(float(line_split[2]))
            weight_gc_list.append(float(line_split[3]))

    frag_a_np = np.array(frag_a_list)
    frag_b_np = np.array(frag_b_list)
    weight_length_np = np.array(weight_length_list)
    print 'drawing histogram'
    print n_fragments
    H, xedges, yedges = np.histogram2d(frag_a_np, frag_b_np,weights=weight_length_np,
        bins=( range(0,n_fragments+1), range(0,n_fragments+1) ) )
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    plt.imshow(np.log2(H), extent=extent, interpolation='nearest', cmap = 'gist_ncar')
    file_mat_unbiased = os.path.join(output_folder,'np_mat_unbiased.txt')
    print 'Shape matrix = ' + str(H.shape)
    np.savetxt(file_mat_unbiased,H,delimiter='\t')
    file_graph = os.path.join(output_folder,'unbiased_matrix.tif')
    plt.savefig(file_graph,dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1)
    plt.colorbar()
    plt.show()
    plt.close()

def contact_vs_gen_distance(output_folder,dict_fragments,fragments_contacts,dict_contigs):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    fragments_contacts_from_file = open(fragments_contacts,'r')
    contact_vs_distance = []
    out_file = os.path.join(output_folder,'contacts_vs_genomic_distance.txt')
    ## computing normalizing factor ###
    norm_factor_dist = []
    step_graph = 2000
    out_norm_weight_file = os.path.join(output_folder,'weight_norm_step_'+str(step_graph)+'.txt')
    out_norm_weight = open(out_norm_weight_file,'w')
    length_contigs = []
    for chr in dict_contigs.keys():
#        print chr +' '+ str(dict_contigs[chr])
        length_contigs.append(dict_contigs[chr]['length(kb)'])
    length_contigs = np.array(length_contigs)
    length_max_contigs = length_contigs.max()
    print 'length max = ' + str(length_contigs.max())

    bins_graph = range(step_graph,length_max_contigs+step_graph,step_graph)
    out_norm_weight.write("%s\t%s\n" %('distance','weight'))
    for i in bins_graph:
        tmp = 0
        for contigs in dict_contigs.keys():
            tmp0 = (dict_contigs[contigs] - i)*np.sqrt(2)
            if tmp0>=0:
                tmp = tmp0 + tmp
        out_norm_weight.write("%s\t%s\n" %(str(i),str(tmp)))
        norm_factor_dist.append(tmp)
    ## computing normalizing factor  done ###
    if not(os.path.exists(out_file)):
        out_contact_vs_distance = open(out_file,'w')
        for line in fragments_contacts_from_file:
            data = line.split()
            if (data[1] == data[3]) and (data[0] != data[2]):
                chr = data[1]
                id_a = data[0]+'-'+data[1]
                id_b = data[2]+'-'+data[3]
                start_a = int(dict_fragments[id_a]['start'])
                start_b = int(dict_fragments[id_b]['start'])
                dist = abs(start_a - start_b)
                contact_vs_distance.append(dist)
                out_contact_vs_distance.write("%s\t%s\n" % (str(dist),chr))
    else:
        out_contact_vs_distance = open(out_file,'r')
        for line in out_contact_vs_distance:
            data = line.split()
            contact_vs_distance.append(int(data[0]))

    # the histogram of the data

    hist,bin_edges = np.histogram(contact_vs_distance, bins=bins_graph)
    norm_data = np.array(hist)/np.array(norm_factor_dist[1:])
    plt.loglog(bins_graph[1:],norm_data)
#    print hist
    plt.xlabel('Genomic distance')
    plt.ylabel('Contact frequency')
    plt.title(r'$\mathrm{Contact frequency\ vs\ genomic distance}$')
    #plt.axis([40, 160, 0, 0.03])
    plt.grid(True)
    file_graph = os.path.join(output_folder,'contact_vs_gen_dist.png')
    plt.savefig(file_graph,dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1)

    plt.close()
    print 'freq contact vs genomic distance done'



def gc_size_bias(output_folder,dict_fragments,fragments_contacts):
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    file_numpy_mat_gc = os.path.join(output_folder,'gc_bias_mat_np.txt')
    file_numpy_mat_length = os.path.join(output_folder,'length_bias_mat_np.txt')
    file_numpy_steps_gc = os.path.join(output_folder,'bins_gc.txt')
    file_numpy_steps_length = os.path.join(output_folder,'bins_length.txt')
    if not(os.path.exists(file_numpy_mat_gc) & os.path.exists(file_numpy_mat_length)):
        fragments_contacts_from_file = open(fragments_contacts,'r')
        tutti_gc_tmp = []
        data_size_2d = []
        data_gc_2d = []
        if not(os.path.exists(os.path.join(output_folder,'bias_size_contact.txt'))):
            print 'data do not exist yet. computing ...'
            file_size_2d = open(os.path.join(output_folder,'bias_size_contact.txt'),'w')
            file_gc_2d = open(os.path.join(output_folder,'bias_gc_contact.txt'),'w')
            for line in fragments_contacts_from_file:
                data = line.split()
                id_a = data[0]+'-'+data[1]
                id_b = data[2]+'-'+data[3]
                ### CODE DE PORC !!!!! A  MODIFIER!!!!!!!!!!!!!!!!!! ##### voir pourquoi les reads extremes deconnent!!!
                if dict_fragments.has_key(id_a) and dict_fragments.has_key(id_b) :
                    gc_read_a = dict_fragments[id_a]['gc_content']
                    gc_read_b = dict_fragments[id_b]['gc_content']
                    size_read_a = dict_fragments[id_a]['size']
                    size_read_b = dict_fragments[id_b]['size']
                    tutti_gc_tmp.append(float(gc_read_a))
                    tutti_gc_tmp.append(float(gc_read_b))
                    if not(id_a == id_b):
                        file_gc_2d.write("%s\t%s\n" %(gc_read_a,gc_read_b))
                        file_size_2d.write("%s\t%s\n" %(size_read_a,size_read_b))
                        data_size_2d.append([int(size_read_a),int(size_read_b)])
                        data_gc_2d.append([float(gc_read_a),float(gc_read_b)])
                        tutti_gc_tmp.append(float(gc_read_a))
                        tutti_gc_tmp.append(float(gc_read_b))
            file_gc_2d.close()
            file_size_2d.close()
            print 'done'
        else:
            print 'data exist. loading...'
            file_size_2d = open(os.path.join(output_folder,'bias_size_contact.txt'),'r')
            file_gc_2d = open(os.path.join(output_folder,'bias_gc_contact.txt'),'r')
            while 1:
                line_a = file_size_2d.readline()
                line_b = file_gc_2d.readline()
                if not line_a:
                    file_gc_2d.close()
                    file_size_2d.close()
                    print 'done'
                    break
                size_data = line_a.split()
                gc_data = line_b.split()
                data_size_2d.append([int(size_data[0]),int(size_data[1])])
                data_gc_2d.append([float(gc_data[0]),float(gc_data[1])])
                tutti_gc_tmp.append(float(gc_data[0]))
                tutti_gc_tmp.append(float(gc_data[1]))

        numpy_size = np.array(data_size_2d)
        numpy_gc = np.array(data_gc_2d)
        tutti_gc = np.array(tutti_gc_tmp)
        print 'data collect_done'
        print " numpy size shape = ", numpy_size.shape
        print " numpy gc shape = ", numpy_gc.shape
    ########### collecting data #########################################################
        vect_size = []
        vect_gc = []
        for frag in dict_fragments.keys():
            vect_size.append(float(dict_fragments[frag]['size']))
            vect_gc.append(float(dict_fragments[frag]['gc_content']))
        np_vect_gc = np.array(vect_gc)
        np_vect_size = np.array(vect_size)
        print " np vect size shape = ", np_vect_size.shape
        print " np vect gc shape = ", np_vect_gc.shape
    ########### theoretical size matrix ##################################################
        step_size = 500
        size_min = numpy_size.min()
        size_max = numpy_size.max()
        size_bins = range(size_min,size_max+step_size,step_size)
        size_bins = np.array(size_bins)
        theo_array_size = np.zeros((len(size_bins)-1,len(size_bins)-1))
#        for ind_x in range(0,len(size_bins)-1):
#            bin_min_x = size_bins[ind_x]
#            bin_max_x = size_bins[ind_x+1]
#            print "bin min x = ", bin_min_x
#            print "bin max x = ", bin_max_x
#            new_index = np.nonzero((np_vect_size>=bin_min_x) & (np_vect_size<bin_max_x) )[0]
#            print " index = ", new_index
#            print " index shape = ", new_index.shape
#            n_x = np_vect_size[:, new_index]
#            tot_x = len(n_x)
#            for ind_y in range(0,len(size_bins)-1):
#                bin_min_y = size_bins[ind_y]
#                bin_max_y = size_bins[ind_y+1]
#                n_y = np_vect_size[:,np.nonzero( (np_vect_size>=bin_min_y) & (np_vect_size<bin_max_y) )[0]]
#                tot_y = len(n_y)
#                out_mat = tot_x*tot_y
#                theo_array_size[ind_x,ind_y] = out_mat
#
#        theo_array_size[np.nonzero(theo_array_size == 0)[0]] = 1
#        print 'theo size matrix computed...'
    ########### theoretical gc matrix ##################################################
        step_gc = 10000
#        print 'tutti gc'
#        print tutti_gc
        gc_max = np.round(tutti_gc*100 * step_gc).max()
        gc_min = np.floor(tutti_gc*100 * step_gc).min()
        gc_bins = np.array(range(int(gc_min),int(gc_max),step_gc))/float( 100 * step_gc)
        print tutti_gc.min()
        gc_bins[0] = tutti_gc.min()
        gc_bins[-1] = tutti_gc.max()
        theo_array_gc = np.zeros((len(gc_bins)-1,len(gc_bins)-1))
        gc_theo_bins =[]
#        for ind_x in range(0,len(gc_bins)-1):
#            bin_min_x = gc_bins[ind_x]
#            bin_max_x = gc_bins[ind_x+1]
#            n_x = np_vect_gc[:,np.nonzero((np_vect_gc>=bin_min_x) & (np_vect_gc<bin_max_x) )[0]]
#            tot_x = len(n_x)
#            gc_theo_bins.append(tot_x)
#            for ind_y in range(0,len(gc_bins)-1):
#                bin_min_y = gc_bins[ind_y]
#                bin_max_y = gc_bins[ind_y+1]
#                n_y = np_vect_gc[:,np.nonzero( (np_vect_gc>=bin_min_y) & (np_vect_gc<bin_max_y) )[0]]
#                tot_y = len(n_y)
#                out_mat = tot_x*tot_y
#                theo_array_gc[ind_x,ind_y] = out_mat
#        theo_array_gc[np.nonzero(theo_array_gc == 0)[0]] = 1
#        print 'theo gc matrix computed...'


    ########################## size content map ########################################
        print 'drawing histograms'
        H_size, xedges_size, yedges_size = np.histogram2d(numpy_size[:,0], numpy_size[:,1], bins=[size_bins,size_bins], normed=True)
        extent_size = [yedges_size[0], yedges_size[-1], xedges_size[-1], xedges_size[0]]
        plt.figure(1, figsize=(16,16))
        plt.subplot(1,2,1)
#        mat_size = (H_size/theo_array_size)
        mat_size = H_size
        mat_size[mat_size ==-np.inf] = 0
        plt.imshow(mat_size, extent=extent_size, interpolation='nearest')
        plt.axis([min(size_bins),max(size_bins),min(size_bins),max(size_bins)])
        cb = plt.colorbar(orientation = 'horizontal')
        cb.set_label('contact enrichment')
        plt.title('fragment length : enrichment map')
    ########################## gc content map ########################################
        H_gc, xedges_gc, yedges_gc = np.histogram2d(numpy_gc[:,0], numpy_gc[:,1],  bins=[gc_bins,gc_bins], normed=True)
        extent_gc = [yedges_gc[0], yedges_gc[-1], xedges_gc[-1], xedges_gc[0]]
        plt.subplot(1,2,2)
#        mat_gc = (H_gc/(theo_array_gc))
        mat_gc = H_gc
        mat_gc[mat_gc == -np.inf] = 0
        plt.imshow(mat_gc, extent=extent_gc, interpolation='nearest')
        plt.axis([min(gc_bins),max(gc_bins),min(gc_bins),max(gc_bins)])
        cb = plt.colorbar(orientation = 'horizontal')
        cb.set_label('contact enrichment')
        plt.title('fragment gc content : enrichment map')
        file_graph = os.path.join(output_folder,'gc_size_enrichment.png')
        plt.savefig(file_graph,dpi=None, facecolor='w', edgecolor='w',
            orientation='portrait', papertype=None, format=None,
            transparent=False, bbox_inches=None, pad_inches=0.1)
        plt.close()

        np.savetxt(file_numpy_mat_gc,mat_gc,delimiter='\t')
        np.savetxt(file_numpy_steps_gc,gc_bins,delimiter='\t')
        np.savetxt(file_numpy_mat_length,mat_size,delimiter='\t')
        np.savetxt(file_numpy_steps_length,size_bins,delimiter='\t')
        return mat_gc, mat_size, gc_bins, size_bins
    else:
        mat_gc = np.load(file_numpy_mat_gc)
        mat_size = np.load(file_numpy_mat_length)
        gc_bins = np.load(file_numpy_steps_gc)
        size_bins = np.load(file_numpy_steps_length)
        return mat_gc, mat_size, gc_bins, size_bins


def fragments_contacts_2_weighted_contacts(dict_cumul_length,dict_fragments,fragments_abs_contacts_files_weighted,fragments_contacts_files_weighted,fragments_contacts_file,mat_gc,mat_length,steps_gc,steps_size):
    import numpy as np
    import os
    mat_gc_sym = (mat_gc + mat_gc.transpose()) / float(2)
    norm_gc_mat = mat_gc_sym.sum()
    GC = (1/norm_gc_mat) * mat_gc_sym

    mat_length_sym = (mat_length + mat_length.transpose()) / float(2)
    norm_length_mat = mat_length_sym.sum()
    LGTH = (1/norm_length_mat) * mat_length_sym
    if not(os.path.exists(fragments_contacts_files_weighted) & os.path.exists(fragments_abs_contacts_files_weighted)):
        input_contact = open(fragments_contacts_file,'r')
        output_contact = open(fragments_contacts_files_weighted,'w')
        output_contact_abs = open(fragments_abs_contacts_files_weighted,'w')
        output_contact_abs.write("%s\t%s\t%s\t%s\n" %('id_read_a','id_read_b','w_length','w_gc'))
        output_contact.write("%s\t%s\t%s\t%s\t%s\t%s\n" %('id_read_a','contig_a','id_read_b','contig_b','w_length','w_gc'))

        while 1:
            line_a = input_contact.readline()
            if not line_a:
                input_contact.close()
                output_contact.close()
                output_contact_abs.close()
                print "done!"
                break


            data = line_a.split()
            id_read_a = data[0]
            contig_a = data[1]
            id_read_b = data[2]
            contig_b = data[3]

            ### RUSTINE a CHANGER !!!! CODE DE PORC ######
            frag_a = id_read_a+'-'+contig_a
            frag_b = id_read_b+'-'+contig_b

            if (frag_a !=frag_b) and dict_fragments.has_key(frag_a) and dict_fragments.has_key(frag_b):
                gc_a = np.float(dict_fragments[frag_a]["gc_content"])
                gc_b = np.float(dict_fragments[frag_b]["gc_content"])
                size_a = np.float(dict_fragments[frag_a]["size"])
                size_b = np.float(dict_fragments[frag_b]["size"])
    ##############################################
                tmp_gc_a = gc_a - steps_gc
                tmp_gc_b = gc_b - steps_gc
                gc_bin_a = np.nonzero(tmp_gc_a>=0)[0][-1]
                gc_bin_b = np.nonzero(tmp_gc_b>=0)[0][-1]
                if gc_bin_a == len(steps_gc)-1:
                    gc_bin_a = len(steps_gc)-2
                if gc_bin_b == len(steps_gc)-1:
                    gc_bin_b = len(steps_gc)-2
    ###############################################
                tmp_length_a = size_a - steps_size
                tmp_length_b = size_b - steps_size
                length_bin_a = np.nonzero(tmp_length_a>=0)[0][-1]
                length_bin_b = np.nonzero(tmp_length_b>=0)[0][-1]
                if length_bin_a == len(steps_size)-1:
                    length_bin_a = len(steps_size)-2
                if length_bin_b == len(steps_size)-1:
                    length_bin_b = len(steps_size)-2
    ################################################
                w_length = LGTH[length_bin_a,length_bin_b]
                w_gc = GC[gc_bin_a,gc_bin_b]
                output_contact.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(id_read_a,contig_a,id_read_b,contig_b,str(w_length),str(w_gc)))
    ################################################
                abs_a = int(id_read_a) + dict_cumul_length[contig_a]
                abs_b = int(id_read_b) + dict_cumul_length[contig_b]
                output_contact_abs.write("%s\t%s\t%s\t%s\n" %(str(abs_a),str(abs_b), str(w_length), str(w_gc)))
