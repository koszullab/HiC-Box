__author__ = 'hervemn'

from fragment import fragment

class contig():

    def __init__(self ,):
        "standart init "

    @classmethod
    def initiate(cls,init_contig,id_frag_start,id_frag_end,orientation):
        obj = cls()
        obj.name = obj.init_contig_name(init_contig,id_frag_start,id_frag_end,orientation)
        obj.length_kb = 0
        obj.n_frags = 0
        obj.kb_metric = [0]
        obj.FragmentS = []
        obj.orientation = orientation
        obj.status = 'active'
        obj.old_2_new = dict()
        obj.old_2_new['0:-1'] = {'contig_id':'self', 'position':1, 'orientation':'w'}
        obj.new_compo = []
        obj.new_compo.append({'contig_id':'self','position':1,'orientation':'w','start':1,'end':-1,'dup':False})
        obj.composition = []
        obj.contig_id = 0
        ####################
        init_compo = dict()

        init_compo['init_contig'] = init_contig
        init_compo['id_frag_start'] = id_frag_start
        init_compo['id_frag_end'] = id_frag_end
        init_compo['orientation'] = orientation
        init_compo['pos_id'] = 1
        #####################
        obj.composition.append(init_compo)
        return obj


    def to_string(self, reverse):
        import string
        out = ""
        if reverse:
            for i in range(len(self.FragmentS) - 1, -1, -1):
                f = self.FragmentS[i]
                out = out + f.init_name + f.orientation.translate(string.maketrans('wc', 'cw'))
        else:
            for f in self.FragmentS:
                out = out + f.init_name + f.orientation
        return out

    def init_contig_name(self,name,id_frag_start,id_frag_end,orientation):
        new_name = '>'+name+':'+str(id_frag_start)+':'+str(id_frag_end)+':'+ orientation +'>'
        return new_name

    def update_name(self,contig_id):
        self.contig_id = contig_id
        if self.old_2_new.has_key('0:-1'):
            if self.old_2_new['0:-1']['contig_id'] == 'self':
                self.old_2_new.pop('0:-1')
                start = str(1)
                end = str(self.n_frags)
                name_0 = start+':'+end
                self.old_2_new[name_0] = dict()
                self.old_2_new[name_0]['contig_id'] = contig_id
                self.old_2_new[name_0]['position'] = 1
                self.old_2_new[name_0]['orientation'] = 'w'
        for ele in self.new_compo:
            if ele['contig_id'] == 'self':
                ele['contig_id'] = contig_id
                ele['end'] = self.n_frags
        #update fragments name
        for f in self.FragmentS:
            f.update_name(contig_id)

    def update_repeats(self,old_instance_2_new_instance):
        """udpate repeats name in contig"""
        from colorama import Fore, Back, Style
        from colorama import init
        init(autoreset=True)
        print " -------------------------------"
        print Fore.BLUE + " update contig name:"+self.name
        for ele in self.new_compo:
            if ele["dup"]:
                old_contig_id = ele["contig_id"]
                new_contig_id = old_instance_2_new_instance[old_contig_id]["new_instance"]
                ele["contig_id"] = new_contig_id
                print "new_compo" + Fore.RED + str(old_contig_id) + " --->" +Fore.GREEN + str(new_contig_id)

    def append_fragment(self,fragment,orientation):
        import numpy as np
        self.kb_metric.append(self.kb_metric[-1] + fragment.length_kb/2)
        ind_frag = self.n_frags + 1
        fragment.curr_id = ind_frag
        fragment.orientation = orientation

        self.length_kb +=  fragment.length_kb

        fragment.pos_kb = np.copy(self.kb_metric[-1])
        self.FragmentS.append(fragment)
        self.n_frags = ind_frag

    def reset_old_2_new_new_compo(self):
        self.old_2_new = dict()
        access_bloc = "1:" + str(self.n_frags)
        self.old_2_new[access_bloc] = dict()
        self.old_2_new[access_bloc]["orientation"] = "w"
        self.old_2_new[access_bloc]["contig_id"] = self.contig_id
        self.old_2_new[access_bloc]["position"] = 1

        if self.new_compo[0]["dup"]:
            is_duplicate = True
        else:
            is_duplicate = False
        self.new_compo = []
        tmp_new_compo = dict()

        tmp_new_compo["contig_id"] = self.contig_id
        tmp_new_compo["orientation"] = "w"
        tmp_new_compo["position"] = 1
        tmp_new_compo["start"] = 1
        tmp_new_compo["end"] = self.n_frags
        tmp_new_compo["dup"] = is_duplicate
        self.new_compo.append(tmp_new_compo)
    def have_to_split(self, pos):

        if (pos == 0) or (pos >= self.n_frags - 1): # added 23 /02/2013 nfrags - 1
            output = False
        else:
            output = True
        return output

    def return_extremity(self,extremity,orientation):
        ## extremity = left or right
        sub_parts = self.old_2_new.keys()
        first_coords = []
        for ele in sub_parts:
            tmp = ele.split(':')
            first_coords.append(int(tmp[0]))
        new_index = sorted(range(len(first_coords)), key=lambda k: first_coords[k])
        print "contig name = ",self.name
        print new_index
        print sub_parts
        ext = extremity == 'right'
        ori = orientation == 'w'
        test = (ext and ori) or (not(ext) and not(ori))
        if test:
            return sub_parts[new_index[-1]]
        else:
            return sub_parts[new_index[0]]

    def vect_orientation(self,):
        orientation = []
        for ele in self.composition:
            n_frag = ele['id_frag_end'] - ele['id_frag_start'] + 1
            ori = ele['orientation']
            for i in range(0,n_frag):
                orientation.append(ori)
        return orientation

    def reverse_name(self,):
        import string
        data = self.name.split('-')
        n = len(data)
        reverso = range(n-1,-1,-1)
        cumul_init = []
        for ele in reverso:
            tmp = data[ele]
            bloc = tmp[:-2]+tmp[-2].translate(string.maketrans('wc', 'cw')) + '>'
            cumul_init.append(bloc)
        output = string.join(cumul_init,'-')
        return output

    def sub_content(self,pos_a,pos_b):
        # a modifier comme le split !!!!!!!!!!!!!!!!!!!!!!!!
        from copy import deepcopy as deepcopy
        ## composition split
        n_compo = len(self.composition)
        compo1 = []
        id_frag = 0
        f = deepcopy(self.FragmentS)
        frag1 = []
        stop = False
        position = 0
        n_frags_prec = 0
        for ele in range(0,n_compo):
            out_cond = True
            data = self.composition[ele]
            orientation = data['orientation']
            start = data['id_frag_start']
            start0 = data['id_frag_start']
            end = data['id_frag_end']
            tmp_compo = dict()
            tmp_compo['orientation'] = data['orientation']
            n_frags_local = end - start0 + 1

            if orientation == 'w':
            #                index_frag = range(start,end+1)
                index_frag = range(1,n_frags_local + 1)
            else:
            #                index_frag = range(end,start-1,-1)
                index_frag = range(n_frags_local,0,-1)

            tmp_compo['init_contig'] = data['init_contig']
            n_frags = end - start + 1
            id_rel = 0
            id_rel_2 = 0
            for i in index_frag:
                id_rel += 1
                id_frag += 1
                if ((id_frag>=pos_a) and (id_frag<=pos_b)):
                    id_rel_2 += 1
                    if id_rel_2 ==1:
                        position +=1
                        tmp_compo['pos_id'] = position
                    ## splitting fragments
                    frag1.append(f[n_frags_prec + i-1])

                    if id_frag == pos_a:

                        start = i
                        end = i
                    if orientation == 'w':
                        tmp_compo['id_frag_start'] = start
                        tmp_compo['id_frag_end'] = i + (start0 -1)
                    else:
                        tmp_compo['id_frag_start'] = i + (start0 -1)
                        tmp_compo['id_frag_end'] = end
                if id_frag == pos_b:
                    if (id_rel  == n_frags):
                        out_cond = False
                    stop = True
                    compo1.append(tmp_compo)
            if (out_cond and not(stop)):
                compo1.append(tmp_compo)
            n_frags_prec += n_frags_local
        ## new composition split
        new_n_compo = len(self.new_compo)
        new_compo1 = []
        id_frag = 0
        stop = False
        position = 0
        for ele in range(0,new_n_compo):
            out_cond = True
            data = self.new_compo[ele]
            orientation = data['orientation']
            start = data['start']
            start0 = data['start']
            end = data['end']
            tmp_compo = dict()
            tmp_compo['orientation'] = data['orientation']
            # the bloc has been duplicated force it here ??
            tmp_compo['dup'] = data['dup']
            if orientation == 'w':
            #                index_frag = range(start,end+1)
                index_frag = range(1,n_frags_local + 1)
            else:
            #                index_frag = range(end,start-1,-1)
                index_frag = range(n_frags_local,0,-1)
            tmp_compo['contig_id'] = data['contig_id']
            n_frags = end - start + 1
            id_rel = 0
            id_rel_2 = 0
            for i in index_frag:
                id_rel += 1
                id_frag += 1
                if ((id_frag>=pos_a) and (id_frag<=pos_b)):
                    id_rel_2 +=1
                    if id_rel_2 ==1:
                        position +=1
                        tmp_compo['position'] = position
                    if id_frag == pos_a:
                        start = i
                        end = i
                    if orientation == 'w':
                        tmp_compo['start'] = start
                        tmp_compo['end'] = i + (start0 -1)
                    else:
                        tmp_compo['start'] = i + (start0 -1)
                        tmp_compo['end'] = end
                if (id_frag == pos_b):
                    if (id_rel  == n_frags):
                        out_cond = False
                    new_compo1.append(tmp_compo)
                    stop = True
            if (out_cond and not(stop)):
                new_compo1.append(tmp_compo)

        return frag1, compo1, new_compo1

    def split_content(self,pos_cut):
        import string
        import numpy as np
        #### To correct !!! ####

        from copy import deepcopy as deepcopy
        ## composition split
        n_compo = len(self.composition)
        compo1 = []
        compo2 = []
        id_frag = 0
        active_compo = compo1
        f = deepcopy(self.FragmentS)
        frag1 = []
        frag2 = []
        active_frag = 1
        position = 0
        n_frags_prec = 0
        already_split = False
        for ele in range(0,n_compo):
            out_cond = True
            data = self.composition[ele]
            orientation = data['orientation']
            start = data['id_frag_start']
            start0 = data['id_frag_start']

            end = data['id_frag_end']
            tmp_compo = dict()
            tmp_compo['orientation'] = data['orientation']
            n_frags_local = end - start0 + 1

            if orientation == 'w':
            #                index_frag = range(start,end+1)
                index_frag = range(1,n_frags_local + 1)
            else:
            #                index_frag = range(end,start-1,-1)
                index_frag = range(n_frags_local,0,-1)


            tmp_compo['init_contig'] = data['init_contig']
            n_frags = end - start + 1
            id_rel = 0
            for i in index_frag:
                id_rel += 1
                id_frag += 1
                ## splitting fragments
                ## Attention la position de split correspond deja a l'index de coupure dans self.FragmentS
                ## donc pas besoin de changer la facon dde progresser dans la liste fragments
                if active_frag == 1:
#                    frag1.append(f[n_frags_prec + i-1])
                    frag1.append(f[id_frag - 1])
                elif active_frag == 2:
#                    frag2.append(f[n_frags_prec + i-1])
                    frag2.append(f[id_frag - 1])
                if id_rel == 1 :
                    position +=1
                    tmp_compo['pos_id'] = position
                if orientation == 'w':
                    tmp_compo['id_frag_start'] = start
                    tmp_compo['id_frag_end'] = i + (start0 -1)
                else:
                    tmp_compo['id_frag_start'] = i + (start0 -1)
                    tmp_compo['id_frag_end'] = end
                if id_frag == pos_cut:
                    if (id_rel  == n_frags):
                        out_cond = False
                    already_split = True# added rv :: be careful
                    active_compo.append(tmp_compo)
                    active_compo = compo2
                    position = 0
                    tmp_compo = dict()
                    tmp_compo['init_contig'] = data['init_contig']
                    tmp_compo['orientation'] = data['orientation']
                    position +=1
                    tmp_compo['pos_id'] = position
                    if orientation == 'w':
                        start = i + 1 + (start0 -1)
                    else:
                        end = i - 1 + (start0 -1)
                    position = 0 # added rv
                    active_frag = 2

            if out_cond:
                active_compo.append(tmp_compo)
            n_frags_prec += n_frags_local
        ## new composition split
        new_n_compo = len(self.new_compo)
#        print " new n compo = ", new_n_compo
#        print self.new_compo
        new_compo1 = []
        new_compo2 = []
        id_frag = 0
        new_active_compo = new_compo1
        position = 0
        already_split = False
        for ele in range(0,new_n_compo):
            out_cond = True
            data = self.new_compo[ele]
            orientation = data['orientation']
            start = data['start']
            start0 = data['start']
            end = data['end']
            tmp_compo = dict()
            tmp_compo['orientation'] = data['orientation']
            # print 'contitg id = ', data['contig_id']
            # print 'champs dup ', data['dup']
            tmp_compo['dup'] = data['dup']
            n_frags_local = end - start0 + 1
            if orientation == 'w':
            #                index_frag = range(start,end+1)
                index_frag = range(1,n_frags_local+1)
            else:
            #                index_frag = range(end,start-1,-1)
                index_frag = range(n_frags_local,0,-1)
            tmp_compo['contig_id'] = data['contig_id']
            n_frags = end - start + 1
            id_rel = 0
            for i in index_frag:
                id_rel += 1
                id_frag += 1
                if id_rel == 1 : # added rv :: be careful
                    position +=1
                    tmp_compo['position'] = position
                if orientation == 'w':
                    tmp_compo['start'] = start
                    tmp_compo['end'] = i + (start0 -1)
                else:
                    tmp_compo['start'] = i + (start0 -1)
                    tmp_compo['end'] = end
                if id_frag == pos_cut:
                    if (id_rel  == n_frags):
                        out_cond = False
                    already_split = True  # added rv :: be careful
                    new_active_compo.append(tmp_compo)
                    new_active_compo = new_compo2
                    tmp_compo = dict()
                    position = 0
#                    position = 1
                    tmp_compo['contig_id'] = data['contig_id']
                    tmp_compo['orientation'] = data['orientation']
                    tmp_compo['dup'] = data['dup']
                    position +=1
                    tmp_compo['position'] = position
                    if orientation == 'w':
                        start = i + 1 + (start0 -1)
                    else:
                        end = i - 1 + (start0 -1)
                    position = 0 # added rv :: be careful
            if out_cond:
                new_active_compo.append(tmp_compo)
        return frag1, compo1, new_compo1, frag2, compo2, new_compo2


    def split_name(self,pos_cut):
        import string
        data = self.name.split('-')
        n = len(data)
        pos_in_contig = 0
        output = []
        check_cut = 0
        output.append([])
        output.append([])
        for ele in range(0,n):
            tmp_ctg = data[ele]
            info = tmp_ctg[1:-1].split(':') # initials contigs extract info
            init_contig = info[0]
            start = int(info[1])
            end = int(info[2])
            orientation = info[3]
            n_frags = (end - start) +1
            if orientation == 'w':
                index_frag = range(1,n_frags+1)
                start_write = start
            else:
                index_frag = range(n_frags,0,-1)
                start_write = end
            id = 0
            out_cond = True
            for i in index_frag:
                id +=1
                if orientation =='w':
                    name = '>' + string.join([init_contig,str(start_write),str(i + (start - 1) ),orientation],':') + '>'
                else:
                    name = '>' + string.join([init_contig,str(i + (start - 1) ),str(start_write),orientation],':') + '>'

                abs_pos = id + pos_in_contig

                if (abs_pos == pos_cut):
                    if (id == n_frags):
                        out_cond = False
                    output[check_cut].append(name)
                    check_cut = 1
                    output[check_cut] = []
                    if orientation == 'w':
                        start_write = i + 1 + (start - 1)
                    else:
                        start_write = i - 1 + (start - 1)

            if out_cond:
                output[check_cut].append(name)
            pos_in_contig = pos_in_contig + n_frags

        name_1 = string.join(output[0],'-')
        name_2 = string.join(output[1],'-')
        return name_1,name_2


    def sub_name(self,pos_a,pos_b):
        import string
        data = self.name.split('-')
        n = len(data)
        pos_in_contig = 0
        output = []
        stop = False
        for ele in range(0,n):
            tmp_ctg = data[ele]
            info = tmp_ctg[1:-1].split(':') # initials contigs extract info
            init_contig = info[0]
            start = int(info[1])
            end = int(info[2])
            orientation = info[3]
            n_frags = (end - start) +1
            if orientation == 'w':
                index_frag = range(1,n_frags+1)
                start_write = start
            else:
                index_frag = range(n_frags,0,-1)
                start_write = end
            id = 0
            out_cond = True
            for i in index_frag:
                id +=1
                abs_pos = id + pos_in_contig
                if ((abs_pos >= pos_a) and (abs_pos <= pos_b)):
                    if abs_pos == pos_a:
                        start_write = i
                    if orientation =='w':
                        name = '>' + string.join([init_contig,str(start_write),str(i),orientation],':') + '>'
                    else:
                        name = '>' + string.join([init_contig,str(i),str(start_write),orientation],':') + '>'

                if (abs_pos == pos_b):
                    if (id == n_frags):
                        out_cond = False
                    output.append(name)
                    stop = True

            if (out_cond and not(stop)):
                output.append(name)
            pos_in_contig = pos_in_contig + n_frags

        name_1 = string.join(output,'-')
        return name_1

    def reverse_content(self,):
        from copy import deepcopy as deepcopy
        import string
        init_contigs_compo = deepcopy(self.composition)
        new_contigs_compo = deepcopy(self.new_compo)
        FragmentS = deepcopy(self.FragmentS)
        new_contigs_index = 0
        init_contigs_index = 0
        fragment_index = 0
        init_contigs_compo.reverse()
        for ele in init_contigs_compo:
            init_contigs_index += 1
            ele['orientation'] = ele['orientation'].translate(string.maketrans('wc', 'cw'))
            ele['pos_id'] = init_contigs_index

        new_contigs_compo.reverse()
        for ele in new_contigs_compo:
            new_contigs_index += 1
            ele['orientation'] = ele['orientation'].translate(string.maketrans('wc', 'cw'))
            ele['position'] = new_contigs_index

        FragmentS.reverse()
        for ele in FragmentS:
            fragment_index +=1
            ele.curr_id = fragment_index

        result = dict()
        result['init_contig_compo'] = init_contigs_compo
        result['new_contig_compo'] = new_contigs_compo
        result['FragmentS'] = FragmentS
        return result




    def export_abs_frag_vect(self):
        abs_frag_vect = []
        for fragment in self.FragmentS:
            abs_frag_vect.append(fragment.np_id_abs-1)
        return abs_frag_vect

    def export_kb_pos_and_ctg_id(self,):
        kb_pos_and_ctg_id = []
        for fragment in self.FragmentS:
            kb_pos_and_ctg_id.append([fragment.pos_kb,fragment.contig_id])
        return kb_pos_and_ctg_id




    def dist(self,frag_a,frag_b):
        return self.kb_metric

    @classmethod
    def join(cls,contig_a,do_reverse_a,contig_b,do_reverse_b):
        from copy import deepcopy as deepcopy
        import string
        obj = cls()
        tmp_name = []
        if do_reverse_a:
            tmp_name.append(contig_a.reverse_name())
            result_a = contig_a.reverse_content()
            FragmentS_a = result_a['FragmentS']
            Compo_init_a = result_a['init_contig_compo']
            ind_init_compo = len(Compo_init_a)
            Compo_new_a = result_a['new_contig_compo']
            ind_new_compo = len(Compo_new_a)
        else:
            tmp_name.append(contig_a.name)
            FragmentS_a = deepcopy(contig_a.FragmentS)
            Compo_init_a = deepcopy(contig_a.composition)
            ind_init_compo = len(Compo_init_a)
            Compo_new_a = deepcopy(contig_a.new_compo)
            ind_new_compo = len(Compo_new_a)


        if do_reverse_b:
            tmp_name.append(contig_b.reverse_name())
            result_b = contig_b.reverse_content()
            FragmentS_b = result_b['FragmentS']
            Compo_init_b = result_b['init_contig_compo']
            for ele in Compo_init_b:
                ind_init_compo +=1
                ele['pos_id'] = ind_init_compo

            Compo_new_b = result_b['new_contig_compo']
            for ele in Compo_new_b:
                ind_new_compo +=1
                ele['position'] = ind_new_compo
        else:
            tmp_name.append(contig_b.name)
            FragmentS_b = deepcopy(contig_b.FragmentS)
            Compo_init_b = deepcopy(contig_b.composition)
            for ele in Compo_init_b:
                ind_init_compo +=1
                ele['pos_id'] = ind_init_compo
            Compo_new_b = deepcopy(contig_b.new_compo)
            for ele in Compo_new_b:
                ind_new_compo +=1
                ele['position'] = ind_new_compo

        obj.name = string.join(tmp_name,'-')
        obj.length_kb = 0
        obj.n_frags = 0
        obj.kb_metric = [0]
        obj.FragmentS = []
        obj.orientation = 'w'
        obj.status = 'active'
        obj.old_2_new = dict()
        obj.old_2_new['0:-1'] = {'contig_id':'self','position':1,'orientation':'w'}
        obj.new_compo = Compo_new_a
        obj.new_compo.extend(Compo_new_b)
        obj.composition = Compo_init_a
        obj.composition.extend(Compo_init_b)
        FragmentS_a.extend(FragmentS_b)
        orientation = obj.vect_orientation()
        i = 0
        for frag in FragmentS_a:
            new_frag = fragment.copy(frag)
            obj.append_fragment(new_frag,orientation[i])
            i += 1
        return obj

    @classmethod
    def split(cls,contig_a, pos_cut):
        # print "splitting"
        obj1 = cls()
        obj2 = cls()

        name1,name2 = contig_a.split_name(pos_cut)
        frag1, compo1, new_compo1, frag2, compo2, new_compo2 = contig_a.split_content(pos_cut)
        # corriger le split ici ...
        # print 'new_compo1 = ',new_compo1
        # print 'new_compo2 = ',new_compo2

        obj1.name = name1
        obj1.length_kb = 0
        obj1.n_frags = 0
        obj1.kb_metric = [0]
        obj1.FragmentS = []
        obj1.orientation = 'w'
        obj1.status = 'active'
        obj1.old_2_new = dict()
        obj1.old_2_new['0:-1'] = {'contig_id':'self','position':1,'orientation':'w'}
        obj1.new_compo = new_compo1
        obj1.composition = compo1
        orientation = obj1.vect_orientation()
        i = 0

        for frag in frag1:
            obj1.append_fragment(frag,orientation[i])
#            print "contig 1 ultimate frag added = ",obj1.FragmentS[-1].curr_name
            i += 1


        obj2.name = name2
        obj2.length_kb = 0
        obj2.n_frags = 0
        obj2.kb_metric = [0]
        obj2.FragmentS = []
        obj2.orientation = 'w'
        obj2.status = 'active'
        obj2.old_2_new = dict()
        obj2.old_2_new['0:-1'] = {'contig_id':'self','position':1,'orientation':'w'}
        obj2.new_compo = new_compo2
        obj2.composition = compo2
        orientation = obj2.vect_orientation()
        i = 0

        # print '##################################'


        for frag in frag2:
#            print frag.init_name
            obj2.append_fragment(frag,orientation[i])
#            print "contig 2 ultimate frag added = ",obj2.FragmentS[-1].curr_name
            i += 1

        return obj1, obj2


    @classmethod
    def duplicate(cls,contig_a,pos_a,pos_b):
        print "duplicate"
        obj1 = cls()
        name1= contig_a.sub_name(pos_a,pos_b)
        frag1, compo1, new_compo1 = contig_a.sub_content(pos_a,pos_b)
        obj1.name = name1
        obj1.length_kb = 0
        obj1.n_frags = 0
        obj1.kb_metric = [0]
        obj1.FragmentS = []
        obj1.orientation = 'w'
        obj1.status = 'active'
        obj1.old_2_new = dict()
        obj1.old_2_new['0:-1'] = {'contig_id': 'self', 'position': 1, 'orientation':'w', 'dup': False}
        obj1.new_compo = new_compo1
        obj1.composition = compo1
        orientation = obj1.vect_orientation()
        i = 0
        for frag in frag1:
            obj1.append_fragment(frag,orientation[i])
            i += 1
        return obj1



    def has_repeats(self):
        output = False
        for ele in self.new_compo:
            if ele["dup"]:
                output = True
        return output
