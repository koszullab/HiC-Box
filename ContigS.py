__author__ = 'hervemn'

from mutation import new_mutation as mutation
import numpy as np

class repeats():
    def __init__(self):
        """ initialization """
    @classmethod
    def init_from_initial_structure(cls, ContigS, contig_a, n_copy):
        from copy import deepcopy as deepcopy
        # duplicated bloc instance has an id which will be tracked
        obj = cls()
        obj.ContigS = ContigS
        obj.name = contig_a.name
        obj.init_contig = contig_a.composition[0]['init_contig']
        obj.id_frag_start = contig_a.composition[0]['id_frag_start']
        obj.id_frag_end = contig_a.composition[0]['id_frag_end']
        obj.orientation = 'w'
        obj.access_bloc = str(obj.id_frag_start)+':'+str(obj.id_frag_end)
        print " access bloc = ",obj.access_bloc
        obj.list_fragments = deepcopy(contig_a.FragmentS)
        obj.current_positions = dict() # dictionnary of current positions of the duplicates
        obj.incoming_connections = [] # list of new paste ( in  the future split and paste)
        obj.r_x = n_copy # number of repeats
        obj.ContigS.init_contigs[obj.init_contig][obj.access_bloc]['dup'] = True # set this area as duplicated
        obj.ContigS.inactivate(contig_a)
        return obj

    @classmethod
    def init_from_previous_repeat(cls, old_repeat, ContigS):
        from copy import deepcopy as deepcopy
        """ initiate repeats and instances from previous data"""
        obj = cls()
        obj.ContigS = ContigS
        obj.name = old_repeat.name
        obj.init_contig = old_repeat.init_contig
        obj.id_frag_start = old_repeat.id_frag_start
        obj.id_frag_end = old_repeat.id_frag_end
        obj.orientation = old_repeat.orientation
        obj.access_bloc = old_repeat.access_bloc
        obj.list_fragments = deepcopy(old_repeat.list_fragments)
        # building translation dictionary from old instances to new instances
        obj.current_positions = dict() # dictionnary of current positions of the duplicates
        obj.incoming_connections = [] # list of new paste ( in  the future split and paste)
        obj.r_x = old_repeat.r_x # number of repeats
        if not(obj.ContigS.init_contigs.has_key(obj.init_contig)):
            obj.ContigS.init_contigs[obj.init_contig] = dict()
        if not(obj.ContigS.init_contigs[obj.init_contig].has_key(obj.access_bloc)):
            obj.ContigS.init_contigs[obj.init_contig][obj.access_bloc] = dict()
        obj.ContigS.init_contigs[obj.init_contig][obj.access_bloc]['dup'] = True # set this area as duplicated
        return obj


    def add_position(self,id_instance,new_pos):
        self.current_positions[id_instance] = new_pos

    def is_instance_free(self,id_instance):
        if self.current_positions[id_instance] == id_instance:
            output = True
        else:
            output = False
        return output

    def remove_position(self,id_instance):
        new_pos = dict()
        new_pos['orientation'] = 'w'
        new_pos['position'] = 1
        new_pos['@'] = id_instance
        self.current_positions[id_instance] = new_pos

    def update_current_position(self, id_instance, new_pos):
        """update positions in class"""
        self.current_positions[id_instance] = new_pos

    def update_incoming_connections(self,connection):
        """ add incoming connection """
    def create_contig_instance(self,):
        """ create instance from template"""
        from copy import deepcopy as deepcopy
        from contig import contig as contig
        instance = contig.initiate(self.init_contig,self.id_frag_start,self.id_frag_end,self.orientation)
        new_list_frag = deepcopy(self.list_fragments)
        for frag in new_list_frag:
            instance.append_fragment(frag,'w')
#            print "frag current name = ",frag.curr_name
        id_instance = self.ContigS.add(instance)
        new_pos = dict()
        new_pos['orientation'] = 'w'
        new_pos['position'] = 1
        new_pos['@'] = id_instance

        instance.new_compo[0]["dup"] = True # importnat d etre defini apres
        self.add_position(id_instance,new_pos)
        return id_instance

    def free_all_instances(self,):
        """ destroy instances from template"""





class duplicator():
    def __init__(self,ContigS):
        """ initialization """
        self.duplicated_blocs = dict()
        self.ContigS = ContigS
        self.incoming_connections = []
        self.old_instance_2_new_instance = dict()

    def define_repeat(self,C,pos_a,pos_b,n_copy):
        """ """
        from contig import contig as contig
        # 1) definir le bloc duplique

        # 2) l extraire = split x 2 + delete du bloc
        # 3) creer une instance
        # 4) inserer l'instance sa position initiale par paste
        ###################
        out_split_a = self.ContigS.split(C,pos_a)  # modifier split afin de mettre a jour duplicator : pas urgent/presque inutile

#        print "name upstream bloc 1 = ",out_split_a['upstream']["contig"].name
#        print "name downstream bloc 1 = ",out_split_a['downstream']["contig"].name


        up_C = out_split_a['upstream']["contig"]
        tmp_C = out_split_a['downstream']["contig"]
        out_split_b = self.ContigS.split(C,pos_b)
        dup_bloc = out_split_b['upstream']["contig"] # duplicated bloc
#        print "name upstream bloc = ",out_split_b['upstream']["contig"].name
#        print "name downstream bloc = ",out_split_b['downstream']["contig"].name

        down_C = out_split_b['downstream']["contig"]
        ###################
        repeat = repeats.init_from_initial_structure(self.ContigS, dup_bloc, n_copy)
        self.duplicated_blocs[repeat.name] = repeat
        for i in range(0,n_copy):
            instance_repeat_id = repeat.create_contig_instance()
            self.old_instance_2_new_instance[instance_repeat_id] = dict()
            self.old_instance_2_new_instance[instance_repeat_id]["new_instance"] = instance_repeat_id
            self.old_instance_2_new_instance[instance_repeat_id]["remapped"] = True
        instance_repeat = self.ContigS.current_contigs[instance_repeat_id]
        up_C_dup = self.ContigS.paste(up_C,'w',instance_repeat,'w') # modifier paste afin de mettre a jour duplicator
        original = self.ContigS.paste(up_C_dup,'w',down_C,'w') # modifier paste afin de mettre a jour duplicator
        # pop out reference to duplicated part in original contig
        access_bloc_dup = str(pos_a + 1 ) +':' + str(pos_b) # be carefull of the split positions
        C.old_2_new.pop(access_bloc_dup)

    def define_repeat_from_old(self, old_repeat,):
        repeat = repeats.init_from_previous_repeat(old_repeat,self.ContigS)
        self.duplicated_blocs[repeat.name] = repeat
        for ele in old_repeat.current_positions:
            old_instance = ele
            self.old_instance_2_new_instance[old_instance] = dict()
            new_instance = repeat.create_contig_instance()
            self.old_instance_2_new_instance[old_instance]["new_instance"] = new_instance
            self.old_instance_2_new_instance[old_instance]["remapped"] = False




    def delete_duplicated_bloc(self,dup):
        """"""

    def add_modification_to_network(self,incoming_connection):
        self.incoming_connections.append(incoming_connection)

    def path_generator(self):

        """ return paste and split operations to perform on duplicated blocs"""
        """ update current network """

class ContigS():
#    from contig import contig as contig
    def __init__(self):
        import numpy as np
        import contig as CONTI
        self.duplicator = duplicator(self,)
        self.current_contigs = dict() # dictionary of current contigs entry : 1 -- n_contigs
        self.init_contigs = dict() # dictionnary of initial contigs : name given by de-novo assembler
        self.g_sigma = dict() # dictionnary of initial fragments
        self.np_id_2_frag = dict() # dictionary conversion of matrix index to init fragment name
        self.n_contigs = 0 # number of contigs
        self.size_current_genome = 0
        # C-structures sent to GPU !
        self.frag_T = np.dtype([('mult', np.int32), ('pos_min',np.int32), ('pos_max',np.int32), ('alpha',np.int32)], align = True)
        self.coord_T = np.dtype([('chrom', np.int32),
            ('pos_kb', np.float32),
            ('low', np.int32),
            ('high', np.int32),
            ('circular', np.int32),
            ('id_rel', np.int32),
            ('size_contig', np.int32),
            ('cumul_index', np.int32),
            ('alpha', np.float32)], align=True)
        self.alpha_collect_T = np.dtype([('obsd', np.float64), ('expd', np.float64)], align=True)

    def compute_size_current_genome(self):
        from contig import contig
        size_current_genome = 0
        for c in self.current_contigs.keys():
            cont = self.current_contigs[c]
            if isinstance(cont,contig):
                if cont.status == 'active':
                    size_current_genome += cont.n_frags
        self.size_current_genome = size_current_genome

    def clear(self):
        self.current_contigs.clear()
        self.init_contigs.clear()
        self.g_sigma.clear()
        self.np_id_2_frag.clear()
        self.n_contigs.clear()

    def add_new_composite(self,contig):
        from colorama import Fore, Back, Style
        from colorama import init
        init(autoreset=True)
        from contig import contig as class_contig
        """ perform adding of a contig which contain repeats"""
        curr_contig = contig
        pos_cut = 0
        list_id_sub_contig = []
#        for ele in contig.composition:
        test_condition = True
        curr_pos = 0
        jump_bloc = True
        while test_condition:
#            for ele in curr_contig.composition:

            ele = curr_contig.composition[curr_pos]
            # variable pour sauter de contig
            init_contig = ele["init_contig"]
            access_bloc = str(ele["id_frag_start"])+':'+str(ele["id_frag_end"])
            length_bloc = ele["id_frag_end"] - ele["id_frag_start"] + 1
#            if not(self.init_contigs.has_key(init_contig)):
#                self.init_contigs[init_contig] = dict()

            if self.init_contigs[init_contig].has_key(access_bloc):
                if self.init_contigs[init_contig][access_bloc]["dup"]:
                    bloc_upstream, tmp_bloc = class_contig.split(curr_contig,pos_cut)
                    # add bloc_upstream
                    if not(bloc_upstream.name == ''):
#                        print "name upstream = ",bloc_upstream.name
#                        print " new compo = ", bloc_upstream.new_compo

                        if bloc_upstream.has_repeats() and len(bloc_upstream.new_compo) ==1 :# <- check if has repeats!!!
#                            bloc_upstream.update_repeats(self.duplicator.old_instance_2_new_instance)
                            old_instance = bloc_upstream.new_compo[0]["contig_id"]
                            id_upstream = self.duplicator.old_instance_2_new_instance[old_instance]["new_instance"]
                        else:
                            id_upstream = self.add(bloc_upstream)
#                            self.current_contigs[id_upstream].new_compo[0]["contig_id"] = id_upstream
                            self.current_contigs[id_upstream].reset_old_2_new_new_compo()


                        list_id_sub_contig.append(id_upstream)
#                        print " tmp bloc name = ",tmp_bloc.name

                    if not(tmp_bloc.name == ''):
#                        print "name tmp = ",tmp_bloc.name
                        dup_bloc, bloc_downstream = class_contig.split(tmp_bloc, length_bloc)
                        # nouveau contig defini
                        curr_contig = bloc_downstream
                        curr_pos = 0
                        pos_cut = 0
                        if not(dup_bloc.name == ''):

                            if dup_bloc.has_repeats() and len(dup_bloc.new_compo) ==1:# <- check if has repeats!!!
                                old_instance = dup_bloc.new_compo[0]["contig_id"]
                                id_bloc_dup = self.duplicator.old_instance_2_new_instance[old_instance]["new_instance"]
#                                dup_bloc.update_repeats(self.duplicator.old_instance_2_new_instance)
                            else:
                                id_bloc_dup = self.add(dup_bloc)
#                                self.current_contigs[id_bloc_dup].new_compo[0]["contig_id"] = id_bloc_dup
                                self.current_contigs[id_bloc_dup].reset_old_2_new_new_compo()
                            list_id_sub_contig.append(id_bloc_dup)
                    else:
                        # nouveau contig defini
                        curr_contig = tmp_bloc
                        curr_pos = 0
                        pos_cut = 0
                else:
                    pos_cut += length_bloc
                    curr_pos += 1
            else:
                curr_pos += 1
                pos_cut += length_bloc
            len_curr_cont = len(curr_contig.composition)
            if curr_contig.name == '':
                test_condition = False
            else:
                if (curr_pos <len_curr_cont):
                    test_condition = True
                else:
                    test_condition = False


        if not(curr_contig.name == '' ):
#            print Fore.BLUE + " on est dedans!"
#            print "name last contig = ", curr_contig.name
            if curr_contig.has_repeats() and len(curr_contig.new_compo) == 1:# <- check if has repeats!!!
#                curr_contig.update_repeats(self.duplicator.old_instance_2_new_instance)
                old_instance = curr_contig.new_compo[0]["contig_id"]
                id_bloc = self.duplicator.old_instance_2_new_instance[old_instance]["new_instance"]
            else:
                id_bloc = self.add(curr_contig)
#                self.current_contigs[id_bloc].new_compo[0]["contig_id"] = id_bloc
                self.current_contigs[id_bloc].reset_old_2_new_new_compo()
            list_id_sub_contig.append(id_bloc)
            curr_contig = self.current_contigs[id_bloc]

#            print " length contig = ",curr_contig.n_frags
#            print " last contig id = ", curr_contig.contig_id




        output = self.current_contigs[list_id_sub_contig[0]]


        len_sub_contig = len(list_id_sub_contig)
#        print "list sub contig = ", list_id_sub_contig
#        print " len sub contig = ", len_sub_contig

        if len_sub_contig>1:
            for i in range(1,len(list_id_sub_contig)):
                sub_contig = list_id_sub_contig[i]
                print Fore.BLUE + output.name
                print Fore.RED + self.current_contigs[sub_contig].name
                output_tmp = self.paste(output,'w',
                    self.current_contigs[sub_contig],'w')
                output = output_tmp

        return output.contig_id


    def add(self,contig):
        self.n_contigs +=1
        # add contig into current_contigs dictionary
        contig.update_name(self.n_contigs)
        self.current_contigs[self.n_contigs] = contig

        for fragment in contig.FragmentS:
            # update np_id_2_frag
#            fragment.update_name(self.n_contigs)
            self.np_id_2_frag[fragment.np_id_abs] = fragment.init_name
            # update g_sigma
            if fragment.init_name in self.g_sigma:
                self.g_sigma[fragment.init_name].add(fragment.curr_name)
            else:
                self.g_sigma[fragment.init_name] = set()
                self.g_sigma[fragment.init_name].add(fragment.curr_name)
        #update init_contigs
        for ele in contig.composition:
            init_contig = ele['init_contig']
            if init_contig in self.init_contigs:
                position = str(ele['id_frag_start'])+':' + str(ele['id_frag_end'])
                if not(self.init_contigs[init_contig].has_key(position)):
                    self.init_contigs[init_contig][position] = dict()
                    self.init_contigs[init_contig][position]['dup'] = False
                bloc_tag = str(ele['pos_id'])+'@'+str(contig.contig_id)
                self.init_contigs[init_contig][position][bloc_tag] = str(ele['orientation'])
#                is_duplicated = self.init_contigs[init_contig][position]['dup']
#                if is_duplicated: # update duplicator structure
#                    name_repeat = '>'+init_contig+':'+position+':w>'
#                    repeat = self.duplicator.duplicated_blocs[name_repeat]
#                    repeat.add_position(bloc_tag)
            else:
                self.init_contigs[init_contig] = dict()
                position = str(ele['id_frag_start'])+':'+str(ele['id_frag_end'])
                self.init_contigs[init_contig][position] = dict()
                bloc_tag = str(ele['pos_id'])+'@'+str(contig.contig_id)
                self.init_contigs[init_contig][position][bloc_tag] = str(ele['orientation'])
                self.init_contigs[init_contig][position]['dup'] = False # track de la signature de duplication

        for ele in contig.new_compo:

            if ele["dup"]:
                position = ele["position"]
                orientation = ele["orientation"]
                tmp_id_instance = ele["contig_id"]
                # id_instance must be renew only at the beginning
                id_instance = tmp_id_instance
#                if self.duplicator.old_instance_2_new_instance[tmp_id_instance]["remapped"]:
#                    id_instance = tmp_id_instance
#                else:
#                    id_instance = self.duplicator.old_instance_2_new_instance[tmp_id_instance]["new_instance"]
#                    self.duplicator.old_instance_2_new_instance[tmp_id_instance]["remapped"] = True
##                    self.duplicator.old_instance_2_new_instance[id_instance]["new_instance"] = id_instance

                print "id instance = ", id_instance
                ele["contig_id"] = id_instance

                self.inactivate(self.current_contigs[id_instance]) # desactivation : dans tous les cas il est rajoute

                name_repeat = self.current_contigs[id_instance].name # problem when going from genome t to genome ( t+1 )
                repeat = self.duplicator.duplicated_blocs[name_repeat]
                new_pos = dict()
                new_pos['orientation'] = orientation
                new_pos['position'] = position
                new_pos['@'] = self.n_contigs
                repeat.update_current_position(id_instance,new_pos)

        self.compute_size_current_genome()
        return self.n_contigs

    def create_transparent_contig(self,):
        from contig import contig as contig
        transpa = contig.initiate("None",0,0,"w")
        return transpa

    def export_name_contigs(self,):
        out = []
        for ele in self.current_contigs:
            curr_cont = self.current_contigs[ele]
            if curr_cont.status == "active":
                out.append(curr_cont.name)
        return out

    def print_active_contig(self):
        from colorama import Fore, Back, Style
        from colorama import init
        init(autoreset=True)
        print Fore.BLUE + "list active contig"
        print "-----------------"
        for ele in self.current_contigs:
            curr_cont = self.current_contigs[ele]
            if curr_cont.status == "active":
                print Fore.RED + str(curr_cont.contig_id) + Fore.WHITE +"  ==>  " + Fore.GREEN + curr_cont.name
        print "-----------------"
        init(autoreset=True)
    ####################################################################################
    ################# basic mutations paste split duplicate inactivate activate ########

    def inactivate(self, contig):
        if self.current_contigs[contig.contig_id].status == 'active':

            # print " on va inactiver les fragments et les retirer de g sigma. contig id = ",contig.contig_id
            for fragment in contig.FragmentS:
#                print " fragment curr name = " ,fragment.curr_name
#                print " contenu g sigma = ", self.g_sigma[fragment.init_name]
#                print "removing : ", fragment.curr_name
                self.g_sigma[fragment.init_name].remove(fragment.curr_name)
                try:
                    self.np_id_2_frag.pop(fragment.np_id_abs) # attention ....
                except:
                    ""
            for ele in contig.composition:

                init_contig = ele['init_contig']
                position = str(ele['id_frag_start'])+':'+str(ele['id_frag_end'])
                bloc_tag = str(ele['pos_id'])+'@'+str(contig.contig_id)

                self.init_contigs[init_contig][position].pop(bloc_tag)

                if self.init_contigs[init_contig][position] == {'dup':False}:
                    self.init_contigs[init_contig].pop(position)

            for ele in contig.new_compo:
                if ele["dup"]:
                    id_instance = ele["contig_id"]
#                    print "id instance = ", id_instance
                    name_repeat = self.current_contigs[id_instance].name
#                    print " name repeat = ", name_repeat
                    repeat = self.duplicator.duplicated_blocs[name_repeat]
                    repeat.remove_position(id_instance)

            self.current_contigs[contig.contig_id].status = 'inactive'

        self.compute_size_current_genome()

    def reactivate(self,contig):
        if self.current_contigs[contig.contig_id].status == 'inactive':
            contig.status = 'active'
            for fragment in contig.FragmentS:
                # update np_id_2_frag
                fragment.update_name(self.n_contigs)
                self.np_id_2_frag[fragment.np_id_abs] = fragment.init_name
                # update g_sigma
                if self.g_sigma.has_key(fragment.init_name):
                    self.g_sigma[fragment.init_name].add(fragment.curr_name)
                else:
                    self.g_sigma[fragment.init_name] = set()
                    self.g_sigma[fragment.init_name].add(fragment.curr_name)
                #update init_contigs
            for ele in contig.composition:
                init_contig = ele['init_contig']
                if self.init_contigs.has_key(init_contig):
                    position = str(ele['id_frag_start'])+':'+str(ele['id_frag_end'])
                    if not(self.init_contigs[init_contig].has_key(position)):
                        self.init_contigs[init_contig][position] = dict()
                        self.init_contigs[init_contig][position]["dup"] = False # added 23 /02 /2013
                    bloc_tag = str(ele['pos_id'])+'@'+str(contig.contig_id)
                    self.init_contigs[init_contig][position][bloc_tag] = str(ele['orientation'])

                else:
                    self.init_contigs[init_contig] = dict()
                    position = str(ele['id_frag_start'])+':'+str(ele['id_frag_end'])
                    self.init_contigs[init_contig][position] = dict()
                    self.init_contigs[init_contig][position]['dup'] = False
                    bloc_tag = str(ele['pos_id'])+'@'+str(contig.contig_id)
                    self.init_contigs[init_contig][position][bloc_tag] = str(ele['orientation'])
            for ele in contig.new_compo:
                if ele["dup"]:
                    position = ele["position"]
                    orientation = ele["orientation"]
                    id_instance = ele["contig_id"]
                    name_repeat = self.current_contigs[id_instance].name
                    repeat = self.duplicator.duplicated_blocs[name_repeat]
                    new_pos = dict()
                    new_pos['orientation'] = orientation
                    new_pos['position'] = position
                    new_pos['@'] = contig.contig_id
                    repeat.add_position(id_instance,new_pos)

            self.compute_size_current_genome()
            return self.n_contigs

    def duplicate(self,contig_a, pos_a,pos_b):
        from contig import contig as contig
        print "duplicate"
        test_a = pos_a>=1
        test_b = pos_b<=contig_a.n_frags
        if test_a and test_b:
            print "duplication possible"
            bloc_duplicated = contig.duplicate(contig_a,pos_a,pos_b)
            for ele in bloc_duplicated.new_compo:
                ele['dup'] = True

            id_bloc = self.add(bloc_duplicated)
        else:
            print "duplication impossible"
            id_bloc = 0
        return self.current_contigs[id_bloc]


    def sort_blocs_in_old_2_new(self,contig_a):
        sub_parts_tmp = contig_a.old_2_new.keys()
        dico_tmp = dict()
        for ele in sub_parts_tmp:
            bla = int(ele.split(':')[0])
            dico_tmp[bla] = ele
        tmp_sub=[]
        tmp_sub.extend(dico_tmp.keys())
        tmp_sub.sort()
        sub_parts = []
        for ele in tmp_sub:
            sub_parts.append(dico_tmp[ele])
        return sub_parts

    def pos_in_bloc_2_global_pos(self,contig,pos_in_bloc,bloc_id):
        output = dict()
        t = range(0,bloc_id - 1)
        cumul_length = 0
        print "contig new compo = ", contig.new_compo
        print "contig name = ", contig.name
        print "pos in bloc :", pos_in_bloc
        print "bloc id = ", bloc_id
        for h in t:
            start = contig.new_compo[h]['start']
            end = contig.new_compo[h]['end']
            len_bloc = end - start + 1
            cumul_length += len_bloc
        start = contig.new_compo[bloc_id - 1]['start']
        end = contig.new_compo[bloc_id - 1]['end']
        len_bloc = end - start + 1
        orientation = contig.new_compo[bloc_id - 1]["orientation"]
        if orientation == 'w':
            output["reverse_order"] = False
            new_pos = cumul_length + pos_in_bloc
            output["new_pos"] = new_pos
            output["orientation"] = 'w'
        else:
            new_pos = cumul_length + (len_bloc - pos_in_bloc) + 1
            output["reverse_order"] = True
            output["new_pos"] = new_pos - 1 # because split will take the upper part of the contig
            output["orientation"] = 'c'
        return output

    def find_new_position(self, contig_a, pos_a):
        sorted_blocs = self.sort_blocs_in_old_2_new(contig_a)
        output = dict()
        output["upstream"] = dict()
        output["downstream"] = dict()
        output["safe_split"] = True
        output["classic"] = dict()
        i = 0
        sum_size_previous_blocs = 0
        size_previous_blocs = [0]

        print "pos_a = ", pos_a
        print "sorted blocs = "
        print sorted_blocs
        for ele in sorted_blocs:
            print "sub parts = ", ele
            tmp = ele.split(':')
            low_lim = int(tmp[0])
            high_lim = int(tmp[1])
            size_bloc = high_lim - low_lim + 1
            pos_in_bloc = pos_a - sum_size_previous_blocs
            is_good_bloc = (pos_a >= low_lim and pos_a <= high_lim) or ((pos_a == 0) and (i ==0) )
            # test if pos fall into a duplicated bloc ...
            if is_good_bloc:
                is_on_edge = False
                if pos_a == high_lim:
                    print "output complex"
                    output["safe_split"] = False
                    cont_up_id = contig_a.old_2_new[ele]['contig_id']
                    up_bloc_id = contig_a.old_2_new[ele]['position']
                    cont_up = self.current_contigs[cont_up_id]
                    print "pos in bloc = ", pos_in_bloc
                    out_up = self.pos_in_bloc_2_global_pos(cont_up,pos_in_bloc,up_bloc_id)

                    output["upstream"]["bloc"] = ele
                    output["upstream"]["contig"] = cont_up
                    output["upstream"]["cut_info"] = out_up
                    output["upstream"]["exists"] = True
                    if i == len(sorted_blocs) - 1: # l index dnas les sorted blocs va de 0 a len -1
                        output["downstream"]["exists"] = False
                        output["downstream"]["bloc"] = "none"
                        output["downstream"]["contig"] = "none"
                        output["downstream"]["cut_info"] = "none"
                    else:
                        #find downstream
                        down_bloc = sorted_blocs[i + 1]
                        output["downstream"]["exists"] = True
                        output["downstream"]["bloc"] = down_bloc
                        cont_down_id = contig_a.old_2_new[down_bloc]['contig_id']
                        down_bloc_id = contig_a.old_2_new[down_bloc]['position']
                        cont_down = self.current_contigs[cont_down_id]
#                        pos_in_bloc_downstream = (pos_a + 1)  - (sum_size_previous_blocs + size_bloc)
                        pos_in_bloc_downstream = 0
                        out_down = self.pos_in_bloc_2_global_pos(cont_down,pos_in_bloc_downstream,down_bloc_id)

                        output["downstream"]["bloc"] = down_bloc
                        output["downstream"]["contig"] = cont_down
                        output["downstream"]["cut_info"] = out_down
                elif pos_a == 0:
                    print "output complex"
                    output["safe_split"] = False
                    cont_down_id = contig_a.old_2_new[ele]['contig_id']
                    down_bloc_id = contig_a.old_2_new[ele]['position']
                    cont_down = self.current_contigs[cont_down_id]
                    print "pos in bloc = ",pos_in_bloc
                    out_down = self.pos_in_bloc_2_global_pos(cont_down,pos_in_bloc,down_bloc_id)

                    output["downstream"]["bloc"] = ele
                    output["downstream"]["contig"] = cont_down
                    output["downstream"]["cut_info"] = out_down
                    output["downstream"]["exists"] = True

                    output["upstream"]["exists"] = False
                    output["upstream"]["bloc"] = "none"
                    output["upstream"]["contig"] = "none"
                    output["upstream"]["cut_info"] = "none"

                else:
                    print "output classic"
                    data_bloc = contig_a.old_2_new[ele]
                    bloc_id = data_bloc['position']
                    contig_id = data_bloc['contig_id']
                    contig_to_split = self.current_contigs[contig_id]
                    output["classic"]["contig"] = contig_to_split
                    print "ele"
                    output["classic"]["bloc"] = ele
                    status = contig_to_split.status
#                    if status == 'active':
                    print " successful candidate to split = ", contig_to_split.name
                    out_classic = self.pos_in_bloc_2_global_pos(contig_to_split,pos_in_bloc,bloc_id)
                    output["classic"]["cut_info"] = out_classic
            i +=1
            sum_size_previous_blocs += size_bloc
            size_previous_blocs.append(size_bloc)
        return output


    def post_split_update(self,contig_a, contig_1, contig_2):

        # print "inactivate contig = ",contig_a.contig_id
        self.inactivate(contig_a)
        id_1 = self.add(contig_1)
        id_2 = self.add(contig_2)

        for ele in self.current_contigs[id_1].new_compo:

            if (ele['dup']): # update duplicator and repeat object
                position = ele["position"]
                orientation = ele["orientation"]
                id_instance = ele["contig_id"]
                name_repeat = self.current_contigs[id_instance].name
                repeat = self.duplicator.duplicated_blocs[name_repeat]
                new_pos = dict()
                new_pos['orientation'] = orientation
                new_pos['position'] = position
                new_pos['@'] = id_1
                repeat.update_current_position(id_instance,new_pos)

            start = ele['start']
            end = ele['end']
            pos = str(start)+':'+str(end)

            self.current_contigs[ele['contig_id']].old_2_new[pos] = {'contig_id':id_1,'position':ele['position'],'orientation':ele['orientation']}
            tmp_n_frags = self.current_contigs[ele['contig_id']].n_frags
            pos_all = '1:'+str(tmp_n_frags)

        for ele in self.current_contigs[id_2].new_compo:

            if (ele['dup']): # update duplicator and repeat object
                position = ele["position"]
                orientation = ele["orientation"]
                id_instance = ele["contig_id"]
                name_repeat = self.current_contigs[id_instance].name
                repeat = self.duplicator.duplicated_blocs[name_repeat]
                new_pos = dict()
                new_pos['orientation'] = orientation
                new_pos['position'] = position
                new_pos['@'] = id_2
                repeat.update_current_position(id_instance,new_pos)

            start = ele['start']
            end = ele['end']
            pos = str(start)+':'+str(end)
            self.current_contigs[ele['contig_id']].old_2_new[pos] = {'contig_id':id_2,'position':ele['position'],'orientation':ele['orientation']}
            tmp_n_frags = self.current_contigs[ele['contig_id']].n_frags
            pos_all = '1:'+str(tmp_n_frags)

        return id_1, id_2


    def genome_2_abs_frag_vect(self):

        abs_frag_vect = []
        list_key = self.current_contigs.keys()
        list_key.sort()
        for contig in list_key:
            if self.current_contigs[contig].status == 'active':
                # be carefull with duplicated elements
                abs_frag_vect.extend(self.current_contigs[contig].export_abs_frag_vect())
        return abs_frag_vect

    def export_kb_pos_and_ctg_id(self):
        kb_pos_and_ctg_id = []
        list_key = self.current_contigs.keys()
        list_key.sort()
        for contig in list_key:
            if self.current_contigs[contig].status == 'active':
                kb_pos_and_ctg_id.extend(self.current_contigs[contig].export_kb_pos_and_ctg_id())
        return kb_pos_and_ctg_id

    def extract_info_frag(self,fragment_name):
        import numpy as np
#        print fragment_name
        data = fragment_name.split('-')
        id_frag = np.int32(data[0])
        contig_id =np.int32(int(data[1]))
        pos_kb = np.int32(int(float(data[2]) ))
        return id_frag, contig_id,pos_kb

    def export_pos_frag(self,):
        "deprecated"
        import numpy as np
        ind = -1
        pos_frag = []
        pos_frag.extend(self.np_id_2_frag.keys())
        pos_frag.sort()
        index_frag = np.zeros(len(pos_frag),dtype=np.int32)
        frag_multi_pos = []
        for np_frag in pos_frag:
            init_frag = self.np_id_2_frag[np_frag]
            multi_frags = self.g_sigma[init_frag]
            for ele in multi_frags:
                id_frag,contig_id,pos_kb = self.extract_info_frag(ele)
                frag_multi_pos.extend([int(pos_kb),int(contig_id)])
                ind +=2
        #            print np_frag
            index_frag[np_frag-1] = ind
        frag_multi_pos_np = np.array(frag_multi_pos,dtype=np.int32)

        return index_frag,frag_multi_pos_np

    def export_dico_frags_per_contig(self,):
        import numpy as np
        output_dico = dict()
        for contg in self.current_contigs.keys():
            tmp_ctg = self.current_contigs[contg]
            tmp_frag_np = []
            for frag in tmp_ctg.FragmentS:
                tmp_frag_np.append(frag.np_id_abs - 1)
            output_dico[contg] = np.array(tmp_frag_np,dtype = np.int32)
        return output_dico

    def export_index_matrix(self,):
        out_np_index = []
        out_np_ori = []
        out_dup = []
        n_frags = len(self.np_id_2_frag)
        rev_np_index = np.zeros((n_frags), dtype=np.int32)
        out_np_index = np.zeros((n_frags), dtype=np.int32)
        new_id = 0
        for ele in self.current_contigs:
            cont = self.current_contigs[ele]
            if cont.status == "active":
                k = 0
                for ele in cont.FragmentS:
                    np_id = ele.np_id_abs - 1
                    orientation = ele.orientation
                    bloc = self.return_bloc(cont, k + 1)
                    init_contig = ele.init_contig
                    # print "init contig = ",init_contig
                    is_duplicated = self.init_contigs[init_contig][bloc]['dup']
                    if is_duplicated:
                        out_dup.append(1)
                    else:
                        out_dup.append(0)
                    out_np_index[new_id] = np_id
                    out_np_ori.append(orientation)
                    rev_np_index[np_id] = new_id
                    k += 1
                    new_id += 1

        return out_np_index, rev_np_index, out_np_ori

    def clean_current_contigs(self,):
        list_contigs = self.current_contigs.keys()
        for c in list_contigs:
            if self.current_contigs[c].status =='inactive':
                self.current_contigs.pop(c)

    def export_pos_frag_2_np_cl(self,):
        import numpy as np
        from contig import contig as contig

        dico_contigs = dict()
        length_contigs = []
        rev_dict = dict()
        list_contigs = self.current_contigs.keys()
        # list_contigs.sort()
        for c in list_contigs:
            cont = self.current_contigs[c]
            if isinstance(cont,contig):
                if cont.status == 'active':
                    length_contigs.append(cont.n_frags)
                    dico_contigs[cont.contig_id] = dict()
                    dico_contigs[cont.contig_id]['length'] = cont.n_frags
                    # rel pos = position relative du contig dans la liste des contigs actifs
                    rel_pos = len(length_contigs) - 1
                    dico_contigs[cont.contig_id]['rel_pos'] = rel_pos
                    rev_dict[rel_pos] = cont.contig_id
                # else:
                #     self.current_contigs.pop(c)


        tmp_cumulu = ( (length_contigs[0] +1)*length_contigs[0])/2
        dico_contigs[rev_dict[0]]['cumul_index'] = 0
        dico_contigs[rev_dict[0]]['cumul_index_max'] = tmp_cumulu - 1
        prec = tmp_cumulu
        for i in xrange(1,len(length_contigs)):
            dico_contigs[rev_dict[i]]['cumul_index'] = prec

            cumulu = ((length_contigs[i] +1)*length_contigs[i]) / 2
            tmp_cumulu = prec + cumulu
            dico_contigs[rev_dict[i]]['cumul_index_max'] =tmp_cumulu -1
            prec = tmp_cumulu
        # on considere que les virtua contigs sont ordonnes ...


        pos_frag = []
        pos_frag.extend(self.np_id_2_frag.keys())
        pos_frag.sort()
        frag_coord_pos = []
        frag_list = []
        # construction dtype collect alpha contigs to do later!!!
        alpha_collect = []
        for i in xrange(0,len(length_contigs)):
            cumulu = np.int32(((length_contigs[i] +1)*length_contigs[i])) / np.int32(2)
            alpha_collect.append(cumulu)
        id = 0
        for np_frag in pos_frag:
            init_frag = self.np_id_2_frag[np_frag]
            multi_frags = self.g_sigma[init_frag]
            n_multi_frag = len(multi_frags)
            pos_data_0 = id
            for ele in multi_frags:
                id += 1
                curr_id_frag, contig_id, pos_kb = self.extract_info_frag(ele)
                low = self.current_contigs[contig_id].FragmentS[0].np_id_abs - 1
                high = self.current_contigs[contig_id].FragmentS[-1].np_id_abs - 1
                id_frag_rel = curr_id_frag - 1
                circular = 0
                size_contig = dico_contigs[contig_id]['length']
#                print "size contig = ",size_contig
                cumul_index = dico_contigs[contig_id]['cumul_index']
                alpha = 0
                frag_coord_pos.append((np.int32(contig_id),
                                       np.float32(pos_kb),
                                       np.int32(low),
                                       np.int32(high),
                                       np.int32(circular),
                                       np.int32(id_frag_rel),
                                       np.int32(size_contig),
                                       np.int32(cumul_index),
                                       np.float32(alpha)))
            pos_data_1 = id - 1
            frag_list.append((np.int32(n_multi_frag), np.int32(pos_data_0), np.int32(pos_data_1),np.int32(0)))

        frag_multi_pos_np = np.array(frag_coord_pos,dtype=self.coord_T)
        frag_list_np = np.array(frag_list,dtype = self.frag_T)


        alpha_collect = np.array(alpha_collect)
        n_alpha = alpha_collect.sum()
        alpha_collect_np = np.zeros( (n_alpha,), dtype=self.alpha_collect_T)
        # print "n element alpha collect = ", alpha_collect.sum()
        # print "size np alpha collect = ", alpha_collect_np.shape
        # print "n alpha = ", n_alpha
        # for ele in range(0,n_alpha):
        #     alpha_collect_np[ele]['obsd'] = np.float64(0)
        #     alpha_collect_np[ele]['expd'] = np.float64(0)
        return frag_list_np, frag_multi_pos_np, alpha_collect_np, dico_contigs



    def candidate_deletion(self,contig_a):
        # modifiable : augmente duplicated blocs en rajoutant les combinaison possible de bloc existant !!!
        # homogeneite : meme contig initial et meme niveau de multiplication
        # version actuelle: destruction des blocs homogenes de contigs
        # pour la mutation inverse il suffit d aller dans le dictionnaire des contigs initiaux et de prendre un candidat para ceux existant !!!
        print "scan contig to find possible deletion"
        pos_duplication = []
        FragS = contig_a.FragmentS
        for frag in FragS:
            data = self.g_sigma[frag.init_name]
            n_repeats = data.__len__()
            init_contig = frag.init_contig
            id_init = frag.id_init
            orientation = frag.orientation
            pos_duplication.append([n_repeats, init_contig, id_init,orientation])
        tmp = []
        duplicated_blocs = []
        prec_n_repeats = pos_duplication[0][0]
        prec_init_contig = pos_duplication[0][1]
        prec_id_init = pos_duplication[0][2]
        prec_orientation = pos_duplication[0][3]
        for i in range(0,len(pos_duplication)):
            n_repeats = pos_duplication[i][0]
            init_contig =  pos_duplication[i][1]
            id_init =  pos_duplication[i][2]
            orientation = pos_duplication[i][3]

            test1 = n_repeats>1
            test2 = init_contig == prec_init_contig
            test3 = orientation == prec_orientation
            test4 = orientation == 'w'
            test5 = not(n_repeats == prec_n_repeats)
            if test1:
                if test2:
                    dif_id = id_init - prec_id_init
                    if (dif_id == 0) or abs(dif_id>1) or not(test3) or (test5):
                        if len(tmp)>0:
                            duplicated_blocs.append([min(tmp) + 1 , max(tmp) + 1 ])
                            tmp = []
                        else:
                            tmp.append(i)
                    elif (dif_id == -1) and test3 and test4:
                        if len(tmp)>0:
                            duplicated_blocs.append([min(tmp) + 1 , max(tmp) + 1 ])
                            tmp = []
                        else:
                            tmp.append(i)
                    elif (dif_id) == 1 and test3 and not(test4):
                        if len(tmp)>0:
                            duplicated_blocs.append([min(tmp) + 1 , max(tmp) + 1 ])
                            tmp = []
                        else:
                            tmp.append(i)
                    else:
                        tmp.append(i)
                else:
                    if len(tmp)>0:
                        duplicated_blocs.append([min(tmp) + 1 , max(tmp) + 1 ])
                        tmp = []
                    else:
                        tmp.append(i)
            else:
                if len(tmp)>0:
                    duplicated_blocs.append([min(tmp) + 1 , max(tmp) + 1 ])
                    tmp = []

            prec_n_repeats = n_repeats
            prec_id_init = id_init
            prec_init_contig = init_contig
            prec_id_init =  id_init
            prec_orientation = orientation
        return duplicated_blocs

########################################################################################################################
########################################################################################################################
    def opposite_or(self,orientation):
        if orientation == 'downstream':
            return 'upstream'
        else:
            return 'downstream'
    def is_new_compo_contig(self,contig):
        # test contig depth
        id = contig.contig_id
        new_compo = contig.new_compo
        if id == new_compo[0]["contig_id"]:
            output = True
        else:
            output = False
        return output
    def find_coord_in_new_contig(self,contig,pos):
        # return the corresponding new contig and updated position
        output = dict()
        if self.is_new_compo_contig(contig): # if new contig no problem
            contig_id = contig.contig_id
            new_pos = pos
        else:
            new_compo = contig.new_compo
            local_bound = 1
            for ele in new_compo: # go through the new compo of the contig to retrieve the actual position
                tmp_contig_id = ele['contig_id']
                start = ele['start']
                end = ele['end']
                orientation = ele['orientation']
                n_frags = end - start + 1
                if pos >= local_bound and pos <= local_bound + n_frags - 1: # retrieve the corresponding bloc
                    bloc = str(start)+':'+str(end)
                    contig_id = tmp_contig_id
                    # retrive the position in the corresponding bloc !!
                    pos_in_bloc = pos - ( local_bound - 1 )
                    if orientation == 'w':
                        index = range(start,end+1)
                    else:
                        index = range(end,start-1,-1)
                    new_pos = index[pos_in_bloc -1 ]

                local_bound = local_bound + n_frags

        output['contig'] = self.current_contigs[contig_id]
        output['pos'] = new_pos
        # print "new contig = ",self.current_contigs[contig_id].name
        # print "pos in new contig = ",new_pos
        return output


    def return_bloc(self,contig,pos):
        init_contig = contig.FragmentS[ pos - 1 ].init_contig
        init_id = contig.FragmentS[ pos - 1 ].id_init
        for ele in self.init_contigs[init_contig]:
            if not( ele == 'n_frags'):
                tmp = ele.split(':')
                if init_id>=int(tmp[0]) and init_id<=int(tmp[1]):
                    bloc = ele
        return bloc


    def pixel_n_or_2_contig(self,pixel, orientation):
        # pixel = position in the numpy array
        # orientation = upstream or downstream orientation with respect to the init contig
        # output = position in new contig
        # attention! we need position in the  new contig

        init_frag = self.np_id_2_frag[pixel + 1] # +1 because np array start at 0 while frag dict start at 1

        set_frag = self.g_sigma[init_frag]
        print "init frag = ", init_frag
        # retrieve init frag
        frag = set_frag.pop()
        set_frag.add(frag)
        # retrive actual position of the frag
        pos_tmp,contig_id,pos_kb = self.extract_info_frag(frag)
        contig_tmp = self.current_contigs[contig_id]
        print "name of contig = ", contig_tmp.name
        print "pos in contig =", pos_tmp
        out_pos = self.find_coord_in_new_contig(contig_tmp,pos_tmp)
        pos = out_pos['pos']
        contig = out_pos['contig']
        # retrouver position precedente dans les contig de 1er niveau
        curr_id = contig.FragmentS[pos - 1].curr_id

        ## take into account orientation in contig ...
        if orientation == 'downstream':
            if curr_id == 1:
#                pos_out = 1
                pos_out = 0
            else:
#                pos_out = pos - 1
                pos_out = pos - 1
        else:
            if curr_id == contig.n_frags:
                pos_out = pos
            else:
#                pos_out = pos + 1
                pos_out = pos
        if pos_out > 0:
            init_contig = contig.FragmentS[pos_out - 1].init_contig
        else:
            init_contig = contig.FragmentS[pos_out].init_contig
        bloc = self.return_bloc(contig,pos)
        print " curr id = ",curr_id
        print " pos out -1 = ", pos_out - 1
        print " init contig = ", init_contig
        is_duplicated = self.init_contigs[init_contig][bloc]['dup']
        contig_out_id = contig.contig_id
        return pos_out, contig_out_id, is_duplicated, init_contig, bloc

    def return_correct_orientation(self,expected_or,current_or):
        if expected_or == current_or:
            output = 'w'
        else:
            output = 'c'
        return output

    def find_bloc_ending_with(self,contig,pos_end):
        blocs = contig.old_2_new
        print blocs.keys()
        print "pos end = ",pos_end
        for ele in blocs:
            data_tmp = ele.split(':')
            if int(data_tmp[1]) == pos_end:
                good_bloc = ele
        print "good bloc = ", good_bloc
        return good_bloc


    def find_bloc_starting_with(self,contig,pos_start):
        blocs = contig.old_2_new.keys()
        print blocs
        print "pos start = ",pos_start
        for ele in blocs:
            data_tmp = ele.split(':')
            if int(data_tmp[0]) == pos_start:
                good_bloc = ele
        print "good bloc = ", good_bloc
        return good_bloc


    def update_output_split(self,output_split):
        """update the dictionary yielded by split"""

        output = dict()
        output["upstream"] = dict()
        output["downstream"] = dict()
        print output_split
        if output_split["upstream"]["exists"]:
            old_up_contig  = output_split["upstream"]["contig"]
            old_up_orientation = output_split["upstream"]["orientation"]
            if old_up_contig.status == "active":
                new_up_contig = old_up_contig
                new_up_orientation = old_up_orientation
            else:
                print " find where is old up contig"
                sorted_up_bloc = old_up_contig.new_compo
                if old_up_orientation == 'w':
                    up_good_index = -1
                else:
                    up_good_index = 0
                up_bloc = sorted_up_bloc[up_good_index]
                # find now where is the up_bloc now
                up_start = up_bloc["start"]
                up_end = up_bloc["end"]
                up_string_bloc = str(up_start)+':'+str(up_end)
                up_contig_id = up_bloc["contig_id"]
                up_contig = self.current_contigs[up_contig_id]
                up_good_bloc = self.find_bloc_ending_with(up_contig, up_end)
                print "info bloc in old 2 new = ", up_contig.old_2_new[up_good_bloc]
                info_bloc_up = up_contig.old_2_new[up_good_bloc]
                new_up_contig = self.current_contigs[info_bloc_up["contig_id"]]
                new_up_orientation = info_bloc_up["orientation"]
        else:
            new_up_contig  = output_split["upstream"]["contig"]
            new_up_orientation = output_split["upstream"]["orientation"]


        output["upstream"]["exists"] = output_split["upstream"]["exists"]
        output["upstream"]["contig"] = new_up_contig
        output["upstream"]["orientation"] = new_up_orientation

        if output_split["downstream"]["exists"]:
            old_down_contig  = output_split["downstream"]["contig"]
            old_down_orientation = output_split["downstream"]["orientation"]
            if old_down_contig.status == "active":
                new_down_contig = old_down_contig
                new_down_orientation = old_down_orientation
            else:
                print " find where is old down contig"
                print "new compo downstream = ", old_down_contig.new_compo
                sorted_down_bloc = old_down_contig.new_compo
                if old_down_orientation == 'w':
                    down_good_index = 0
                else:
                    down_good_index = -1
                down_bloc = sorted_down_bloc[down_good_index]
                # find now where is the up_bloc now
                down_start = down_bloc["start"]
                down_end = down_bloc["end"]
                down_string_bloc = str(down_start)+':'+str(down_end)
                down_contig_id = down_bloc["contig_id"]
                down_contig = self.current_contigs[down_contig_id]
                down_good_bloc = self.find_bloc_starting_with(down_contig, down_start)
                print "name new contig = ",down_contig.name
                print "info bloc in old 2 new = ", down_contig.old_2_new[down_good_bloc]
                info_bloc_down = down_contig.old_2_new[down_good_bloc]
                new_down_contig = self.current_contigs[info_bloc_down["contig_id"]]
                new_down_orientation = info_bloc_down["orientation"]

        else:
            new_down_contig  = output_split["downstream"]["contig"]
            new_down_orientation = output_split["downstream"]["orientation"]

        output["downstream"]["exists"] = output_split["downstream"]["exists"]
        output["downstream"]["contig"] = new_down_contig

        output["downstream"]["orientation"] = new_down_orientation

        return output


    def  nSplit(self, samplerSplit):
        """ perform split from the sampler """
        from contig import contig as contig
        update_pos = self.return_pos_split(samplerSplit)
        proba_split = samplerSplit.proba
        output_split = dict()
        output_split["upstream"] = ""
        output_split["downstream"] = ""
        output_split["performed"] = False
        rev_split = []
        if update_pos["is_duplicated"]:
            print "duplicated blocs will be handled later"
        else:
            if update_pos["do_split"]:
                output_split["performed"] = True
                contig_to_split = update_pos["curr_contig"]
                pos_split = update_pos["pos_split"]
                upstream_tmp, downstream_tmp = contig.split(contig_to_split, pos_split)
                id_upstream, id_downstream = self.post_split_update(contig_to_split, upstream_tmp, downstream_tmp)

                upstream = self.current_contigs[id_upstream]
                downstream = self.current_contigs[id_downstream]
                rev_split_m = self.compute_reverse_split(upstream, downstream)
                # print "pos split = ", pos_split
                # print "contig to split = ",contig_to_split.name
                # print "output upstream =", upstream.name
                # print "output_ downstream = ",downstream.name
                rev_split.append(rev_split_m)
        return output_split, rev_split

    def return_pos_split(self, samplerSplit):
        """
        :param samplerSplit:
        return: output : caracteristic mutation
        """
        output = dict()
        pixel_split = samplerSplit.pixel_split
        or_split = samplerSplit.orientation
        init_frag = self.np_id_2_frag[pixel_split + 1] # +1 because np array start at 0 while frag dict start at 1
        set_frag = self.g_sigma[init_frag]
        frag = set_frag.pop()
        set_frag.add(frag)
        # retrieve actual position of the frag
        pos_tmp, contig_id, pos_kb = self.extract_info_frag(frag)
        # print "pos tmp  = ",pos_tmp
        # print "pixel_split =",pixel_split
        contig_split = self.current_contigs[contig_id]
        local_or_split = contig_split.FragmentS[pos_tmp - 1].orientation
        init_contig = contig_split.FragmentS[pos_tmp - 1].init_contig
        bloc = self.return_bloc(contig_split, pos_tmp)
        is_duplicated = self.init_contigs[init_contig][bloc]["dup"]
        nfrags = contig_split.n_frags
        pos_split = -10

        if nfrags == 1:
            do_split = False
        else:
            if pos_tmp == 1:  # frag is located @ beginning of contig
                if or_split == "upstream":
                    if local_or_split == "w":
                        do_split = False
                    else:
                        do_split = True
                        pos_split = pos_tmp
                else:
                    if local_or_split == "w":
                        do_split = True
                        pos_split = pos_tmp
                    else:
                        do_split = False

            elif pos_tmp == nfrags:  # frag is located @ end of contig
                if or_split == "upstream":
                    if local_or_split == "w":
                        do_split = True
                        pos_split = pos_tmp - 1
                    else:
                        do_split = False
                else:
                    if local_or_split == "w":
                        do_split = False
                    else:
                        do_split = True
                        pos_split = pos_tmp - 1
            else:
                if or_split == "upstream":
                    if local_or_split == "w":
                        do_split = True
                        pos_split = pos_tmp - 1
                    else:
                        do_split = True
                        pos_split  = pos_tmp

                else:
                    if local_or_split == "w":
                        do_split = True
                        pos_split = pos_tmp

                    else:
                        do_split = True
                        pos_split = pos_tmp - 1


        output["do_split"] = do_split
        output["pixel_split"] = pixel_split
        output["or_split"] = or_split
        output["pos_in_curr_contig"] = pos_tmp
        output["curr_contig"] = contig_split
        output["local_or"] = local_or_split
        output["is_duplicated"] = is_duplicated
        output["init_contig"] = init_contig
        output["bloc"] = bloc
        output["pos_split"] = pos_split

        return output

    def return_pos_paste_and_split(self, samplerPaste):

        """
        :param samplerPaste: mutation paste
        :return: caracteristic of mutation
        """

        pixel_left = samplerPaste.pixel_left
        pixel_right = samplerPaste.pixel_right
        or_left = samplerPaste.or_left
        or_right = samplerPaste.or_right
        #############################################################################
        init_frag_left = self.np_id_2_frag[pixel_left + 1] # +1 because np array start at 0 while frag dict start at 1
        set_frag_left = self.g_sigma[init_frag_left]
        # print "init frag = ", init_frag_left
        # retrieve init frag
        frag_left = set_frag_left.pop()
        set_frag_left.add(frag_left)
        # retrieve actual position of the left frag
        pos_tmp_left, contig_id_left, pos_kb_left = self.extract_info_frag(frag_left)
        contig_tmp_left = self.current_contigs[contig_id_left]

        local_or_left = contig_tmp_left.FragmentS[pos_tmp_left - 1].orientation
        init_contig_left = contig_tmp_left.FragmentS[pos_tmp_left - 1].init_contig
        bloc_left = self.return_bloc(contig_tmp_left, pos_tmp_left)
        is_duplicated_left = self.init_contigs[init_contig_left][bloc_left]["dup"]

        pos_left = pos_tmp_left
        contig_left = contig_tmp_left
        # find previous position in 1st level contigs
        curr_id = contig_left.FragmentS[pos_left - 1].curr_id
        if or_left == "upstream":
            if local_or_left == "w":
                pos_split_left = pos_left
                flip_left = False
                take_left = "left_post_split"
            else:
                pos_split_left = pos_left - 1
                flip_left = True
                take_left = "right_post_split"
        else:
            if local_or_left == "w":
                pos_split_left = pos_left - 1
                flip_left = True
                take_left = "right_post_split"
            else:
                pos_split_left = pos_left
                flip_left = False
                take_left = "left_post_split"

        #############################################################################
        init_frag_right = self.np_id_2_frag[pixel_right + 1] # +1 because np array start at 0 while frag dict start at 1
        set_frag_right = self.g_sigma[init_frag_right]
        # print "init frag = ", init_frag_right
        # retrieve init frag
        frag_right = set_frag_right.pop()
        set_frag_right.add(frag_right)
        # retrieve actual position of the right frag
        pos_tmp_right, contig_id_right, pos_kb_right = self.extract_info_frag(frag_right)
        contig_tmp_right = self.current_contigs[contig_id_right]
        local_or_right = contig_tmp_right.FragmentS[pos_tmp_right - 1].orientation
        init_contig_right = contig_tmp_right.FragmentS[pos_tmp_right - 1].init_contig
        bloc_right = self.return_bloc(contig_tmp_right, pos_tmp_right)
        is_duplicated_right = self.init_contigs[init_contig_right][bloc_right]["dup"]

        pos_right = pos_tmp_right
        contig_right = contig_tmp_right
        # find previous position in 1st level contigs
        curr_id = contig_right.FragmentS[pos_right - 1].curr_id
        if or_right == "upstream":
            if local_or_right == "w":
                pos_split_right = pos_right
                flip_right = True
                take_right = "left_post_split"
            else:
                pos_split_right = pos_right - 1
                flip_right = False
                take_right = "right_post_split"
        else:
            if local_or_right == "w":
                pos_split_right = pos_right - 1
                flip_right = False
                take_right = "right_post_split"
            else:
                pos_split_right = pos_right
                flip_right = True
                take_right = "left_post_split"

        out = dict()
        out["or_left"] = or_left
        out["or_right"] = or_right
        out["local_or_left"] = local_or_left
        out["local_or_right"] = local_or_right
        out["pos_split_left"] = pos_split_left
        out["pixel_left"] = pixel_left
        out["flip_left"] = flip_left
        out["take_left"] = take_left
        out["init_contig_left"] = init_contig_left
        out["bloc_left"] = bloc_left
        out["is_duplicated_left"] = is_duplicated_left
        out["new_contig_left"] = contig_left
        out["curr_contig_left"] = contig_tmp_left
        out["pos_in_current_left"] = pos_tmp_left
        out["pos_split_right"] = pos_split_right
        out["pixel_right"] = pixel_right
        out["flip_right"] = flip_right
        out["take_right"] = take_right
        out["init_contig_right"] = init_contig_right
        out["bloc_right"] = bloc_right
        out["is_duplicated_right"] = is_duplicated_right
        out["new_contig_right"] = contig_right
        out["curr_contig_right"] = contig_tmp_right
        out["pos_in_current_right"] = pos_tmp_right
        return out



    def nPaste(self, samplerPaste):
    #        from colorama import Fore, Back, Style
    #        from colorama import init
    #        init(autoreset=True)
        import mutation as mutation
        from contig import contig as contig
        """ new paste taking as input paste operation from the sampler """

        output = dict()
        M_rev = []
        # print " Start translocation!!!!!!!!!!!"

        out = self.return_pos_paste_and_split(samplerPaste)
        is_duplicated_left = out["is_duplicated_left"]
        is_duplicated_right = out["is_duplicated_right"]
        pixel_left = out["pixel_left"]
        pixel_right = out["pixel_right"]

        init_contig_left = out["init_contig_left"]
        init_contig_right = out["init_contig_right"]
        bloc_left = out["bloc_left"]
        bloc_right = out["bloc_right"]
        contig_left = out["curr_contig_left"]
        curr_pos_left = out["pos_in_current_left"]
        or_left = out["or_left"]
        or_right = out["or_right"]
        # print "pixel left = ", pixel_left
        # print "or left = ", or_left
        # print "pixel right = ", pixel_right
        # print "or right = ", or_right
        if is_duplicated_left or is_duplicated_right:
            """ update duplicator object"""
            print "update incoming connection"
            incoming_connection = dict()
            incoming_connection['left'] = dict()
            incoming_connection['right'] = dict()

            if is_duplicated_left:
                incoming_connection['left']['dup'] = True
                name_repeat_left = '>' + init_contig_left + ':'+bloc_left + ':w>'
                incoming_connection['left']['obj'] = name_repeat_left
            else:
                incoming_connection['left']['dup'] = False
                incoming_connection['left']['obj'] = pixel_left

            if is_duplicated_right:
                incoming_connection['right']['dup'] = True
                name_repeat_right = '>' + init_contig_right + ':' + bloc_right + ':w>'
                incoming_connection['right']['obj'] = name_repeat_right
            else:
                incoming_connection['right']['dup'] = False
                incoming_connection['right']['obj'] = pixel_right
            self.duplicator.add_modification_to_network(incoming_connection)

            output['achieved'] = False

        else:
            """ perform translocation of the selected blocs"""
            # print "perform translocation..."
            # print "contig left name = ", contig_left.name
            # print "pos cut left = ", curr_pos_left
            # print "contig right name = ", contig_right.name
            # print "pos cut right = ", curr_pos_right
            global_rev = []
            output1 = self.new_split_4_translocation(samplerPaste)

            output = self.new_split_4_translocation(samplerPaste)
            output["perform_split_left"] = output1["perform_split_left"]
            output["perform_split_right"] = output1["perform_split_right"]

            M_rev = output1["reverse_mutations"]
            out_left = output["left"]
            out_right = output["right"]
            reverse_off = False
            reverse_on = True
            l_name = out_left.to_string(reverse_off)
            r_name = out_right.to_string(reverse_off)
            r_name_rev = out_right.to_string(reverse_on)
            check_circular = l_name == r_name or l_name == r_name_rev
            if check_circular:
                print "contig circular being created"
                id_AB = self.reactivate(out_left)
                global_rev.extend(M_rev)
            else:
                # print "contig_left = ", out_left.name
                # print "contig_right = ", out_right.name
                rev_paste = self.compute_reverse_paste(out_left, output["flip_left"], out_right, output["flip_right"])

                global_rev.append(rev_paste)
                global_rev.extend(M_rev)
                self.inactivate(out_left)
                self.inactivate(out_right)
                AB = contig.join(out_left, output["flip_left"], out_right, output["flip_right"])
                id_AB = self.add(AB)
                # print "result = ", self.current_contigs[id_AB].name
            output['achieved'] = True
            output['paste'] = id_AB

        # print " End translocation!!!!!!!!!!!"
        # print " ###@@@####"
        return output, global_rev


    def new_split_4_translocation(self, samplerPaste):
        """
        :param samplerPaste: mutation paste
        :return: caracteristic of mutation
        """
        from contig import contig as contig
        M_rev = []
        #############################################################################
        pixel_left = samplerPaste.pixel_left
        or_left = samplerPaste.or_left
        init_frag_left = self.np_id_2_frag[pixel_left + 1] # +1 because np array start at 0 while frag dict start at 1
        set_frag_left = self.g_sigma[init_frag_left]
        frag_left = set_frag_left.pop()
        set_frag_left.add(frag_left)
        pos_left, contig_id_left, pos_kb_left = self.extract_info_frag(frag_left)
        contig_left = self.current_contigs[contig_id_left]
        local_or_left = contig_left.FragmentS[pos_left - 1].orientation
        init_contig_left = contig_left.FragmentS[pos_left - 1].init_contig
        bloc_left = self.return_bloc(contig_left, pos_left)
        is_duplicated_left = self.init_contigs[init_contig_left][bloc_left]["dup"]
        nfrags_left = contig_left.n_frags

        if nfrags_left == 1:
            candidate_left = contig_left
            perform_split_left = False
            if or_left == "upstream":
                if local_or_left == "w":
                    flip_left = False
                else:
                    flip_left = True
            else:
                if local_or_left == "w":
                    flip_left = True
                else:
                    flip_left = False
        else:
            if pos_left == 1:
                if or_left == "upstream":
                    if local_or_left == "w":
                        perform_split_left = True
                        pos_split_left = pos_left
                        flip_left = False
                        take_left = "left_post_split"
                    else:
                        perform_split_left = False
                        flip_left = True
                        candidate_left = contig_left
                else:
                    if local_or_left == "w":
                        perform_split_left = False
                        candidate_left = contig_left
                        flip_left = True
                    else:
                        perform_split_left = True
                        pos_split_left = pos_left
                        flip_left = False
                        take_left = "left_post_split"
            elif pos_left == nfrags_left:
                if or_left == "upstream":
                    if local_or_left == "w":
                        perform_split_left = False
                        flip_left = False
                        candidate_left = contig_left
                    else:
                        perform_split_left = True
                        pos_split_left = pos_left - 1
                        flip_left = True
                        take_left = "right_post_split"
                else:
                    if local_or_left == "w":
                        perform_split_left = True
                        pos_split_left = pos_left - 1
                        flip_left = True
                        take_left = "right_post_split"
                    else:
                        perform_split_left = False
                        flip_left = False
                        candidate_left = contig_left
            else:
                if or_left == "upstream":
                    if local_or_left == "w":
                        perform_split_left = True
                        pos_split_left = pos_left
                        flip_left = False
                        take_left = "left_post_split"
                    else:
                        perform_split_left = True
                        pos_split_left = pos_left - 1
                        flip_left = True
                        take_left = "right_post_split"
                else:
                    if local_or_left == "w":
                        perform_split_left = True
                        pos_split_left = pos_left - 1
                        flip_left = True
                        take_left = "right_post_split"
                    else:
                        perform_split_left = True
                        pos_split_left = pos_left
                        flip_left = False
                        take_left = "left_post_split"
            if perform_split_left:
                L_left, L_right = contig.split(contig_left, pos_split_left)
                ## compute reverse split
                rev_split_left = self.compute_reverse_split(L_left, L_right)
                M_rev.append(rev_split_left)
                id_left_post_split, id_right_post_split = self.post_split_update(contig_left,
                                                                                 L_left,
                                                                                 L_right)
                if take_left == "right_post_split":
                    candidate_left = self.current_contigs[id_right_post_split]
                else:
                    candidate_left = self.current_contigs[id_left_post_split]
        #############################################################################
        pixel_right = samplerPaste.pixel_right
        or_right = samplerPaste.or_right
        init_frag_right = self.np_id_2_frag[pixel_right + 1] # +1 because np array start at 0 while frag dict start at 1
        set_frag_right = self.g_sigma[init_frag_right]
        frag_right = set_frag_right.pop()
        set_frag_right.add(frag_right)
        pos_right, contig_id_right, pos_kb_right = self.extract_info_frag(frag_right)
        contig_right = self.current_contigs[contig_id_right]
        local_or_right = contig_right.FragmentS[pos_right - 1].orientation
        init_contig_right = contig_right.FragmentS[pos_right - 1].init_contig
        bloc_right = self.return_bloc(contig_right, pos_right)
        is_duplicated_right = self.init_contigs[init_contig_right][bloc_right]["dup"]
        nfrags_right = contig_right.n_frags

        if nfrags_right == 1:
            candidate_right = contig_right
            perform_split_right = False
            if or_right == "upstream":
                if local_or_right == "w":
                    flip_right = True
                else:
                    flip_right = False
            else:
                if local_or_right == "w":
                    flip_right = False
                else:
                    flip_right = True
        else:
            if pos_right == 1:
                if or_right == "upstream":
                    if local_or_right == "w":
                        perform_split_right = True
                        pos_split_right = pos_right
                        flip_right = True
                        take_right = "left_post_split"
                    else:
                        perform_split_right = False
                        flip_right = False
                        candidate_right = contig_right
                else:
                    if local_or_right == "w":
                        perform_split_right = False
                        candidate_right = contig_right
                        flip_right = False
                    else:
                        perform_split_right = True
                        pos_split_right = pos_right
                        flip_right = True
                        take_right = "left_post_split"
            elif pos_right == nfrags_right:
                if or_right == "upstream":
                    if local_or_right == "w":
                        perform_split_right = False
                        flip_right = True
                        candidate_right = contig_right
                    else:
                        perform_split_right = True
                        pos_split_right = pos_right - 1
                        flip_right = True
                        take_right = "right_post_split"
                else:
                    if local_or_right == "w":
                        perform_split_right = True
                        pos_split_right = pos_right - 1
                        flip_right = False
                        take_right = "right_post_split"
                    else:
                        perform_split_right = False
                        flip_right = True
                        candidate_right = contig_right
            else:
                if or_right == "upstream":
                    if local_or_right == "w":
                        perform_split_right = True
                        pos_split_right = pos_right
                        flip_right = True
                        take_right = "left_post_split"
                    else:
                        perform_split_right = True
                        pos_split_right = pos_right - 1
                        flip_right = False
                        take_right = "right_post_split"
                else:
                    if local_or_right == "w":
                        perform_split_right = True
                        pos_split_right = pos_right - 1
                        flip_right = False
                        take_right = "right_post_split"
                    else:
                        perform_split_right = True
                        pos_split_right = pos_right
                        flip_right = True
                        take_right = "left_post_split"
            if perform_split_right:
                R_left, R_right = contig.split(contig_right, pos_split_right)
                rev_split_right = self.compute_reverse_split(R_left, R_right)
                M_rev.append(rev_split_right)
                id_left_post_split, id_right_post_split = self.post_split_update(contig_right,
                                                                                 R_left,
                                                                                 R_right)
                if take_right == "right_post_split":
                    candidate_right = self.current_contigs[id_right_post_split]
                else:
                    candidate_right = self.current_contigs[id_left_post_split]

        output = dict()
        output["left"] = candidate_left
        output["flip_left"] = flip_left
        output["is_duplicated_left"] = is_duplicated_left
        output["perform_split_left"] = perform_split_left
        output["right"] = candidate_right
        output["flip_right"] = flip_right
        output["is_duplicated_right"] = is_duplicated_right
        output["perform_split_right"] = perform_split_right
        output["reverse_mutations"] = M_rev
        return output

    def compute_reverse_split(self, contig_left, contig_right):

        """

        :param contig_left:
        :param contig_right:
        :return:
        """
        frag_left = contig_left.FragmentS[-1]
        local_or_left = frag_left.orientation
        pixel_left = frag_left.np_id_abs - 1

        frag_right = contig_right.FragmentS[0]
        local_or_right = frag_right.orientation
        pixel_right = frag_right.np_id_abs - 1

        if local_or_left == "w":
            orientation_left = "upstream"
        else:
            orientation_left = "downstream"

        if local_or_right == "w":
            orientation_right = "downstream"
        else:
            orientation_right = "upstream"

        output_mut = mutation.paste(pixel_left, orientation_left, pixel_right, orientation_right,
                                    np.nan, np.nan, np.nan)

        return output_mut

    def compute_reverse_paste(self, contig_left, flip_left, contig_right, flip_right):
        """

        :param contig_left:
        :param flip_left:
        :param contig_right:
        :param flip_right:
        """
        import string
        if flip_left:
            frag_left = contig_left.FragmentS[0]
            or_left = frag_left.orientation.translate(string.maketrans('wc', 'cw'))
            pixel_left = frag_left.np_id_abs - 1
        else:
            frag_left = contig_left.FragmentS[-1]
            or_left = frag_left.orientation
            pixel_left = frag_left.np_id_abs - 1
        if flip_right:
            frag_right = contig_right.FragmentS[-1]
            or_right = frag_right.orientation.translate(string.maketrans('wc', 'cw'))
            pixel_right = frag_right.np_id_abs - 1
        else:
            frag_right = contig_right.FragmentS[-1]
            or_right = frag_left.orientation
            pixel_right = frag_right.np_id_abs - 1

        if or_left == "w":
            output_split = mutation.split(pixel_left, "downstream", np.nan)
        else:
            output_split = mutation.split(pixel_left, "upstream", np.nan)

        return output_split