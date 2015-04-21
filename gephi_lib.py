__author__ = 'hervemn'
def data2gexf(output_folder,fragments_list,fragments_contacts):
    from xml.dom.minidom import Document
    from collections import Counter
    import os
    doc = Document()
    ####################Create the <gexf> base element
    gexf = doc.createElement("gexf")
    gexf.setAttribute("xmlns", "http://www.gexf.net/1/2draft")
    gexf.setAttribute("version","1.1")
    doc.appendChild(gexf)

    ##################### Create the main <graph> element
    maingraph = doc.createElement("graph")
    maingraph.setAttribute("defaultedgetype" , "undirected")
    maingraph.setAttribute("mode","static")
    gexf.appendChild(maingraph)

    #################### define nodes##################################
    attributes_node = doc.createElement("attributes")
    attributes_node.setAttribute("class" , "node")
    attributes_node.setAttribute("mode","static")
    maingraph.appendChild(attributes_node)

    attribute_size = doc.createElement("attribute")
    attribute_size.setAttribute("id" , "size")
    attribute_size.setAttribute("title","size(kb)")
    attribute_size.setAttribute("type","integer")
    attributes_node.appendChild(attribute_size)

    attribute_gc = doc.createElement("attribute")
    attribute_gc.setAttribute("id" , "gc")
    attribute_gc.setAttribute("title","gc_content")
    attribute_gc.setAttribute("type","integer")
    attributes_node.appendChild(attribute_gc)

    attribute_coord = doc.createElement("attribute")
    attribute_coord.setAttribute("id" , "coord")
    attribute_coord.setAttribute("title","coordinates")
    attribute_coord.setAttribute("type","string")
    attributes_node.appendChild(attribute_coord)
    #################### define edges##################################
    attributes_edge = doc.createElement("attributes")
    attributes_edge.setAttribute("class" , "edge")
    attributes_edge.setAttribute("mode","static")
    maingraph.appendChild(attributes_edge)

    attribute_weight = doc.createElement("attribute")
    attribute_weight.setAttribute("id" , "weight")
    attribute_weight.setAttribute("title","weight")
    attribute_weight.setAttribute("type","float")
    attributes_edge.appendChild(attribute_coord)
    ########################### writing nodes #########################
    nodes_list = doc.createElement("nodes")
    maingraph.appendChild(nodes_list)
    nodes_list_from_file = open(fragments_list,'r')
    line_a = nodes_list_from_file.readline()
    id_fragment = 0
    dict_fragments = dict()
    while 1:
        line_a = nodes_list_from_file.readline()
        if not line_a:
            nodes_list_from_file.close()
            break
        id_fragment = id_fragment +1
        node_ele = doc.createElement("node")
        data = line_a.split()
        id = data[0]
        chrom = data[1]
        start_pos = data[2]
        end_pos = data[3]
        size = data[4]
        gc_content = data[5]
        dict_fragments[id+chrom] = id_fragment
        node_ele.setAttribute("id",str(id_fragment))
        node_ele.setAttribute("size",str(size))
        node_ele.setAttribute("coordinates",str(start_pos)+','+str(end_pos))
        node_ele.setAttribute("label",chrom)
        node_ele.setAttribute("gc_content",gc_content)
        nodes_list.appendChild(node_ele)

    ###########################  writing edges #########################s
    edges_list = doc.createElement("edges")
    edges_list_from_file = open(fragments_contacts,'r')
    maingraph.appendChild(edges_list)
    pool_contacts = Counter()
    print 'Pool contacts in dictionnary'
    while 1:
        line_a = edges_list_from_file.readline()
        if not line_a:
            edges_list_from_file.close()
            break
        data = line_a.split()
        source = min( int(data[0]),int(data[1]) )
        target = max( int(data[0]),int(data[1]) )
        pool_contacts[str(source)+'-'+str(target)] += 1
    i = 0
    for contact in pool_contacts.keys():
        i = i+1
        weight = pool_contacts[contact]
        partners = contact.split('-')
        edge_ele = doc.createElement("edge")
        edge_ele.setAttribute("id",str(i))
        edge_ele.setAttribute("source",partners[0])
        edge_ele.setAttribute("target",partners[1])
        edge_ele.setAttribute("weight",str(weight))
        edges_list.appendChild(edge_ele)

    output = open(os.path.join(output_folder,'graph.gexf'),'w')
    print 'write gephi file'
    doc.writexml(output,indent="\t",addindent="\t",newl="\n")
    output.close()
    print 'done.'