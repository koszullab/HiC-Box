__author__ = 'hervemn'


def remove_self_fragments_contacts(input_file,output_file):
    a = open(input_file,'r')
    b = open(output_file,'w')
    for line in a:
        dat = line.split()
        if not((dat[0] == dat[2]) and (dat[1] == dat[3])):
            b.write(line)
