__author__ = 'hervemn'

def int_to_roman (integer):

    returnstring=''
    table=[['M',1000],['CM',900],['D',500],['CD',400],['C',100],['XC',90],['L',50],['XL',40],['X',10],['IX',9],['V',5],['IV',4],['I',1]]

    for pair in table:

        while integer-pair[1]>=0:

            integer-=pair[1]
            returnstring+=pair[0]

    return returnstring

def rom_to_int(string):

    table=[['M',1000],['CM',900],['D',500],['CD',400],['C',100],['XC',90],['L',50],['XL',40],['X',10],['IX',9],['V',5],['IV',4],['I',1]]
    returnint=0
    for pair in table:


        continueyes=True

        while continueyes:
            if len(string)>=len(pair[0]):

                if string[0:len(pair[0])]==pair[0]:
                    returnint+=pair[1]
                    string=string[len(pair[0]):]

                else: continueyes=False
            else: continueyes=False

    return returnint


input_a = open('/Users/hervemn/Dev/bowtie-0.12.7/S288C_reference_sequence_R64-1-1_20110203.fsa','r')
output_a = open('/Users/hervemn/Dev/bowtie-0.12.7/new_ref_genome.fsa','w')

while 1:
    line_a = input_a.readline()
    if not line_a:
        input_a.close()
        output_a.close()
        break
    if line_a[0] == '>':
        data = line_a.split(' ')
        chrtmp0 = data[5]
        chrtmp1 = chrtmp0.replace('[chromosome=','')
        chrtmp1 = chrtmp1.replace('[location=','')
        chr = chrtmp1.replace(']\n','')
        chr = chr.replace(']','')
        print chr
        if chr == 'mitochondrion':
            output_a.write('>Scmito\n')
        else:
            n = rom_to_int(chr)
            chrstr = str(n)
            out_str = '>Scchr'+chrstr.zfill(2)
            output_a.write(out_str+'\n')
    else:
        output_a.write(line_a)
