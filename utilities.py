def colored(seq):
    bcolors = {
        'A' : '\033[92m',
        'C' : '\033[94m',
        'G' : '\033[93m',
        'T' : '\033[91m',
        'U' : '\033[91m',
        'reset': '\033[0;0m'
    }
    rmpStr = ''

    for nuc in seq:
        if nuc in bcolors:
            tmpStr = bcolors[nuc] + nuc
        else:
            tmpStr += bcolors['reset'] + nuc

    return tmpStr + '\033[0;0m'

def fasta_dict(name):
    file = open(name,'r')
    data = file.read()
    file.close()
    data1 = data.split('\n')
    fasta_dict = {}
    for line in data1:
        if '>' in line:
            label = line
            fasta_dict[label] = ''
        else:
            fasta_dict[label] += line
    return fasta_dict