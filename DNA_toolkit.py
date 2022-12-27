from structures import *
from collections import Counter
def validateSeq(dna_seq):
    Nucleotides = ['A','C','G','T']
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

def countNucFrequency(seq):
    
    for nuc in seq:
        for nuc in seq:
            tmpFreqDict[nuc] += 1
    return tmpFreqDict

def transcription(seq):

    return seq.replace ('T','U')


def ReverseComp(seq):
    mapping = str.maketrans('ATCG','TAGC')
    return seq.translate(mapping)[::-1]

def gc_content(seq):

    return round((seq.count('C') + seq.count('G')) / len(seq) * 100)

def translate_seq(seq,init_pos=0):

    return [DNA_codons[seq[pos:pos + 3]]for pos in range(init_pos, len(seq)-2,3)]

def gen_ORFs(seq):
    frames = []
    frames.append(translate_seq(seq,0))
    frames.append(translate_seq(seq,1))
    frames.append(translate_seq(seq,2))
    frames.append(translate_seq(ReverseComp(seq),0))
    frames.append(translate_seq(ReverseComp(seq),1))
    frames.append(translate_seq(ReverseComp(seq),2))
    return frames

def proteins_from_rf(aa_seq):
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == '_':

            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            if aa == 'M':
                current_prot.append('')
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins


def all_proteins_from_orfs(seq, StartReadpos=0,EndReadpos=0,ordered=False):
    if EndReadpos > StartReadpos:
        rfs = gen_ORFs(seq[StartReadpos:EndReadpos])
    else:
        rfs = gen_ORFs(seq)
    
    res = []

    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)
    if ordered:
        return sorted(res,key=len, reverse=True)
    return res


def find_superstring(substrings, result = ""):
    #checks one by one each sequence, looking for matches an adding them together in the result.

    
    if len(substrings) == 0:    #termination step, when all the substrings have been removed, the function will return the result
        return result

    
    elif len(result) == 0:
        result = substrings.pop(0)        #starting step, the first substring will become our starting sequence so we can build up on it
        return find_superstring(substrings, result)

    else:
        for i in range(len(substrings)): #start iterating through all substrings
            a = substrings[i]
            l = len(a)

            for p in range(int(l/2)):   #rosalind gives us a hint, the matching will at least be 50%
                q = l - p               #so we start cutting the string in two, and trying both halfs for a match

                if result.startswith(a[p:]): #if a[p:] matches, then we must add a[:p] since we dont want to repeat both times the same region
                    substrings.pop(i)        # if the highest match possible is found, there is no need to keep analyzing the substring so we pop it
                    return find_superstring(substrings, a[:p] + result)

                if result.endswith(a[:q]): #same idea here just repeated
                    substrings.pop(i) 
                    return find_superstring(substrings, result + a[q:])

