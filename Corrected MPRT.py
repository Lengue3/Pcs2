import requests


url = 'https://www.uniprot.org/uniprot/'
 
with open("rosalind_mprt.txt") as file:
    seqIDs = file.read().replace("\n", " ").split()
sequences = {}

for proID in seqIDs:
    trim = proID.split('_')
    valid_ID = trim[0]
    goToURL = url+valid_ID+".fasta"
    response = requests.get(goToURL)
    sequences[proID] = (response.text.split("\n"))
    sequences[proID] = "".join(sequences[proID][1::])
#copied this part of the problem because I did not know how to do it, the rest is mine

for key in sequences.keys():
    seq = sequences[key]
    positions = []
    for i in range(len(seq)-2):
        test = 'AAAA'
        if seq[i] == 'N':
            flag1 = False
            try:
                test = seq[i:i+4]
                check = test[3]
            except:
                test = 'AAAA'

        if test[1] == 'P' or test[3] == 'P':
            flag1 = True
        

        if test[2] == 'S' or test[2] == 'T' and flag1 != True:
            positions.append(str(i+1))

    if len(positions) > 0:
        print(key)
        print(' '.join(positions))