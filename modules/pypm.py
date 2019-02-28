from Bio import SeqIO
import numpy as np
import math


def counting_nucleotides(dataset):
    st = dataset.upper()  # pour être sûr que toutes les nt soient en majuscule
    return st.count('A'), st.count('C'), st.count('G'), st.count('T')



def rna(s):
    if s == 'T':
        return 'U'
    else: return s

def dna_to_rna(s):
    print("".join(dna_to_rna(s)))  # La méthode "".join() transforme une liste en chaine de caractères, on indique le séparateur que l'on veut entre les " "
    return map(rna, s)  # Execute une fonction sur chaque item d'un élément itérable, ici la fonction rna sur l'item s



# Renvoie le complémentaire de la séquence
def complementaire(s):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return complement[s]

# Renvoie le complémentaire inversé de la séquence
def complementaire_reverse(s):
    return ''.join(reversed(list(map(complementaire, s))))  # reverse() inverse l'ordre d'une liste
# La méthode "".join() transforme une liste en chaine de caractères, on indique le séparateur que l'on veut entre les " "



def fib(Fn, Fn_1, k):
    return Fn + k * Fn_1

def fibonnacci(n, k):
    Fn = 1
    Fn_1 = 1
    for i in range(n - 2):
        Fn, Fn_1 = fib(Fn, Fn_1, k), Fn
    return Fn

# Autre méthode
# b = b + a
# a = b  => ça c'est faux

# # il faut faire
# c = b + a
# a = b
# b = c
# # ou / équivalent :
# a, b = b, a + b




# % GC d'une séquence
def percentage_GC(seq):
    nG = seq.count('G')
    nC = seq.count('C')
    n_total = len(seq)
    result = (nG + nC) / (1. * n_total) * 100.
    return result



# Différence entre 2 séquences
def hamming_distance(seq_1, seq_2):
    d = 0
    for x, y in zip(seq_1, seq_2):  # Permet de regrouper sous la forme d'un tuple les items de listes => [a, b] et [1, 2] donne [(a, 1), (b, 2)]
        if x != y:
            d += 1
    return d



# The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele
def mendelian(DD, Dr, rr):
    T = DD + Dr + rr

    # starting with k we count all
    DD_start = (DD / T) * ((DD - 1) / (T - 1) + Dr / (T - 1) + rr / (T - 1))
    # starting with m we count 1 with k 3/4 with m 1/2 with n
    Dr_start = (Dr / T) * (DD / (T - 1) + 0.75 * (Dr - 1) / (T - 1) + 0.5 * rr / (T - 1))
    # starting with n we count 1 with k 1/2 with m and 0 with n
    rr_start = (rr / T) * (DD / (T - 1) + 0.5 * Dr / (T - 1))

    frequency = DD_start + Dr_start + rr_start
    return frequency



def rna_to_protein(seq_RNA):  # rajouter une fonction try en vérifiant si seq_RNA est un multiple de 3
    aa_prots = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
                "UAU": "Y", "UAC": "Y", "UAA": "", "UAG": "", "UGU": "C", "UGC": "C", "UGA": "", "UGG": "W", "CUU": "L",
                "CUC": "L", "CUA": "L", "CUG": "L", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
                "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AUU": "I",
                "AUC": "I", "AUA": "I", "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAU": "N",
                "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V",
                "GUC": "V", "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A", "GAU": "D",
                "GAC": "D", "GAA": "E", "GAG": "E", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
    s = ""
    nrna = len(seq_RNA)
    for i in range(0, nrna, 3):
        aa = seq_RNA[i:i + 3]
        s += aa_prots[aa]
    return s



# Renvoie les positions d'un motif particulier dans une séquence
def find_motif(seq_DNA, motif):
    offset = 0
    stop = False
    s = ""
    while stop == False:
        n = seq_DNA.find(motif, offset)  # str.find(motif, beg=0, end=len(string)) envoie la première position ou se trouve le motif
        if n == -1:
            stop = True
        else:
            s += "%d " % (n + 1)
            offset = n + 1
    return s



def consensus(handle):
    sequences = []
    for record in SeqIO.parse(handle, 'fasta'):
         sequence = []
         for nt in record.seq:
              sequence.extend(nt)
         # sequence = [nt for nt in record.seq]   ça devrait être similaire mais tenir en une ligne
         sequences.append(sequence)
    handle.close()
    profile = [[0]*len(sequences)]*4
    profile = np.zeros((4, len(sequences[0])), dtype=np.int)
    for i,line in enumerate(sequences):
         for j, nt in enumerate(line):
              if nt == 'A':
                   profile[0][j] += 1
              elif nt == 'C':
                   profile[1][j] += 1
              elif nt == 'G':
                   profile[2][j] += 1
              elif nt == 'T':
                   profile[3][j] += 1

    consensus = ''
    for A,C,G,T in zip(profile[0],profile[1],profile[2],profile[3]):
         if A >= C and A >= G and A >= T:
              consensus += 'A'
         elif C >= A and C >= G and C >= T:
              consensus += 'C'
         elif G >= A and G >= C and G >= T:
              consensus += 'G'
         elif T >= A and T >= C and T >= G:
              consensus += 'T'



def parse_fasta(fasta):
    results = []
    strings = fasta.strip().split('>')  # str.strip() renvoie une copie du str dans laquelle des charactères précis ont été enlevé au début et à la fin de la string, par défaut la fonction enlève les espaces
    for s in strings:
        if len(s):
            parts = s.split()  # str.split(séparateur) transforme un str en list, si le séparateur n'est pas définit, ça enlève les espaces et les retours à la ligne
            k = parts[0]
            v = ''.join(parts[1:])
            results.append((k, v))  # renvoie une liste de tuple, avec le nom = k puis la séquence = v
    return results

def overlap_graph(fasta, n):
    results = []
    dna = parse_fasta(fasta)  # dna => liste de tuple, avec le nom = k puis la séquence = v
    for k1, v1 in dna:
        print("k1: ", k1, " v1 : ", v1)
        for k2, v2 in dna:
            print("k2 : ", k2, " v2 : ", v2)
            if k1 != k2 and v1.endswith(v2[:n]):  # str.endswitch(suffix[, start[, end]]) retourne Vrai si la chaine se termine par le suffise spécifié, sinon revoie Faux
                results.append((k1, k2))
    return results



def protein_weight(protein):
    aa_mass_dict = {'A': '71.03711', 'C': '103.00919', 'D': '115.02694', 'E': '129.04259', 'F': '147.06841',
                    'G': '57.02146', 'H': '137.05891', 'I': '113.08406', 'K': '128.09496', 'L': '113.08406',
                    'M': '131.04049', 'N': '114.04293', 'P': '97.05276', 'Q': '128.05858', 'R': '156.10111',
                    'S': '87.03203', 'T': '101.04768', 'V': '99.06841', 'W': '186.07931', 'Y': '163.06333'}
    molecular_weight = 0.
    for letter in protein:  # iterate through each letter and add weight of aa
        molecular_weight += float(aa_mass_dict[letter])
    return molecular_weight



def prime_numbers(number):
    primes = []
    for x in range(2, number + 1):
        primes.append(x)
    x = 2
    while(x <= int(math.sqrt(number))):
        if x in primes:
            for i in range(x*2, number+1, x):
                    if i in primes:
                        primes.remove(i)
        x = x+1
    return primes