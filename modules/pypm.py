def counting_nucleotides(dataset):
    st = dataset.upper()
    return st.count('A'), st.count('C'), st.count('G'), st.count('T')

test 2

def rna(s):
    if s == 'T':
        return 'U'
    else: return s

def dna_to_rna(s):
    print("".join(dna_to_rna(s)))
    return map(rna, s)



# Renvoie le complémentaire de la séquence
def complementaire(s):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return complement[s]

# Renvoie le complémentaire inversé de la séquence
def complementaire_reverse(s):
    return ''.join(reversed(list(map(complementaire, s))))



def fib(Fn, Fn_1, k):
    return Fn + k * Fn_1

def fibonnacci(n, k):
    Fn = 1
    Fn_1 = 1
    for i in range(n - 2):
        Fn, Fn_1 = fib(Fn, Fn_1, k), Fn
    return Fn



# % GC d'une séquence
def compute_GC(seq):
    nG = seq.count('G')
    nC = seq.count('C')
    n_total = len(seq)
    result = (nG + nC) / (1. * n_total) * 100.
    return result



# Différence entre 2 séquences
def hamming_distance(seq_1, seq_2):
    d = 0
    for x, y in zip(seq_1, seq_2): 
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



def rna_to_protein(seq_RNA):  # rajouter la fonction try en vérifiant si seq_RNA est un multiple de 3
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
        n = seq_DNA.find(motif, offset)
        if n == -1:
            stop = True
        else:
            s += "%d " % (n + 1)
            offset = n + 1
    return s



def parse_fasta(fasta):
    results = []
    strings = fasta.strip().split('>')
    for s in strings:
        if len(s):
            parts = s.split()
            k = parts[0]
            v = ''.join(parts[1:])
            results.append((k, v))
    return results

def overlap_graph(fasta, n):
    results = []
    dna = parse_fasta(fasta)
    for k1, v1 in dna:
        for k2, v2 in dna:
            if k1 != k2 and v1.endswith(v2[:n]):
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
