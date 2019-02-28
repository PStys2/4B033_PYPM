import modules.pypm as pypm

seq_DNA_1 = 'ACTTGCATATCGATAGCTATGACTATCGAGCATGCAGCGATG'
seq_DNA_2 = 'CCTTGCATATCGATAGCTATGACTATCGAGCATGCAGCGATG'
seq_RNA_1 = 'ACUUGCAUAUCGAUAGCUAUGACUAUCGAGCAUGCAGCGAUG'
motif = "TAT"
protein_1 = "TCISIAMTIEHAAM"

print(pypm.counting_nucleotides(seq_DNA_1))
print(" ".join(["%d" % x for x in pypm.counting_nucleotides(seq_DNA_1)]))

# pypm.dna_to_rna(seq_DNA_1)

# pypm.complementaire(seq_DNA_1)
print(pypm.complementaire_reverse(seq_DNA_1))

print(pypm.fibonnacci(30, 2))

print(pypm.compute_GC(seq_DNA_1))

print(pypm.hamming_distance(seq_DNA_1, seq_DNA_2))

print(pypm.mendelian(float(20), float(29), float(16)))

# pypm.rna_to_protein(pypm.dna_to_rna(seq_DNA_1))
print(pypm.rna_to_protein(seq_RNA_1))

print(pypm.find_motif(seq_DNA_1, motif))

dataset = open('rosalind_grph.txt').read()
for i in pypm.overlap_graph(dataset, 3):
    print(i[0] + " " + i[1])

print(pypm.protein_weight(protein_1))
