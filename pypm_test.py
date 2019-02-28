import modules.pypm as pypm

seq_DNA_1 = 'ACTTGCATATCGATAGCTATGACTATCGAGCATGCAGCGATG'
seq_DNA_2 = 'CCTTGCATATCGATAGCTATGACTATCGAGCATGCAGCGATG'
seq_RNA_1 = 'ACUUGCAUAUCGAUAGCUAUGACUAUCGAGCAUGCAGCGAUG'
motif = "TAT"
protein_1 = "TCISIAMTIEHAAM"

# print(pypm.counting_nucleotides(seq_DNA_1))
# print(" ".join(["%d" % x for x in pypm.counting_nucleotides(seq_DNA_1)]))
#
# # pypm.dna_to_rna(seq_DNA_1)
#
# # pypm.complementaire(seq_DNA_1)
# print(pypm.complementaire_reverse(seq_DNA_1))
#
# print(pypm.fibonnacci(30, 2))
#
# print(pypm.percentage_GC(seq_DNA_1))
#
# print(pypm.hamming_distance(seq_DNA_1, seq_DNA_2))
#
# print(pypm.mendelian(float(20), float(29), float(16)))
#
# # pypm.rna_to_protein(pypm.dna_to_rna(seq_DNA_1))
# print(pypm.rna_to_protein(seq_RNA_1))
#
# print(pypm.find_motif(seq_DNA_1, motif))
#
consensus, profile = pypm.consensus(open('Rosalind/00- datasets/rosalind_cons.txt', 'r'))
print(consensus)
print('A: ' + ' '.join(str(e) for e in profile[0]))
print('C: ' + ' '.join(str(e) for e in profile[1]))
print('G: ' + ' '.join(str(e) for e in profile[2]))
print('T: ' + ' '.join(str(e) for e in profile[3]))

# dataset = open('Rosalind/00- datasets/rosalind_grph.txt').read()
# for edge in pypm.overlap_graph(dataset, 3):
#     print(edge[0] + " " + edge[1])

# print(pypm.protein_weight(protein_1))

# print(pypm.prime_numbers(100))