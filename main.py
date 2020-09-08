# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def get_gene(begin, end):
    sequence = open("sequence.fasta", "r").read()
    sequence = sequence.replace("\n", "")
    return sequence[begin - 1: end]


def pcr():
    nsp2_gene = get_gene(2720, 8554)
    c_dna = get_complimentary(nsp2_gene)
    cc_dna = nsp2_gene
    print(f"cDNA:  {c_dna}\nccDNA: {cc_dna}\n")
    
    
def get_complimentary(dna):
    dna = dna.replace("A", "X")
    dna = dna.replace("T", "A")
    dna = dna.replace("X", "T")
    dna = dna.replace("G", "X")
    dna = dna.replace("C", "G")
    dna = dna.replace("X", "C")
    return dna


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    pcr()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
