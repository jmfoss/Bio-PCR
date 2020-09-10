# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

def get_primers(dna):
    # Needs finished
    return 0

def get_gene(begin, end):
    sequence = open("sequence.fasta", "r").read()
    sequence = sequence.replace("\n", "")
    return sequence[begin - 1: end]

def pcr(dna, fall_off_rate, num_cycles):
    primers = get_primers(dna)

    return 0

def main():
    nsp2_gene = get_gene(2720, 8554)
    c_dna = get_reverse_complement(nsp2_gene)
    rcc_dna = nsp2_gene
    dna = (c_dna, rcc_dna)

    results = pcr(dna, 50, 20)
    display_results(results)

def display_results(results):
    # Needs finished
    return 0

def get_reverse_complement(dna):
    dna = dna.replace("A", "X")
    dna = dna.replace("T", "A")
    dna = dna.replace("X", "T")
    dna = dna.replace("G", "X")
    dna = dna.replace("C", "G")
    dna = dna.replace("X", "C")
    return dna[::-1]


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
