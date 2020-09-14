# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
test = "1234554321"
r = [(test[i], test[i]) for i in range(len(test)) if test[i] == '5']


def main():
    nsp2_gene = get_gene(2720, 8554)
    dna = (nsp2_gene, get_reverse_complement(nsp2_gene))
    results = pcr(dna, 50, 20)
    display_results(results)


# param: start and end index of dna segment
# return: string of dna segment 5" to 3"
def get_gene(begin, end):
    sequence = open("sequence.fasta", "r").read()
    sequence = sequence.replace("\n", "")
    return sequence[begin - 1: end]


# param: a single strand 5" to 3" dna
# return: a single strand 5" to 3" dna that is reverse complement to the input strand
def get_reverse_complement(dna):
    dna = dna.replace("A", "X")
    dna = dna.replace("T", "A")
    dna = dna.replace("X", "T")
    dna = dna.replace("G", "X")
    dna = dna.replace("C", "G")
    dna = dna.replace("X", "C")
    return dna[::-1]


# param: gene to be copied (a tuple of 2 strings), fall of rate of DNA polymerase (int),
# and num_cycles to run PCR (int)
# return: a list of double stranded dna segments
def pcr(dna, fall_off_rate, num_cycles):
    primers = get_primers(dna)
    print(primers)
    # Needs finished
    return 0


# param: a double strand dna, a tuple of 2 strings, representing 2 segments of dna from 5" to 3"
# return: a tuple of 2 strings representing the pair of primers (Forward, Reverse).
# (5" -> 3", GC content > 40%, bases btw the 2 primers: ~200)
def get_primers(dna):
    # Finds first pair of two primers, forward and reverse (both 5' to 3'), 200 bases apart
    # each of size of 20 bases, with GC content > 40%
    return [(dna[0][iterator: iterator + 20], dna[1][::-1][iterator + 220: iterator + 240][::-1])
            for iterator in range(len(dna[0]) - 240)
            if len(dna[0][iterator: iterator + 20].replace('A', '').replace('T', '')) > 8 and
            len(dna[1][::-1][iterator + 220: iterator + 240].replace('A', '').replace('T', '')) > 8][0]


# param: a list of tuples of 2 strings, representing double stranded dna segments
# return: a list of single strand dna segments
def denaturation(dna_segments):
    # Needs finished
    return 0


# param: a list of single strand dna segments, each segment is from 5" to 3"
# return: a list of tuples of 2 strings (2 dna segments from 5" to 3")
def annealing_elongation(single_strand_dna, primers, fall_of_rate):
    # Needs finished
    return 0


# param: a list of double stranded dna segments
# displays: stats of results
def display_results(results):
    # Needs finished
    return 0


# Causes the program to run main when script is ran
if __name__ == '__main__':
    main()
