def main():
    nsp2_gene = get_gene(2720, 8554)  # RNA 3'-5'
    dna = (nsp2_gene, get_reverse_complement(nsp2_gene))  # dna both 5'-3' dna[0]: postive, dna[1]: negative
    results = pcr(dna, 50, 20)
    display_results(results)


# param: start and end index of dna segment
# return: string of dna segment 3' to 5'
def get_gene(begin, end):
    sequence = open("sequence.fasta", "r").read()
    sequence = sequence.replace("\n", "")
    return sequence[begin - 1: end]


# param: a single strand 3" to 5" rna
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
    forward = [dna[0][iterator: iterator + 20] for iterator in range(len(dna[0]) - 20)
               if check_gc_content(dna[0][iterator: iterator + 20], 8) and
               check_primer_bind(dna[0][iterator: iterator + 20], dna[1]) and
               check_gc_clamp(dna[0][iterator: iterator + 20])]

    reverse = [dna[1][iterator: iterator + 20] for iterator in range(len(dna[1]) - 20)
               if check_gc_content(dna[1][iterator: iterator + 20], 8) and
               check_primer_bind(dna[1][iterator: iterator + 20], dna[0]) and
               check_gc_clamp(dna[1][iterator: iterator + 20])]
    return [(i, j) for i in forward for j in reverse if 190 < find_product_length(dna, i, j) < 200][0]


# Helper function for get_primers
def check_primer_bind(primer, complementary):
    return complementary.count(get_reverse_complement(primer)) == 1


# Helper function for get_primers
def find_product_length(dna, forward, reverse):
    return dna[1][::-1].index(reverse[::-1]) - dna[0].index(forward) + 20


# Helper function for get_primers
def check_gc_clamp(primer):
    return 0 < len(primer[14: 20].replace('A', '').replace('T', '')) < 3


# Helper function for get_primers
def check_gc_content(primer, content):
    return len(primer.replace('A', '').replace('T', '')) > content


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
