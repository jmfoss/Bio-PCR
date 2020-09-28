import random
import matplotlib.pyplot as plt

def main():
    nsp2_gene = get_gene(2720, 8554)

    # 5'-3' [positive, negative]
    dna = (nsp2_gene, get_reverse_complement(nsp2_gene))
    results = pcr(dna, 50, 20)
    display_results(results)


# Returns a random number between around 200
def get_fall_off(e):
    return 200 + random.randint(-e, e)


# param: start and end index of dna segment
# return: string of rna segment 3'-5'
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
    # get primers
    primers = get_primers(dna)
    # dna_strands=denaturation(dna)
    # dna_double_strands=annealing_elongation(dna_strands, primers, fall_off_rate)
    # put dna tuple into a list
    dna_double_strands = [dna]
    # run a loop for the cycles
    i = 0
    while i < num_cycles:
        i += 1
        temp_dna = []
        # loop through all the dna
        for x in dna_double_strands:
            # extend the dna i
            temp_dna.extend(annealing_elongation(denaturation([x]), primers, fall_off_rate))
        # put the dna Back in the list
        dna_double_strands = temp_dna

    return dna_double_strands


# param: a double strand dna, a tuple of 2 strings, representing 2 segments of dna from 5" to 3"
# return: a tuple of 2 strings representing the pair of primers [Forward, Reverse].
# (5" -> 3", GC content > 40%, bases btw the 2 primers: ~200, and GC clamp < 3)
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
    # use list comprehension to create a list with an element for each strand in the tuples of dna_segments
    dna_list = [strand for double_strand in dna_segments for strand in double_strand]
    return dna_list


# param: a list of single strand dna segments, each segment is from 5" to 3"
# return: a list of tuples of 2 strings (2 dna segments from 5" to 3")
# length of list is depended on successful bindings. Maximum size of 2 tuples
def annealing_elongation(strands, primers, fall_of_rate):
    results = list()
    # Finds new -/Template strand
    if annealing(strands[0], primers[1]):
        index = strands[1].index(primers[1])
        neg = elongation(strands[1], index, get_fall_off(fall_of_rate))
        results.append((strands[0], neg))

    # Finds new +/Coding strand
    if annealing(strands[1], primers[0]):
        index = strands[0].index(primers[0])
        pos = elongation(strands[0], index, get_fall_off(fall_of_rate))
        results.append((pos, strands[1]))

    return results


# Helper function for annealing_elongation
# Checks if the given strand of dna and primer bind to each other.
# Returns true or false
def annealing(single_strand_dna, primer):
    return single_strand_dna.count(get_reverse_complement(primer)) == 1


# Helper function for annealing_elongation
# Returns elongated strand based on given strand, primer_index, and fall_off
def elongation(single_strand_dna, primer_index, fall_off):
    return single_strand_dna[primer_index: primer_index + fall_off + 20]


# param: a list of double stranded dna segments
# displays: stats of results
def display_results(results):
    # Prints several statistics
    single_strand = [i[1 if len(i[0]) > len(i[1]) else 0] for i in results]
    print("Number of strands: ")
    print(len(single_strand))
    print("Shortest strand length: ")
    print(len(min(single_strand, key=len)))
    print("Longest strand length: ")
    print(len(max(single_strand, key=len)))
    print("Average strand length: ")
    print(round(average(single_strand)))

    # Breaks strands to individual codons and then gets the average GC codon count.
    codons = [i for ele in single_strand for i in ele]
    codon_freq = gc(codons)
    print("Count: \n" + str(codon_freq))
    gc_content = ((codon_freq["G"] + codon_freq["C"]) / (codon_freq["A"] + codon_freq["G"] + codon_freq["C"] + codon_freq["T"]))
    print("Average GC Count: \n" + str(round(gc_content, 3)) + "/1 of total codons")

    # Graphs DNA segments and their lengths
    single_strand_lengths=[len(i) for i in single_strand]
    plt.style.use('ggplot')
    plt.bar(single_strand, single_strand_lengths, color='green')
    plt.xlabel("DNA Strands")
    plt.ylabel("Strand Length")
    plt.title("NSP2 COVID19")
    plt.ylim([165, 220])
    plt.show()

    return 0

# param: a list of single strands
# returns: an average of single strand lengths
def average(list):
    lengths = [len(i) for i in list]
    return 0 if len(lengths) == 0 else (float(sum(lengths)) / len(lengths))


# param: a list of single strands
# returns: a list of individual codons
def gc(list):
    all_freq = {}
    for char in list:
        if char in all_freq:
            all_freq[char] +=1
        else:
            all_freq[char] = 1
    return all_freq


# Causes the program to run main when script is ran
if __name__ == '__main__':
    main()
