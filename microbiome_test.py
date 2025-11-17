# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 17:32:48 2018

@author: jdkan
"""
import numpy
import time
import alignment
import copy
import random
import matplotlib.pyplot as plt

def Load16SFastA(path, fraction = 1.0):
    # from a file, read in all sequences and store them in a dictionary
    # sequences can be randomly ignored for testing by adjusting the fraction
    random.seed(11)

    infile = open(path, 'r')
    sequences_16s = {}
    c = 0
    my_seq = ""
    for line in infile:
        if ">" in line:
            my_id = line[1:-1]
            if random.random() < fraction:
                sequences_16s[my_id] = ""


        else:
            if my_id in sequences_16s:
                sequences_16s[my_id] += line[:-1]


    return sequences_16s

def split_dataset(seqs, lib_len, quer_len):
    items = list(seqs.items())
    random.shuffle(items)
    library = dict(items[quer_len:(quer_len + lib_len)])
    queries = dict(items[:quer_len])
    return library, queries

def ConvertLibraryToKmerSets(library, k=2):
    new_lib = {}
    c = 0
    for key in library.keys():
        new_lib[key] = set()
        seq = library[key]
        for i in range(0, len(seq) - k + 1):
            new_lib[key].add(seq[i:i+k]) # sets automatically remove duplicate values

    return new_lib

def convert_seq_to_kmer_set(seq, k=2):
    kmers = set()
    for i in range(0, len(seq) - k + 1):
        kmers.add(seq[i:i + k])
    return kmers

def JaccardIndex(s1, s2):
    numerator = float(len(s1.intersection(s2)))
    denominator = float(len(s1.union(s2)))
    return numerator/denominator

def KmerMatch(sequence_kmer_set, library_kmer_set):
    best_score = 0.0
    best_match = None

    for key in library_kmer_set:
        score = JaccardIndex(sequence_kmer_set, library_kmer_set[key])
        if score > best_score:
            best_score = score
            best_match = key

    return best_score, best_match


def AlignmentMatch(sequence, library):
    best_score = -10000000000
    best_match = None

    for x in library:
        score, loc, matrix = alignment.local_align(sequence, library[x])
        if score > best_score:
            best_score = score
            best_match = x # or do we want to return the actual sequence, library[x] ?

    return best_score, best_match

def find_agreement(queries, library):
    '''
    See whether AlignmentMatch and KmerMatch agree on the best library sequence match for
    (k=1,3,5,7,9,11,13,15,17,19)

    input:
        queries (dict of sequences)
        library (dict of sequences)
    output:
        agreements (list of fraction of match agreements for k=1,3,5,7,9,11,13,15,17,19)
    '''
    agreements = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    total = len(queries)

    # Make kmer-sets for each k
    library_k1 = ConvertLibraryToKmerSets(library, k=1)
    library_k3 = ConvertLibraryToKmerSets(library, k=3)
    library_k5 = ConvertLibraryToKmerSets(library, k=5)
    library_k7 = ConvertLibraryToKmerSets(library, k=7)
    library_k9 = ConvertLibraryToKmerSets(library, k=9)
    library_k11 = ConvertLibraryToKmerSets(library, k=11)
    library_k13 = ConvertLibraryToKmerSets(library, k=13)
    library_k15 = ConvertLibraryToKmerSets(library, k=15)
    library_k17 = ConvertLibraryToKmerSets(library, k=17)
    library_k19 = ConvertLibraryToKmerSets(library, k=19)

    queries_k1 = ConvertLibraryToKmerSets(queries, k=1)
    queries_k3 = ConvertLibraryToKmerSets(queries, k=3)
    queries_k5 = ConvertLibraryToKmerSets(queries, k=5)
    queries_k7 = ConvertLibraryToKmerSets(queries, k=7)
    queries_k9 = ConvertLibraryToKmerSets(queries, k=9)
    queries_k11 = ConvertLibraryToKmerSets(queries, k=11)
    queries_k13 = ConvertLibraryToKmerSets(queries, k=13)
    queries_k15 = ConvertLibraryToKmerSets(queries, k=15)
    queries_k17 = ConvertLibraryToKmerSets(queries, k=17)
    queries_k19 = ConvertLibraryToKmerSets(queries, k=19)

    # Make dict for best match using local alignment
    alignment_match = {}
    for i in queries:
        seq = queries[i]
        best_score, best_match = AlignmentMatch(seq, library)
        alignment_match[i] = best_match

    # For each sequence, get the best match using kmer_alignment
    for i in queries:
        score_1, match_1 = KmerMatch(queries_k1[i], library_k1)
        score_3, match_3 = KmerMatch(queries_k3[i], library_k3)
        score_5, match_5 = KmerMatch(queries_k5[i], library_k5)
        score_7, match_7 = KmerMatch(queries_k7[i], library_k7)
        score_9, match_9 = KmerMatch(queries_k9[i], library_k9)
        score_11, match_11 = KmerMatch(queries_k11[i], library_k11)
        score_13, match_13 = KmerMatch(queries_k13[i], library_k13)
        score_15, match_15 = KmerMatch(queries_k15[i], library_k15)
        score_17, match_17 = KmerMatch(queries_k17[i], library_k17)
        score_19, match_19 = KmerMatch(queries_k19[i], library_k19)

        # Update agreements based on matches
        if match_1 == alignment_match[i]:
            agreements[0] += 1
        if match_3 == alignment_match[i]:
            agreements[1] += 1
        if match_5 == alignment_match[i]:
            agreements[2] += 1
        if match_7 == alignment_match[i]:
            agreements[3] += 1
        if match_9 == alignment_match[i]:
            agreements[4] += 1
        if match_11 == alignment_match[i]:
            agreements[5] += 1
        if match_13 == alignment_match[i]:
            agreements[6] += 1
        if match_15 == alignment_match[i]:
            agreements[7] += 1
        if match_17 == alignment_match[i]:
            agreements[8] += 1
        if match_19 == alignment_match[i]:
            agreements[9] += 1

    ###
    for i in range(0, len(agreements)):
        agreements[i] /= total

    return agreements

ALL_LETTERS = ["A", "C", "G", "T"]

#Task 5: Minimizer

def get_minimizer(sequence, m, k):
    minimizers = set()
    n = len(sequence)

    win_num = n // (m - k)

    for i in range(win_num):
        windom = sequence[i*(m-k): i*(m-k) + m]

        kmer_list = []

        for j in range(m - k + 1):
            kmer = windom[j:j + k]
            kmer_list.append(kmer)

        minimizers.add(min(kmer_list))

    return minimizers

def generate_min_library(input_library, m, k):
    min_library = {}
    for key in input_library.keys():
        seq = input_library[key]
        min_library[key] = get_minimizer(seq, m, k)
    return min_library

def find_matched_keys(query_seq, min_library, m, k):
    query_minimizer = get_minimizer(query_seq, m, k)
    matched_keys = []
    for key in min_library.keys():
        lib_minimizers = min_library[key]
        if (lib_minimizers & query_minimizer):
            matched_keys.append(key)
    return matched_keys

def find_minimizer_agreement(queries, min_library, seq_lib, m, k):
    agreements = 0
    total = len(queries)

    for key in queries.keys():
        new_seq_lib = dict()
        query_seq = queries[key]
        matched_keys = find_matched_keys(query_seq, min_library, m, k)
        for matched_key in matched_keys:
            new_seq_lib[matched_key] = seq_lib[matched_key]
        _, best_match = AlignmentMatch(query_seq, seq_lib)
        seq_kmers = convert_seq_to_kmer_set(query_seq, k)
        lib_kmers = ConvertLibraryToKmerSets(new_seq_lib, k)
        _, kmer_best_match = KmerMatch(seq_kmers, lib_kmers)
        if best_match == kmer_best_match:
            agreements += 1
    
    return agreements / total

def generate_agreement(queries, library):
    win = [i for i in range(20, 66, 5)]
    agreements = []
    for w in win:
        min_library = generate_min_library(library, m=w, k=7)
        agreement = find_minimizer_agreement(queries, min_library, library, m=w, k=7)
        agreements.append(agreement)
    return agreements

def mutate_string(query, rate):
    new_query = ""
    for char in query:
        change = numpy.random.binomial(1, rate)
        if change:
            try:
                ALL_LETTERS.remove(char)
                new_query += random.choice(ALL_LETTERS)
                ALL_LETTERS.append(char)
            except:
                new_query += random.choice(ALL_LETTERS)
        else:
            new_query += char
    
    return new_query



if __name__ == "__main__":
    fn = "bacterial_16s_genes.fa"
    sequences_16s = Load16SFastA(fn, fraction = 0.015)

    print("Loaded %d 16s sequences." % len(sequences_16s))

    library, queries = split_dataset(sequences_16s, 10, 5)
    print("Library of length %d, Queries of length %d" % (len(library), len(queries)))

    agreements = find_agreement(library, queries)
    print(agreements)

    ks = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
    plt.plot(ks, agreements)
    plt.title("Agreement Curve")
    plt.xlabel("K-mer length")
    plt.ylabel("Agreement")
    plt.xticks([1, 3, 5, 7, 9, 11, 13, 15, 17, 19])
    plt.show()

    #Task 6: Agreement Curve with minimizer
    min_agreements = generate_agreement(queries, library)
    plt.plot([i for i in range(20, 66, 5)], min_agreements)
    plt.title("Minimizer Agreement Curve")
    plt.xlabel("Window size (m)")
    plt.ylabel("Agreement")
    plt.xticks([i for i in range(20, 66, 5)])
    plt.show()


    print("=============== Testing with 1% mutated queries ===============")
    library, queries = split_dataset(sequences_16s, 10, 5)
    print("Library of length %d, Queries of length %d" % (len(library), len(queries)))

    for key in queries.keys():
        queries[key] = mutate_string(queries[key], 1.0/100)

    agreements = find_agreement(library, queries)
    print(agreements)


    plt.plot(ks, agreements)
    plt.title("Agreement Curve- 1% mutated")
    plt.xlabel("K-mer length")
    plt.ylabel("Agreement")
    plt.xticks([1, 3, 5, 7, 9, 11, 13, 15, 17, 19])
    plt.show()


    print("=============== Testing with 10% mutated queries ===============")
    library, queries = split_dataset(sequences_16s, 10, 5)
    print("Library of length %d, Queries of length %d" % (len(library), len(queries)))

    for key in queries.keys():
        queries[key] = mutate_string(queries[key], 1.0/10)

    plt.plot(ks, agreements)
    plt.title("Agreement Curve- 10% mutated")
    plt.xlabel("K-mer length")
    plt.ylabel("Agreement")
    plt.xticks([1, 3, 5, 7, 9, 11, 13, 15, 17, 19])
    plt.show()


    # print(queries)

    # seq1 = queries['RS_GCF_001553625.1~NZ_CP014230.1 d__Bacteria;p__Desulfobacterota;c__Desulfovibrionia;o__Desulfovibrionales;f__Desulfomicrobiaceae;g__Desulfomicrobium;s__ 1548 2783374']
    # best_score, best_match = AlignmentMatch(seq1, library)
    # print("----------------------\n"
    #       "LOCAL ALIGNMENT\n "
    #       "Best score: %d \n"
    #       "Best match: %s" % (best_score, best_match))
    # kmer_16s_sequences = ConvertLibraryToKmerSets(sequences_16s, k=6)
    # print(convert_seq_to_kmer_set('ACTGA', 4))

    # kmers = convert_seq_to_kmer_set(seq1)
    # library_kmers = ConvertLibraryToKmerSets(library)
    # best_score, best_match = KmerMatch(kmers, library_kmers)
    # print("----------------------\n"
    #       "KMER ALIGNMENT\n "
    #       "Best score: %d \n"
    #       "Best match: %s" % (best_score, best_match))











