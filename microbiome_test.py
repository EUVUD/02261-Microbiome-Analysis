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
from collections import defaultdict
import concurrent.futures
import matplotlib.pyplot as plt
import json

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

def load_datasets():
    with open("library_len200.json", "r") as fp:
        library = json.load(fp)

    with open("queries_len50.json", "r") as fp:
        queries = json.load(fp)

    with open("alignment_matches.json", "r") as fp:
        alignment_match = json.load(fp)

    return library, queries, alignment_match

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
    print("Finding agreements.................")

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
    print("Finished converting library to kmer sets")

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

    print("Finished converting queries to kmer sets")
    print("Starting local alignments")
    # Make dict for best match using local alignment
    alignment_match = {}
    for i in queries:
        seq = queries[i]
        best_score, best_match = AlignmentMatch(seq, library)
        alignment_match[i] = best_match
    print("Finished local alignments")
    with open('alignment_matches.json', 'w') as file:
        json.dump(alignment_match, file, indent=4)

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

def find_agreement_fast(queries, library, alignment_match):
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
    print("Finding agreements.................")

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
    print("Finished converting library to kmer sets")

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

    print("Finished converting queries to kmer sets")

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

    if not sequence or m <= 0 or k <= 0 or n < k:
        return minimizers

    step = max(1, m - k)
    win_num = n // step
    window_kmer_len = m - k + 1
    total_kmers = n - k + 1

    if window_kmer_len <= 0 or total_kmers <= 0 or win_num == 0:
        return minimizers

    kmers = [sequence[i:i + k] for i in range(total_kmers)]
    max_kmer_start = total_kmers - 1

    for i in range(win_num):
        window_start = i * step
        if window_start > max_kmer_start:
            break

        end_idx = min(window_start + window_kmer_len, total_kmers)
        if end_idx <= window_start:
            continue

        window_min = min(kmers[idx] for idx in range(window_start, end_idx))
        minimizers.add(window_min)

    return minimizers

def generate_min_library(input_library, m, k):
    min_library = {}
    minimizer_index = defaultdict(set)

    for key, seq in input_library.items():
        mins = get_minimizer(seq, m, k)
        min_library[key] = mins
        for minimizer in mins:
            minimizer_index[minimizer].add(key)

    return min_library, minimizer_index


def generate_min_library_parallel(input_library, m, k, n_workers=2):
    """Parallel version of generate_min_library using processes.

    Returns (min_library, minimizer_index)
    """
    items = list(input_library.items())
    min_library = {}
    minimizer_index = defaultdict(set)

    if n_workers is None or n_workers <= 1:
        return generate_min_library(input_library, m, k)

    def _compute_min(pair):
        key, seq = pair
        return key, get_minimizer(seq, m, k)

    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as exc:
        for key, mins in exc.map(_compute_min, items):
            min_library[key] = mins
            for minimizer in mins:
                minimizer_index[minimizer].add(key)

    return min_library, minimizer_index

def find_matched_keys(query_seq, min_library, m, k, minimizer_index=None, query_minimizers=None):
    if query_minimizers is None:
        query_minimizers = get_minimizer(query_seq, m, k)

    if not query_minimizers:
        return []

    if minimizer_index is not None:
        matched = set()
        for minimizer in query_minimizers:
            matched.update(minimizer_index.get(minimizer, ()))
        return list(matched)

    return [
        key
        for key, lib_minimizers in min_library.items()
        if lib_minimizers & query_minimizers
    ]

def find_minimizer_agreement(
    queries,
    min_library,
    seq_lib,
    m,
    k,
    minimizer_index=None,
    alignment_cache=None,
    query_kmer_cache=None,
    library_kmer_cache=None,
):
    total = len(queries)
    if total == 0:
        return 0.0

    if alignment_cache is None:
        alignment_cache = {
            key: AlignmentMatch(seq, seq_lib)[1]
            for key, seq in queries.items()
        }

    if query_kmer_cache is None:
        query_kmer_cache = {
            key: convert_seq_to_kmer_set(seq, k)
            for key, seq in queries.items()
        }

    if library_kmer_cache is None:
        library_kmer_cache = ConvertLibraryToKmerSets(seq_lib, k)

    agreements = 0

    for key, query_seq in queries.items():
        query_min = get_minimizer(query_seq, m, k)
        if not query_min:
            continue

        matched_keys = find_matched_keys(
            query_seq,
            min_library,
            m,
            k,
            minimizer_index=minimizer_index,
            query_minimizers=query_min,
        )

        if not matched_keys:
            continue

        subset_kmers = {
            mk: library_kmer_cache[mk]
            for mk in matched_keys
            if mk in library_kmer_cache
        }

        if not subset_kmers:
            continue

        _, kmer_best_match = KmerMatch(query_kmer_cache[key], subset_kmers)
        if alignment_cache.get(key) == kmer_best_match:
            agreements += 1

    return agreements / total

def generate_agreement(queries, library, n_workers=1):
    """Compute agreement for a range of window sizes. Optionally parallelize heavy work by
    specifying n_workers>1 (uses process pool).
    """
    win = [i for i in range(20, 66, 5)]
    agreements = []
    if not queries or not library:
        return agreements

    fixed_k = 7

    # Build caches (alignment and kmer sets). Alignment is CPU-heavy so parallelize when requested.
    if n_workers and n_workers > 1:
        def _align_item(item):
            key, seq = item
            return key, AlignmentMatch(seq, library)[1]

        items = list(queries.items())
        alignment_cache = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as exc:
            for key, best_match in exc.map(_align_item, items):
                alignment_cache[key] = best_match
    else:
        alignment_cache = {key: AlignmentMatch(seq, library)[1] for key, seq in queries.items()}

    query_kmer_cache = {key: convert_seq_to_kmer_set(seq, fixed_k) for key, seq in queries.items()}
    library_kmer_cache = ConvertLibraryToKmerSets(library, fixed_k)

    for w in win:
        if n_workers and n_workers > 1:
            min_library, minimizer_index = generate_min_library_parallel(library, m=w, k=fixed_k, n_workers=n_workers)
        else:
            min_library, minimizer_index = generate_min_library(library, m=w, k=fixed_k)

        agreement = find_minimizer_agreement(
            queries,
            min_library,
            library,
            m=w,
            k=fixed_k,
            minimizer_index=minimizer_index,
            alignment_cache=alignment_cache,
            query_kmer_cache=query_kmer_cache,
            library_kmer_cache=library_kmer_cache,
        )
        agreements.append(agreement)
    return agreements

def fast_generate_agreement(queries, library, alignment_cache, n_workers=1):
    """Compute agreement for a range of window sizes. Optionally parallelize heavy work by
    specifying n_workers>1 (uses process pool).
    """
    win = [i for i in range(20, 66, 5)]
    agreements = []
    if not queries or not library:
        return agreements

    fixed_k = 7



    query_kmer_cache = {key: convert_seq_to_kmer_set(seq, fixed_k) for key, seq in queries.items()}
    library_kmer_cache = ConvertLibraryToKmerSets(library, fixed_k)

    for w in win:
        if n_workers and n_workers > 1:
            min_library, minimizer_index = generate_min_library_parallel(library, m=w, k=fixed_k, n_workers=n_workers)
        else:
            min_library, minimizer_index = generate_min_library(library, m=w, k=fixed_k)

        agreement = find_minimizer_agreement(
            queries,
            min_library,
            library,
            m=w,
            k=fixed_k,
            minimizer_index=minimizer_index,
            alignment_cache=alignment_cache,
            query_kmer_cache=query_kmer_cache,
            library_kmer_cache=library_kmer_cache,
        )
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

# ============================================================
# Task 8: Benchmarking Runtime vs Agreement
# ============================================================

# Best-performing parameters found earlier
BEST_K = 7
BEST_M = 20
BEST_K_AGREEMENT = 0.76
BEST_MIN_AGREEMENT = 0.76

def _sample_items(d, n):
    items = list(d.items())
    random.shuffle(items)
    return items[:min(n, len(items))]


def benchmark_local_alignment(queries, library):
    """
    Measure the time for a single local alignment.
    Task 8 only times the alignment step itself, not a full library search.
    """
    if not queries or not library:
        return 0.0

    q_seq = next(iter(queries.values()))
    l_seq = next(iter(library.values()))

    start = time.time()
    alignment.local_align(q_seq, l_seq)
    return time.time() - start


def benchmark_kmer_best(queries, library, k=BEST_K, n_samples=50):
    """
    Measure k-mer (AFSM) runtime using the best k.
    K-mer sets are precomputed and not included in the timing.
    """
    library_kmer = ConvertLibraryToKmerSets(library, k)
    sampled = _sample_items(queries, n_samples)
    if not sampled:
        return 0.0

    query_kmers = [convert_seq_to_kmer_set(seq, k) for _, seq in sampled]

    start = time.time()
    for qk in query_kmers:
        KmerMatch(qk, library_kmer)
    return (time.time() - start) / len(sampled)


def build_minimizer_structures(library, m=BEST_M, k=BEST_K):
    """
    Precompute minimizers and library k-mers (not timed).
    """
    min_library, min_index = generate_min_library(library, m, k)
    library_kmer_cache = ConvertLibraryToKmerSets(library, k)
    return min_library, min_index, library_kmer_cache


def benchmark_minimizer_best(
    queries, min_library, minimizer_index, library_kmer_cache,
    m=BEST_M, k=BEST_K, n_samples=50
):
    """
    Measure minimizer + AFSM runtime using best (m, k).
    For each query:
      - compute minimizers
      - filter library candidates
      - run AFSM only on candidates
    """
    sampled = _sample_items(queries, n_samples)
    if not sampled:
        return 0.0

    query_data = []
    for _, seq in sampled:
        qk = convert_seq_to_kmer_set(seq, k)
        qm = get_minimizer(seq, m, k)
        if qm:
            query_data.append((qk, qm))

    if not query_data:
        return 0.0

    start = time.time()
    for qk, qm in query_data:
        candidates = set()
        for mini in qm:
            candidates.update(minimizer_index.get(mini, ()))

        subset = {key: library_kmer_cache[key] for key in candidates}
        if subset:
            KmerMatch(qk, subset)

    return (time.time() - start) / len(query_data)


def run_task8_benchmarks(queries, library):
    """
    Run runtime + agreement comparisons for:
      - Local Alignment
      - AFSM (best k)
      - Minimizers + AFSM (best m,k)
    """
    results = []

    # Local Alignment
    local_rt = benchmark_local_alignment(queries, library)
    results.append({
        "Method": "Local Alignment",
        "Params": "-",
        "Runtime (sec/query)": local_rt,
        "Agreement": 1.0,
    })

    # AFSM
    kmer_rt = benchmark_kmer_best(queries, library, k=BEST_K)
    results.append({
        "Method": "AFSM (k-mer)",
        "Params": f"k = {BEST_K}",
        "Runtime (sec/query)": kmer_rt,
        "Agreement": BEST_K_AGREEMENT,
    })

    # Minimizers + AFSM
    min_library, min_index, lib_kmer_cache = build_minimizer_structures(
        library, BEST_M, BEST_K
    )
    min_rt = benchmark_minimizer_best(
        queries, min_library, min_index, lib_kmer_cache,
        m=BEST_M, k=BEST_K
    )
    results.append({
        "Method": "Minimizers + AFSM",
        "Params": f"m = {BEST_M}, k = {BEST_K}",
        "Runtime (sec/query)": min_rt,
        "Agreement": BEST_MIN_AGREEMENT,
    })

    return results

def plot_task8_bars(results):
    """
    Bar plots for runtime and agreement for the three methods.
    """
    methods = [r["Method"] for r in results]
    runtimes = [r["Runtime (sec/query)"] for r in results]
    agreements = [r["Agreement"] for r in results]

    # Runtime plot
    plt.figure(figsize=(10, 5))
    bars = plt.bar(methods, runtimes, color=["#4e79a7", "#59a14f", "#9c755f"])
    plt.title("Task 8: Runtime per Query")
    plt.ylabel("Seconds per Query")
    for bar, rt in zip(bars, runtimes):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                 f"{rt:.4f}s", ha="center", va="bottom")
    plt.tight_layout()
    plt.show()

    # Agreement plot
    plt.figure(figsize=(10, 5))
    bars = plt.bar(methods, agreements, color=["#f28e2b", "#76b7b2", "#e15759"])
    plt.title("Task 8: Agreement")
    plt.ylabel("Agreement Fraction")
    plt.ylim(0, 1.1)
    for bar, ag in zip(bars, agreements):
        plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                 f"{ag:.3f}", ha="center", va="bottom")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    fn = "bacterial_16s_genes.fa"
    sequences_16s = Load16SFastA(fn, fraction = 0.015)

    # print("Loaded %d 16s sequences." % len(sequences_16s))

    # library, queries = split_dataset(sequences_16s, 10, 5)
    library, queries, alignment_match = load_datasets()
    print("Library of length %d, Queries of length %d" % (len(library), len(queries)))

    # agreements = find_agreement(queries, library)
    agreements = find_agreement_fast(queries, library, alignment_match)
    print(agreements)
    print("WE DID IT")

    ks = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
    plt.plot(ks, agreements)
    plt.title("Agreement Curve")
    plt.xlabel("K-mer length")
    plt.ylabel("Agreement")
    plt.xticks([1, 3, 5, 7, 9, 11, 13, 15, 17, 19])
    plt.show()

    #Task 6: Agreement Curve with minimizer

    # library, queries = split_dataset(sequences_16s, 200, 50)
    # min_agreements = generate_agreement(queries, library)
    # plt.plot([i for i in range(20, 66, 5)], min_agreements)
    # plt.title("Minimizer Agreement Curve")
    # plt.xlabel("Window size (m)")
    # plt.ylabel("Agreement")
    # plt.xticks([i for i in range(20, 66, 5)])
    # plt.show()

    # #Task 7
    ks = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
    print("=============== Testing with 1% mutated queries ===============")
    library, queries, alignment_match = load_datasets()
    print("Library of length %d, Queries of length %d" % (len(library), len(queries)))

    for key in queries.keys():
        queries[key] = mutate_string(queries[key], 1.0/100)

    agreements = find_agreement_fast(queries, library, alignment_match)
    print(agreements)

    plt.plot(ks, agreements)
    plt.title("Agreement Curve- 1% mutated")
    plt.xlabel("K-mer length")
    plt.ylabel("Agreement")
    plt.xticks(ks)
    plt.show()


    print("=============== Testing with 10% mutated queries ===============")
    library, queries, alignment_match = load_datasets()
    print("Library of length %d, Queries of length %d" % (len(library), len(queries)))

    for key in queries.keys():
        queries[key] = mutate_string(queries[key], 1.0/10)

    agreements = find_agreement_fast(queries, library, alignment_match)
    print(agreements)

    plt.plot(ks, agreements)
    plt.title("Agreement Curve- 10% mutated")
    plt.xlabel("K-mer length")
    plt.ylabel("Agreement")
    plt.xticks(ks)
    plt.show()

    print("============================== TASK 7 - 12 units =====================================")
    print("=============== Testing with 1% mutated queries ===============")
    library, queries, alignment_match = load_datasets()
    print("Library of length %d, Queries of length %d" % (len(library), len(queries)))

    for key in queries.keys():
        queries[key] = mutate_string(queries[key], 1.0 / 100)

    min_agreements = fast_generate_agreement(queries, library, alignment_match)
    plt.plot([i for i in range(20, 66, 5)], min_agreements)
    plt.title("Minimizer Agreement Curve for 1% Mutated")
    plt.xlabel("Window size (m)")
    plt.ylabel("Agreement")
    plt.xticks([i for i in range(20, 66, 5)])
    plt.show()

    print("=============== Testing with 10% mutated queries ===============")
    library, queries, alignment_match = load_datasets()
    print("Library of length %d, Queries of length %d" % (len(library), len(queries)))

    for key in queries.keys():
        queries[key] = mutate_string(queries[key], 1.0 / 10)

    min_agreements = fast_generate_agreement(queries, library, alignment_match)
    plt.plot([i for i in range(20, 66, 5)], min_agreements)
    plt.title("Minimizer Agreement Curve for 10% Mutated")
    plt.xlabel("Window size (m)")
    plt.ylabel("Agreement")
    plt.xticks([i for i in range(20, 66, 5)])
    plt.show()



    print("============= TASK 8 BENCHMARKING ONLY =============")
    library, queries, alignment_match = load_datasets()
    results = run_task8_benchmarks(queries, library)
    plot_task8_bars(results)
