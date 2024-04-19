
import numpy as np
import ahocorasick
import cython
from Bio import SeqIO
import os

# Generate the 18-mer array
@cython.boundscheck(False)
@cython.wraparound(False)
cdef get_kmers_arr(str individual_seq, int k):
    cdef int i
    cdef Py_ssize_t num_kmers = len(individual_seq) - k + 1
    kmers = np.zeros(num_kmers, dtype=np.dtype('U' + str(k)))
    for i in range(num_kmers):
        kmers[i] = individual_seq[i:i+k]
    return kmers

# Create index to store 18-mer arrays mapped to each UID
@cython.boundscheck(False)
@cython.wraparound(False)
cdef create_uid_index(str db_file, int k):
    cdef dict uid_index = {}
    cdef str s
    with open(db_file, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            s = str(record.seq)
            if len(s) > 0:
                uid_index[record.id] = get_kmers_arr(s, k)
            else:
                continue
    return uid_index


# Get length of the original miRNAs
cpdef uid_length_index(str db_file):
    cdef dict matches = {}
    cdef str s
    cdef int x
    with open(db_file, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            s = str(record.seq)
            x = len(s)
            temp_vector = np.zeros(x)
            matches[record.id] = temp_vector
    return matches

# Create dictionary of UIDs and original strings for those RNAs
cpdef original_strings(db_file):
    cdef dict records = {}
    cdef str s
    with open(db_file, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            s = str(record.seq)
            if len(s) > 0:
                records[record.id] = s
            else:
                continue
    return records

cpdef get_file_extension(file_path):
    return os.path.splitext(file_path)[1]

@cython.boundscheck(False)
@cython.wraparound(False)
# PREPROCESSING: Separate the sequences that definitely contain a match from those who definitely don't
cpdef separate_seqs(str db_file, str fasta_file, str contains_match_file, int k):
    ac = ahocorasick.Automaton()
    cdef int read_counts = 0
    cdef str line 
    with open(db_file, 'r') as f:
        for record in SeqIO.parse(f, "fasta"):
            line = str(record.seq)
            if len(line) > 0:
                arr = get_kmers_arr(line, k)
                for idx, mer in enumerate(arr):
                    ac.add_word(mer, (idx, mer))
            else:
                continue                 
    ac.make_automaton()
    in_file = open(fasta_file, "r")
    out_file = open(contains_match_file, "w")
    for line in in_file:
        read_counts += 1
        if line.startswith('A') or line.startswith('G') or line.startswith('T') or line.startswith('C'):
            for idx, mer in ac.iter(line):
                out_file.write(line)
                break
    ac = None
    in_file.close()
    out_file.close()
    ext = get_file_extension(fasta_file)
    if ext == '.fastq':
        return read_counts/4
    elif ext == '.fasta':
        return read_counts/2
    else:
        return 1


@cython.boundscheck(False)
@cython.wraparound(False)
# Aligns all sequences and returns hashmap with expression levels mapped to the UIDs
cpdef seq_aligner(str db_file, str known_match_file, int k):

    cdef int i
    cdef str key, test_string
    
    cdef dict uid_index = create_uid_index(db_file, k) 
    cdef dict matches = uid_length_index(db_file) 

    ac = ahocorasick.Automaton()
    ac1 = ahocorasick.Automaton()
    ac2 = ahocorasick.Automaton()

    for key, ss_arr in uid_index.items():
        for i, mer in enumerate(ss_arr):
            if not mer in ac:
                ac.add_word(mer, (i, key))
            elif not mer in ac1: 
                ac1.add_word(mer, (i, key))
            else:
                ac2.add_word(mer, (i, key))
    
    ac.make_automaton()
    ac1.make_automaton()
    ac2.make_automaton()

    with open(known_match_file, 'r') as in_file: 
        for test_string in in_file:
            if len(ac) != 0:
                for indx, (i, key) in ac.iter(test_string):
                    matches[key][i] += 1
            if len(ac1) != 0:
                for indx, (i, key) in ac1.iter(test_string):
                    matches[key][i] += 1
            if len(ac2) != 0:
                for indx, (i, key) in ac2.iter(test_string):
                    matches[key][i] += 1
    ac = None
    ac1 = None
    ac2 = None
    return matches
