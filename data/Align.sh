#!/bin/sh

# align the reads
../bin/MosaikAligner -in sequence_archives/c_elegans_chr2_test.dat -out sequence_archives/c_elegans_chr2_test_aligned.dat -ia reference/c.elegans_chr2.dat -hs 14 -act 17 -mm 2 -m unique

# sort the reads
../bin/MosaikSort -in sequence_archives/c_elegans_chr2_test_aligned.dat -out sequence_archives/c_elegans_chr2_test_sorted.dat