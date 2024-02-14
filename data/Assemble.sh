#!/bin/sh

# create a multiple sequence alignment (reference-guided assembly)
../bin/MosaikAssembler -in sequence_archives/c_elegans_chr2_test_sorted.dat -out assembly/c.elegans_chr2_test -ia reference/c.elegans_chr2.dat -f ace
