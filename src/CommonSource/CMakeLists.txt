# Define all common source files that executables can link against

# Assembly Formats sources
set(ASSEMBLY_SOURCES
    AssemblyFormats/AbstractAssemblyFormat.cpp
    AssemblyFormats/Ace.cpp
    AssemblyFormats/GigaBayesFormat.cpp
)

# Data Structures sources
set(DATA_STRUCTURES_SOURCES
    DataStructures/AbstractDnaHash.cpp
    DataStructures/DnaHash.cpp
    DataStructures/HashRegionTree.cpp
    DataStructures/JumpDnaHash.cpp
    DataStructures/MosaikString.cpp
    DataStructures/MultiDnaHash.cpp
    DataStructures/NaiveAlignmentSet.cpp
    DataStructures/UbiqDnaHash.cpp
)

# External Read Formats sources
set(EXTERNAL_FORMATS_SOURCES
    ExternalReadFormats/Fasta.cpp
    ExternalReadFormats/Fastq.cpp
    ExternalReadFormats/SRF.cpp
    ExternalReadFormats/BamWriter.cpp
)

# MOSAIK Read Format sources
set(MOSAIK_READ_SOURCES
    MosaikReadFormat/AlignmentReader.cpp
    MosaikReadFormat/AlignmentWriter.cpp
    MosaikReadFormat/ReadReader.cpp
    MosaikReadFormat/ReadWriter.cpp
    MosaikReadFormat/ReferenceSequenceReader.cpp
)

# Pairwise Alignment sources
set(PAIRWISE_SOURCES
    PairwiseAlignment/SmithWatermanGotoh.cpp
    PairwiseAlignment/BandedSmithWaterman.cpp
)

# Utilities C++ sources
set(UTILITIES_SOURCES
    Utilities/AlignmentQuality.cpp
    Utilities/Benchmark.cpp
    Utilities/ColorspaceUtilities.cpp
    Utilities/ConsoleUtilities.cpp
    Utilities/FastLZIO.cpp
    Utilities/FileUtilities.cpp
    Utilities/MemoryUtilities.cpp
    Utilities/Options.cpp
    Utilities/PairwiseUtilities.cpp
    Utilities/RegexUtilities.cpp
    Utilities/SequenceUtilities.cpp
    Utilities/SHA1.cpp
    Utilities/TimeSupport.cpp
)

# Utilities C sources
set(UTILITIES_C_SOURCES
    Utilities/fastlz.c
    Utilities/MD5.c
    Utilities/sqlite3.c
)

# Make these available to parent scope
set(COMMON_SOURCES
    ${ASSEMBLY_SOURCES}
    ${DATA_STRUCTURES_SOURCES}
    ${EXTERNAL_FORMATS_SOURCES}
    ${MOSAIK_READ_SOURCES}
    ${PAIRWISE_SOURCES}
    ${UTILITIES_SOURCES}
    ${UTILITIES_C_SOURCES}
    PARENT_SCOPE
)

# Set specific flags for sqlite3.c
set_source_files_properties(Utilities/sqlite3.c PROPERTIES 
    COMPILE_FLAGS "-DSQLITE_OMIT_LOAD_EXTENSION"
)