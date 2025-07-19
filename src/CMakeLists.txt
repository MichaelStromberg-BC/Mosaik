# Include common source definitions
add_subdirectory(CommonSource)

# Define common source files that many executables need
set(COMMON_UTILITY_SOURCES
    "CommonSource/Utilities/Benchmark.cpp"
    "CommonSource/Utilities/ConsoleUtilities.cpp"
    "CommonSource/Utilities/FastLZIO.cpp"
    "CommonSource/Utilities/fastlz.c"
    "CommonSource/Utilities/FileUtilities.cpp"
    "CommonSource/Utilities/MemoryUtilities.cpp"
    "CommonSource/DataStructures/MosaikString.cpp"
    "CommonSource/Utilities/Options.cpp"
    "CommonSource/Utilities/SequenceUtilities.cpp"
    "CommonSource/Utilities/SHA1.cpp"
    "CommonSource/Utilities/TimeSupport.cpp"
    "CommonSource/Utilities/sqlite3.c"
)

# Define read format sources
set(READ_FORMAT_SOURCES
    "CommonSource/MosaikReadFormat/AlignmentReader.cpp"
    "CommonSource/MosaikReadFormat/AlignmentWriter.cpp"
    "CommonSource/MosaikReadFormat/ReadReader.cpp"
    "CommonSource/MosaikReadFormat/ReadWriter.cpp"
    "CommonSource/MosaikReadFormat/ReferenceSequenceReader.cpp"
)

# MosaikBuild
add_executable(MosaikBuild 
    MosaikBuild/BuildMain.cpp 
    MosaikBuild/MosaikBuild.cpp
    ${COMMON_UTILITY_SOURCES}
    "CommonSource/Utilities/ColorspaceUtilities.cpp"
    "CommonSource/ExternalReadFormats/Fasta.cpp"
    "CommonSource/ExternalReadFormats/Fastq.cpp"
    "CommonSource/Utilities/MD5.c"
    "CommonSource/MosaikReadFormat/ReadWriter.cpp"
    "CommonSource/Utilities/RegexUtilities.cpp"
    "CommonSource/ExternalReadFormats/SRF.cpp"
)
target_link_libraries(MosaikBuild ZLIB::ZLIB Threads::Threads)

# MosaikAligner
add_executable(MosaikAligner
    MosaikAligner/AlignerMain.cpp
    MosaikAligner/AlignmentThread.cpp  
    MosaikAligner/MosaikAligner.cpp
    ${COMMON_UTILITY_SOURCES}
    ${READ_FORMAT_SOURCES}
    "CommonSource/Utilities/AlignmentQuality.cpp"
    "CommonSource/DataStructures/AbstractDnaHash.cpp"
    "CommonSource/PairwiseAlignment/BandedSmithWaterman.cpp"
    "CommonSource/Utilities/ColorspaceUtilities.cpp"
    "CommonSource/DataStructures/DnaHash.cpp"
    "CommonSource/DataStructures/HashRegionTree.cpp"
    "CommonSource/DataStructures/JumpDnaHash.cpp"
    "CommonSource/DataStructures/MultiDnaHash.cpp"
    "CommonSource/DataStructures/NaiveAlignmentSet.cpp"
    "CommonSource/Utilities/PairwiseUtilities.cpp"
    "CommonSource/Utilities/RegexUtilities.cpp"
    "CommonSource/PairwiseAlignment/SmithWatermanGotoh.cpp"
    "CommonSource/DataStructures/UbiqDnaHash.cpp"
)
target_link_libraries(MosaikAligner Threads::Threads)

# MosaikSort
add_executable(MosaikSort
    MosaikSort/SortMain.cpp
    MosaikSort/PairedEndSort.cpp
    MosaikSort/SingleEndSort.cpp
    ${COMMON_UTILITY_SOURCES}
    ${READ_FORMAT_SOURCES}
    "CommonSource/DataStructures/NaiveAlignmentSet.cpp"
    "CommonSource/Utilities/AlignmentQuality.cpp"
)
target_link_libraries(MosaikSort Threads::Threads)

# MosaikMerge  
add_executable(MosaikMerge
    MosaikMerge/MergeMain.cpp
    MosaikMerge/MosaikMerge.cpp
    ${COMMON_UTILITY_SOURCES}
    ${READ_FORMAT_SOURCES}
    "CommonSource/DataStructures/NaiveAlignmentSet.cpp"
    "CommonSource/Utilities/AlignmentQuality.cpp"
)
target_link_libraries(MosaikMerge Threads::Threads)

# MosaikAssembler
add_executable(MosaikAssembler
    MosaikAssembler/AssemblerMain.cpp
    MosaikAssembler/MosaikAssembler.cpp
    ${COMMON_UTILITY_SOURCES}
    ${READ_FORMAT_SOURCES}
    "CommonSource/AssemblyFormats/AbstractAssemblyFormat.cpp"
    "CommonSource/AssemblyFormats/Ace.cpp"
    "CommonSource/AssemblyFormats/GigaBayesFormat.cpp"
    "CommonSource/DataStructures/NaiveAlignmentSet.cpp"
    "CommonSource/Utilities/AlignmentQuality.cpp"
)
target_link_libraries(MosaikAssembler Threads::Threads ZLIB::ZLIB)

# MosaikJump
add_executable(MosaikJump
    MosaikJump/JumpMain.cpp
    MosaikJump/JumpCreator.cpp
    ${COMMON_UTILITY_SOURCES}
    "CommonSource/DataStructures/JumpDnaHash.cpp"
    "CommonSource/DataStructures/AbstractDnaHash.cpp"
    "CommonSource/DataStructures/HashRegionTree.cpp"
    "CommonSource/MosaikReadFormat/ReferenceSequenceReader.cpp"
)
target_link_libraries(MosaikJump Threads::Threads)

# MosaikCoverage
add_executable(MosaikCoverage
    MosaikCoverage/CoverageMain.cpp
    MosaikCoverage/MosaikCoverage.cpp
    ${COMMON_UTILITY_SOURCES}
    ${READ_FORMAT_SOURCES}
    "CommonSource/DataStructures/NaiveAlignmentSet.cpp"
    "CommonSource/Utilities/AlignmentQuality.cpp"
)
target_link_libraries(MosaikCoverage ZLIB::ZLIB)

# MosaikText
add_executable(MosaikText
    MosaikText/TextMain.cpp
    MosaikText/MosaikText.cpp
    ${COMMON_UTILITY_SOURCES}
    ${READ_FORMAT_SOURCES}
    "CommonSource/ExternalReadFormats/BamWriter.cpp"
    "CommonSource/Utilities/ColorspaceUtilities.cpp"
    "CommonSource/DataStructures/NaiveAlignmentSet.cpp"
    "CommonSource/Utilities/AlignmentQuality.cpp"
)
target_link_libraries(MosaikText ZLIB::ZLIB Threads::Threads)

# MosaikDupSnoop
add_executable(MosaikDupSnoop
    MosaikDupSnoop/DupSnoopMain.cpp
    MosaikDupSnoop/MosaikDupSnoop.cpp
    ${COMMON_UTILITY_SOURCES}
    ${READ_FORMAT_SOURCES}
    "CommonSource/DataStructures/NaiveAlignmentSet.cpp"
    "CommonSource/Utilities/AlignmentQuality.cpp"
)
target_link_libraries(MosaikDupSnoop Threads::Threads)

# Set specific flags for sqlite3.c in all targets
set_source_files_properties("CommonSource/Utilities/sqlite3.c" PROPERTIES 
    COMPILE_FLAGS "-DSQLITE_OMIT_LOAD_EXTENSION"
)