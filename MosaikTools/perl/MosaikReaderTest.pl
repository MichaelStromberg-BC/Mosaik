#!/usr/bin/perl -w
# ***************************************************************************
# MosaikReaderTest.pl - loads and displays all of the alignments in a MOSAIK
#                       alignment archive.
# ---------------------------------------------------------------------------
# (c) 2006 - 2009 Michael Strömberg
# Marth Lab, Department of Biology, Boston College
# ---------------------------------------------------------------------------
# Dual licenced under the GNU General Public License 2.0+ license or as
# a commercial license with the Marth Lab.
# ***************************************************************************

use strict;
use mosaik;

# displays the info found in the Alignment object
sub DisplayAlignment {

    # retrieve the Alignment object
    my $alignment = shift;
    my $mate      = shift;

    print "mate $mate: Reference: ", mosaikc::Alignment_ReferenceName_get($alignment), 
    ", orientation: ", (mosaikc::Alignment_IsReverseStrand_get($alignment) ? 'R' : 'F'), 
    ", alignment quality: " , mosaikc::Alignment_Quality_get($alignment) ,"\n";

    printf "%10u %s %u (reference)\n", mosaikc::Alignment_ReferenceBegin_get($alignment), mosaikc::Alignment_Reference_get($alignment), mosaikc::Alignment_ReferenceEnd_get($alignment);
    printf "%10u %s %u (query)\n\n", mosaikc::Alignment_QueryBegin_get($alignment), mosaikc::Alignment_Query_get($alignment), mosaikc::Alignment_QueryEnd_get($alignment);
}


# displays the info found in the AlignedRead object
sub DisplayRead {

    # retrieve the AlignedRead object
    my $ar = shift;
    my $readGroupID = shift;

    print "--------------------------------------------------\n";
    print "read name: ", mosaikc::AlignedRead_Name_get($ar), "\n\n";

    my $mate1Alignments = mosaikc::AlignedRead_Mate1Alignments_get($ar);
    my $mate2Alignments = mosaikc::AlignedRead_Mate2Alignments_get($ar);

    my $numMate1Alignments = $mate1Alignments->size();
    my $numMate2Alignments = $mate2Alignments->size();

    print "# of mate 1 alignments: $numMate1Alignments\n";
    print "# of mate 2 alignments: $numMate2Alignments\n";

    # =============================
    # display the mate 1 alignments
    # =============================

    for(my $i = 0; $i < $numMate1Alignments; $i++) {
	DisplayAlignment($mate1Alignments->get($i), 1);
    }

    # =============================
    # display the mate 2 alignments
    # =============================

    for(my $i = 0; $i < $numMate2Alignments; $i++) {
	DisplayAlignment($mate2Alignments->get($i), 2);
    }
}

# grab the command line arguments
my $numArgs = $#ARGV + 1;

if($numArgs != 1) {
    print "USAGE: $0 <MOSAIK alignment archive>\n";
    exit;
}

# localize our arguments
my $filename = shift;

# check if this is a valid MOSAIK alignment file (optional)
mosaik::CAlignmentReader::CheckFile($filename, 1);

# open the MOSAIK alignments file
my $reader = new mosaik::CAlignmentReader();
$reader->Open($filename);

# get some basic statistics
my $numBases = $reader->GetNumBases();
my $numReads = $reader->GetNumReads();

print "# of bases: $numBases\n";
print "# of reads: $numReads\n\n";

# retrieve the reference sequences, read groups, and alignment status
my $referenceSequences = $reader->GetReferenceData();
my $readGroups         = $reader->GetReadGroups();

# test the status bits
my $alignStatus = $reader->GetStatus();

if($alignStatus == $mosaik::AS_UNKNOWN) {
    print "The archive has an unknown alignment status.\n";
} else {
    if(($alignStatus & $mosaik::AS_SINGLE_END_READ) != 0) {
        print "Archive contains single-end reads.\n";
    }
    if(($alignStatus & $mosaik::AS_PAIRED_END_READ) != 0) {
        print "Archive contains paired-end reads.\n";
    }
    if(($alignStatus & $mosaik::AS_UNSORTED_READ) != 0) {
        print "Archive contains unsorted reads.\n";
    }
    if(($alignStatus & $mosaik::AS_SORTED_ALIGNMENT) != 0) {
        print "Archive contains sorted alignments.\n";
    }
    if(($alignStatus & $mosaik::AS_ALL_MODE) != 0) {
        print "Archive was generated using 'all' mode.\n";
    }
    if(($alignStatus & $mosaik::AS_UNIQUE_MODE) != 0) {
        print "Archive was generated using 'unique' mode.\n";
    }
}

# dump the reference sequences used
print "\nDumping reference sequences:\n";

my $numRefSeqs = $referenceSequences->size();
for(my $i = 0; $i < $numRefSeqs; $i++) {
    
    my $refSeq = $referenceSequences->get($i);
    printf "- %s, aligned reads: %u, length: %u", mosaikc::ReferenceSequence_Name_get($refSeq), mosaikc::ReferenceSequence_NumAligned_get($refSeq), mosaikc::ReferenceSequence_NumBases_get($refSeq);

    # grab some optional parameters
    my $genomeAssemblyID = mosaikc::ReferenceSequence_GenomeAssemblyID_get($refSeq);
    my $md5              = mosaikc::ReferenceSequence_MD5_get($refSeq);
    my $species          = mosaikc::ReferenceSequence_Species_get($refSeq);
    my $uri              = mosaikc::ReferenceSequence_URI_get($refSeq);

    if($genomeAssemblyID) {
	print ", genome assembly ID: $genomeAssemblyID";
    }

    if($md5) {
	print ", MD5 checksum: $md5";
    }

    if($species) {
	print ", species: $species";
    }

    if($uri) {
	print ", URI: $uri";
    }

    print "\n";
}

# dump the read groups used
print "\nDumping read groups:\n";

my $numReadGroups = $readGroups->size();
for(my $i = 0; $i < $numReadGroups; $i++) {

    my $rg = $readGroups->get($i);
    printf "%s: name: %s, sample: %s", mosaikc::ReadGroup_ReadGroupCode_get($rg), mosaikc::ReadGroup_ReadGroupID_get($rg), mosaikc::ReadGroup_SampleName_get($rg);

    # grab some optional parameters
    my $centerName           = mosaikc::ReadGroup_CenterName_get($rg);
    my $description          = mosaikc::ReadGroup_Description_get($rg);
    my $libraryName          = mosaikc::ReadGroup_LibraryName_get($rg);
    my $medianFragmentLength = mosaikc::ReadGroup_MedianFragmentLength_get($rg);
    my $platformUnit         = mosaikc::ReadGroup_PlatformUnit_get($rg);
    my $seqTechnology        = mosaikc::ReadGroup_SequencingTechnology_get($rg);

    if($centerName) {
        print ", center: $centerName";
    }

    if($description) {
        print ", description: $description";
    }

    if($libraryName) {
        print ", library: $libraryName";
    }

    if($medianFragmentLength) {
        print ", median fragment length: $medianFragmentLength";
    }

    if($platformUnit) {
        print ", platform unit: $platformUnit";
    }

    if($seqTechnology == $mosaikc::ST_454) {
	print ", technology: 454\n";
    } elsif($seqTechnology == $mosaikc::ST_HELICOS) {
	print ", technology: Helicos\n";
    } elsif($seqTechnology == $mosaikc::ST_ILLUMINA) {
	print ", technology: Illumina\n";
    } elsif($seqTechnology == $mosaikc::ST_PACIFIC_BIOSCIENCES) {
	print ", technology: Pacific Biosciences\n";
    } elsif($seqTechnology == $mosaikc::ST_SOLID) {
	print ", technology: AB SOLiD\n";
    } elsif($seqTechnology == $mosaikc::ST_SANGER) {
	print ", technology: Sanger capillary\n";
    } else {
	print ", technology: unknown\n";
    }
}
print "\n";

# keep reading all of the alignments
my $ar = new mosaik::AlignedRead();
while($reader->LoadNextRead($ar)) {
    DisplayRead($ar);
}

# close the alignments file
$reader->Close();

# destroy our objects
$reader->DESTROY();
