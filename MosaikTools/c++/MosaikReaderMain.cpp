// ***************************************************************************
// MosaikReaderMain.cpp - loads and displays all of the alignments in a
//                        MOSAIK alignment archive.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <iostream>
#include <vector>
#include <cstdlib>
#include "MosaikAlignment.h"

using namespace std;

// function prototypes
void DisplayAlignedRead(Mosaik::AlignedRead& ar, const string& readGroupID);
void DisplayAlignment(string& readName, vector<Mosaik::Alignment>::iterator& alIter);

int main(int argc, char* argv[]) {

	if(argc != 2) {
		cout << "USAGE: " << argv[0] << " <input MOSAIK file>" << endl;
		exit(1);
	}

	// localize our arguments
	const char* inputFilename  = argv[1];

	// Check if this is a valid MOSAIK alignment file (optional)
	Mosaik::CAlignmentReader::CheckFile(inputFilename, true);

	// open our MOSAIK alignment reader
	Mosaik::CAlignmentReader reader;
	reader.Open(inputFilename);

	// get some basic statistics
	uint64_t numBases = reader.GetNumBases();
	uint64_t numReads = reader.GetNumReads();

	cout << "# of bases: " << numBases << endl;
	cout << "# of reads: " << numReads << endl << endl;

	// retrieve the reference sequences, read groups, and alignment status
	Mosaik::RefVector referenceSequences    = reader.GetReferenceData();
	Mosaik::ReadGroupVector readGroups      = reader.GetReadGroups();
	Mosaik::AlignmentStatus alignmentStatus = reader.GetStatus();

	// dump the alignment status
	if(alignmentStatus == AS_UNKNOWN) {
		cout << "The archive has an unknown alignment status." << endl;
	} else {
		if((alignmentStatus & AS_SINGLE_END_READ)  != 0) cout << "Archive contains single-end reads." << endl;
		if((alignmentStatus & AS_PAIRED_END_READ)  != 0) cout << "Archive contains paired-end reads." << endl;
		if((alignmentStatus & AS_UNSORTED_READ)    != 0) cout << "Archive contains unsorted reads." << endl;
		if((alignmentStatus & AS_SORTED_ALIGNMENT) != 0) cout << "Archive contains sorted alignments." << endl;
		if((alignmentStatus & AS_ALL_MODE)         != 0) cout << "Archive was generated using 'all' mode." << endl;
		if((alignmentStatus & AS_UNIQUE_MODE)      != 0) cout << "Archive was generated using 'unique' mode." << endl;
	}

	// dump the reference sequences used
	cout << endl << "Dumping reference sequences:" << endl;
	Mosaik::RefVector::const_iterator rsIter;
	for(rsIter = referenceSequences.begin(); rsIter != referenceSequences.end(); rsIter++) {
		cout << "- " << rsIter->Name << ", aligned reads: " << rsIter->NumAligned << ", length: " << rsIter->NumBases;
		
		if(!rsIter->GenomeAssemblyID.empty()) cout << ", genome assembly ID: " << rsIter->GenomeAssemblyID;
		if(!rsIter->MD5.empty())              cout << ", MD5 checksum: "       << rsIter->MD5;
		if(!rsIter->Species.empty())          cout << ", species: "            << rsIter->Species;
		if(!rsIter->URI.empty())              cout << ", URI: "                << rsIter->URI;
		cout << endl;
	}

	// dump the read groups used
	cout << endl << "Dumping read groups:" << endl;
	Mosaik::ReadGroupVector::const_iterator rgIter;
	for(rgIter = readGroups.begin(); rgIter != readGroups.end(); rgIter++) {
		cout << rgIter->ReadGroupCode << ": name: " << rgIter->ReadGroupID << ", sample: " << rgIter->SampleName;
		
		if(!rgIter->CenterName.empty())       cout << ", center: "                 << rgIter->CenterName;
		if(!rgIter->Description.empty())      cout << ", description: "            << rgIter->Description;
		if(!rgIter->LibraryName.empty())      cout << ", library: "                << rgIter->LibraryName;
		if(rgIter->MedianFragmentLength != 0) cout << ", median fragment length: " << rgIter->MedianFragmentLength;
		if(!rgIter->PlatformUnit.empty())     cout << ", platform unit: "          << rgIter->PlatformUnit;
		
		switch(rgIter->SequencingTechnology) {
			case ST_454:
				cout << ", technology: 454" << endl;
				break;
			case ST_HELICOS:
				cout << ", technology: Helicos" << endl;
				break;
			case ST_ILLUMINA:
				cout << ", technology: Illumina" << endl;
				break;
			case ST_PACIFIC_BIOSCIENCES:
				cout << ", technology: Pacific Biosciences" << endl;
				break;
			case ST_SOLID:
				cout << ", technology: AB SOLiD" << endl;
				break;
			case ST_SANGER:
				cout << ", technology: Sanger capillary" << endl;
				break;
			default:
				cout << ", technology: unknown" << endl;
		}
	}
	cout << endl;

	// dump all of the reads from the input file to the screen
	Mosaik::AlignedRead ar;
	while(reader.LoadNextRead(ar)) {
		const string readGroupID = reader.GetReadGroupFromCode(ar.ReadGroupCode).ReadGroupID;	
		DisplayAlignedRead(ar, readGroupID);
	}

	// close our files
	reader.Close();

	return 0;
}

// displays the data present in an aligned read data structure
void DisplayAlignedRead(Mosaik::AlignedRead& ar, const string& readGroupID) {

	// initialize
	vector<Mosaik::Alignment>::iterator alIter;

	// display aligned read information
	cout << "--------------------------------------------------" << endl;
	cout << "Name: " << ar.Name << ", read group: " << readGroupID << ", # mate 1 alignments: "
		<< ar.Mate1Alignments.size() << ", # mate 2 alignments: " << ar.Mate2Alignments.size() << endl << endl;

	// dump the mate 1 alignments
	for(alIter = ar.Mate1Alignments.begin(); alIter != ar.Mate1Alignments.end(); alIter++) DisplayAlignment(ar.Name, alIter);

	// dump the mate 2 alignments
	for(alIter = ar.Mate2Alignments.begin(); alIter != ar.Mate2Alignments.end(); alIter++) DisplayAlignment(ar.Name, alIter);
}

// displays the data present in an alignment data structure
void DisplayAlignment(string& readName, vector<Mosaik::Alignment>::iterator& alIter) {

	// increment the base qualities
	char* pBQ = (char*)alIter->BaseQualities.data();
	for(unsigned int i = 0; i < alIter->BaseQualities.size(); i++) pBQ[i] += 33;

	// define the mate type
	string mateType = "SE";
	if(alIter->IsPairedEnd) {
		if(alIter->IsResolvedAsPair) {
			if(alIter->IsFirstMate) mateType = "resolved mate 1";
			else mateType = "resolved mate 2";
		} else {
			if(alIter->IsFirstMate) mateType = "mate 1";
			else mateType = "mate 2";
		}
	}

	printf("%s (%u) %u %u; %s (%s) %u %u; %c %u; %u mm\n", alIter->ReferenceName, alIter->ReferenceIndex,
		alIter->ReferenceBegin + 1, alIter->ReferenceEnd + 1, readName.c_str(), mateType.c_str(), alIter->QueryBegin + 1,
		alIter->QueryEnd + 1, (alIter->IsReverseStrand ? '-' : '+'), alIter->Quality, alIter->NumMismatches);
	printf("%s\n", alIter->Reference.c_str());
	printf("%s\n", alIter->Query.c_str());
	printf("%s\n", alIter->BaseQualities.c_str());

	// show if the alignment was rescued
	if(alIter->WasRescued) printf("RESCUED.\n");

	// show mate information
	if(alIter->IsResolvedAsPair) {
		printf("mate reference: %u, begin: %u, end: %u, orientation: %c\n", alIter->MateReferenceIndex, alIter->MateReferenceBegin + 1, alIter->MateReferenceEnd + 1, (alIter->IsMateReverseStrand ? '-' : '+'));
	}

	printf("\n");
}
