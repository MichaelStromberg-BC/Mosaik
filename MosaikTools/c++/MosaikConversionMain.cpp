// ***************************************************************************
// MosaikConversionMain.cpp - loads all of the reads from one MOSAIK alignment
//                            archive and saves them to another.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <iostream>
#include <cstdlib>
#include "MosaikAlignment.h"

using namespace std;

int main(int argc, char* argv[]) {

	if(argc != 3) {
		cout << "USAGE: " << argv[0] << " <input MOSAIK file> <output MOSAIK file>" << endl;
		exit(1);
	}

	// localize our arguments
	const char* inputFilename  = argv[1];
	const char* outputFilename = argv[2];

	// open our MOSAIK alignment reader
	Mosaik::CAlignmentReader reader;
	reader.Open(inputFilename);

	// retrieve the reference sequences, read groups, and alignment status
	Mosaik::RefVector referenceSequences    = reader.GetReferenceData();
	Mosaik::ReadGroupVector readGroups      = reader.GetReadGroups();
	Mosaik::AlignmentStatus alignmentStatus = reader.GetStatus();

	// open the MOSAIK alignment writer
	Mosaik::CAlignmentWriter writer;
	writer.Open(outputFilename, referenceSequences, readGroups, alignmentStatus);

	// copy all of the reads from the input file to the output file
	Mosaik::AlignedRead ar;
	while(reader.LoadNextRead(ar)) writer.SaveAlignedRead(ar);

	// close our files
	reader.Close();
	writer.Close();

	return 0;
}
