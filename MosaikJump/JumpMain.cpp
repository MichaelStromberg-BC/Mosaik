// ***************************************************************************
// JumpMain.cpp - collects command-line parameters for MosaikJump.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include <iostream>
#include <iomanip>
#include "Benchmark.h"
#include "ConsoleUtilities.h"
#include "ConversionUtilities.h"
#include "JumpCreator.h"
#include "Mosaik.h"
#include "Options.h"
#include "ReferenceSequenceReader.h"

using namespace std;

// define our default values
#define MIN_HASH_SIZE     4
#define MAX_HASH_SIZE     32

// create a configuration variable struct
struct ConfigurationSettings {

	// flags
	bool HasJumpFilename;
	bool HasHashPositionsFilename;
	bool HasHashSize;
	bool HasReferenceFilename;
	bool LimitHashPositions;

	// filenames
	string ReferenceFilename;
	string JumpFilename;

	// parameters
	unsigned int HashPositionThreshold;
	unsigned int HashSize;

	// constructor
	ConfigurationSettings()
		: HasJumpFilename(false)
		, HasHashSize(false)
		, HasReferenceFilename(false)
		, LimitHashPositions(false)
	{}
};

int main(int argc, char* argv[]) {

	CConsole::Initialize();
	ConfigurationSettings settings;

	printf("------------------------------------------------------------------------------\n");
	printf("Mosaik"); CConsole::Red(); printf("Jump"); CConsole::Reset();
	printf(" %u.%u.%04u                                                 %s\n", 
		MOSAIK_MAJOR_VERSION, MOSAIK_MINOR_VERSION, MOSAIK_BUILD_VERSION, MOSAIK_VERSION_DATE);
	printf("Michael Stromberg                 Marth Lab, Boston College Biology Department\n");
	printf("------------------------------------------------------------------------------\n\n");

	// =================================
	// configure the command line parser
	// =================================

	// set general info about the program
	COptions::SetProgramInfo("MosaikJump", "produces a jump database from a MOSAIK reference file", "-ia <filename> -out <filename> -hs <hash size>");

	// add the input options
	OptionGroup* pInputOpts = COptions::CreateOptionGroup("Input");
	COptions::AddValueOption("-ia", "MOSAIK reference filename", "the input reference file", "An input MOSAIK reference file",  settings.HasReferenceFilename, settings.ReferenceFilename, pInputOpts);

	// add the output options
	OptionGroup* pOutputOpts = COptions::CreateOptionGroup("Output");
	COptions::AddValueOption("-out", "jump filename stub", "the stub for the output filenames", "A filename stub",  settings.HasJumpFilename, settings.JumpFilename, pOutputOpts);

	// add the options
	OptionGroup* pOpts = COptions::CreateOptionGroup("Options");
	COptions::AddValueOption("-hs",  "hash size",      "the hash size [4 - 32]",                     "The hash size", settings.HasHashSize,        settings.HashSize,              pOpts);
	COptions::AddValueOption("-mhp", "hash positions", "sets the max number of hash positions",      "",              settings.LimitHashPositions, settings.HashPositionThreshold, pOpts);

	// parse the current command line
	COptions::Parse(argc, argv);

	// =============================
	// check for missing information
	// =============================

	bool foundError = false;
	ostringstream errorBuilder;
	const string ERROR_SPACER(7, ' ');

	// check the hash position threshold
	if(settings.LimitHashPositions && (settings.HashPositionThreshold < 1)) {
		errorBuilder << ERROR_SPACER << "The hash position threshold should be larger than 0." << endl;
		foundError = true;
	}

	// check the hash size
	if(settings.HasHashSize && ((settings.HashSize < MIN_HASH_SIZE) || (settings.HashSize > MAX_HASH_SIZE))) {
		errorBuilder << ERROR_SPACER << "Hash size should be between " << MIN_HASH_SIZE << " and " << MAX_HASH_SIZE << ". Please revise with the -hs parameter." << endl;
		foundError = true;
	}

	// print the errors if any were found
	if(foundError) {

		CConsole::Red();
		printf("ERROR: Some problems were encountered when parsing the command line options:\n");
		CConsole::Reset();

		printf("%s\n", errorBuilder.str().c_str());
		printf("For a complete list of command line options, type \"%s -h\"\n", argv[0]);
		exit(1);
	}

	// ===================================================
	// Parse configuration strings and set class variables
	// ===================================================

	// test to see if the specified input files exist
	if(settings.HasReferenceFilename)
		MosaikReadFormat::CReferenceSequenceReader::CheckFile(settings.ReferenceFilename, true);

	// start benchmarking
	CBenchmark bench;
	bench.Start();

	CJumpCreator jc(settings.HashSize, settings.JumpFilename, settings.HashPositionThreshold);

	// hash the reference and store the results in sorted temporary files
	jc.HashReferenceSequences(settings.ReferenceFilename);

	// build the jump database
	jc.BuildJumpDatabase();

	// ==================
	// Show total runtime
	// ==================

	// stop benchmarking
	bench.Stop();

	// show the benchmarking results
	bench.DisplayTime("\nMosaikJump");

	return 0;
}
