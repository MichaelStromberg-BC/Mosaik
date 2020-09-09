// ***************************************************************************
// CJumpCreator - creates a jump database for use with MosaikAligner.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "ConsoleUtilities.h"
#include "FileUtilities.h"
#include "ProgressBar.h"
#include "ReferenceSequenceReader.h"
#include "TimeSupport.h"

using namespace std;

#define KEY_LENGTH                       5
#define DEFAULT_SORTING_MEMORY           2
#define RESERVED_SPACE                  30
#define REFERENCE_SEQUENCE_INDEX_OFFSET 50

struct HashPosition {
	uint64_t Hash;
	unsigned int Position;
	unsigned char Owner;

	// deserialize this object from the supplied file stream
	bool Deserialize(FILE* temp) {
		fread((char*)&Hash,     SIZEOF_UINT64, 1, temp);
		fread((char*)&Position, SIZEOF_INT,    1, temp);

		if(feof(temp)) return false;
		return true;
	}

	// serialize this object to the supplied file stream
	void Serialize(FILE* temp) {
		fwrite((char*)&Hash,     SIZEOF_UINT64, 1, temp);
		fwrite((char*)&Position, SIZEOF_INT,    1, temp);
	}
};

struct ReferenceSequenceIndexEntry {
	off_type BeginOffset;
	off_type EndOffset;
	unsigned int NumHashes;

	// constructor
	ReferenceSequenceIndexEntry(void)
		: BeginOffset(0)
		, EndOffset(0)
		, NumHashes(0)
	{}
};

typedef vector<string>                      TempFilenames_t;
typedef vector<TempFilenames_t>             TempReferenceFilenames_t;
typedef vector<unsigned char>               Buffer_t;
typedef vector<HashPosition>                HashPositions_t;
typedef vector<FILE*>                       FileHandles_t;   
typedef vector<ReferenceSequenceIndexEntry> ReferenceSequenceIndex_t;

class CJumpCreator {
public:
	// constructor
	CJumpCreator(const unsigned char hashSize, const string& filename, const unsigned int hashPositionThreshold);
	// destructor
	~CJumpCreator(void);
	// builds the jump database
	void BuildJumpDatabase(void);
	// hashes the reference sequences in the specified reference sequence archive
	void HashReferenceSequences(const string& referenceFilename);
	// saves the metadata to the jump database
	void WriteMetadata(void);

private:
	// define a comparison function for sorting our hash positions (ascending)
	struct SortHashPositionAsc {
		bool operator()(const HashPosition& hp1, const HashPosition& hp2) {
			if(hp1.Hash == hp2.Hash) return hp1.Position < hp2.Position;
			return hp1.Hash < hp2.Hash;
		}
	};
	// define a comparison function for sorting our hash positions (descending)
	struct SortHashPositionDesc {
		bool operator()(const HashPosition& hp1, const HashPosition& hp2) {
			if(hp1.Hash == hp2.Hash) return hp2.Position < hp1.Position;
			return hp2.Hash < hp1.Hash;
		}
	};
	// checks that the IO buffer can hold the specified number of bytes
	inline void CheckBufferSize(unsigned int numBytes);
	// creates the hash for a supplied fragment
	static void CreateHash(const char* fragment, const unsigned char fragmentLen, uint64_t& key);
	// hashes the reference sequence and stores the results in sorted temporary files
	void HashReferenceSequence(const string& bases);
	// serializes the hash positions to temporary files
	void SerializeHashPositions(void);
	// stores the supplied hash positions in the jump database
	void StoreHash(HashPositions_t& hashPositions);
	// stores the all of the serialized filenames used
	TempReferenceFilenames_t mTempReferenceFilenames;
	// our hash size
	unsigned char mHashSize;
	// our jump database file handles
	FILE* mJumpOut;
	// our output buffer
	Buffer_t mBuffer;
	// toggles whether or not we return all genome positions or just a subset
	bool mLimitPositions;
	// our hash positions buffer
	HashPositions_t mHashPositions;
	int mMaxHashPositions;
	// the number of reference sequences in the reference sequence archive
	unsigned int mNumRefSeqs;
};

// checks that the IO buffer can hold the specified number of bytes
inline void CJumpCreator::CheckBufferSize(unsigned int numBytes) {
	if((unsigned int)mBuffer.size() < numBytes) mBuffer.resize(numBytes);
}