// ***************************************************************************
// CJumpCreator - creates a jump database for use with MosaikAligner.
// ---------------------------------------------------------------------------
// (c) 2006 - 2010 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "JumpCreator.h"

// constructor
CJumpCreator::CJumpCreator(const unsigned char hashSize, const string& filename, const unsigned int hashPositionThreshold)
	: mHashSize(hashSize)
	, mJumpOut(NULL)
	, mLimitPositions(false)
	, mMaxHashPositions(0)
	, mNumRefSeqs(0)
{
	// initialize the file buffer
	mBuffer.resize(4096);

	// open the jump database for output
	mJumpOut = fopen(filename.c_str(), "wb");

	if(!mJumpOut) {
		printf("ERROR: Unable to open the positions file (%s) for writing.\n", filename.c_str());
		exit(1);
	}

	// apply the hash position threshold
	if(hashPositionThreshold > 0) {
		cout << "- setting hash position threshold to " << hashPositionThreshold << endl;
		mLimitPositions   = true;
		mMaxHashPositions = hashPositionThreshold;
	}

	// allocate the hash positions memory
	const double memoryAllocated = DEFAULT_SORTING_MEMORY * 1073741824.0;
	mMaxHashPositions = (int)(memoryAllocated / (double)sizeof(HashPosition));
	mHashPositions.reserve(mMaxHashPositions);
}

// destructor
CJumpCreator::~CJumpCreator(void) {

	// close the jump database files
	fclose(mJumpOut);

	// delete our temporary files
	TempReferenceFilenames_t::const_iterator trfIter;
	TempFilenames_t::const_iterator tfIter;
	for(trfIter = mTempReferenceFilenames.begin(); trfIter != mTempReferenceFilenames.end(); ++trfIter) {
		for(tfIter = trfIter->begin(); tfIter != trfIter->end(); ++tfIter) rm(tfIter->c_str());
	}
}

// builds the jump database
void CJumpCreator::BuildJumpDatabase(void) {

	// ================
	// write the header
	// ================

	// MOSAIK_SIGNATURE[7]	   0  -  6
	// HASH_SIZE[1]            7  -  7
	// ARCHIVE_DATE[8]		   8  - 15
	// NUM_REFERENCE_SEQS[4]   16 - 19
	// RESERVED[30]            20 - 49

	// write the MOSAIK signature
	const unsigned char SIGNATURE_LENGTH = 7;
	const char* MOSAIK_SIGNATURE = "MSKJMP\0";
	fwrite(MOSAIK_SIGNATURE, SIGNATURE_LENGTH, 1, mJumpOut);

	// write the hash size
	putc(mHashSize, mJumpOut);

	// write the archive date
	uint64_t currentTime = CTimeSupport::GetSystemTime();
	fwrite((char*)&currentTime, SIZEOF_UINT64, 1, mJumpOut);

	// write the number of reference sequences
	fwrite((char*)&mNumRefSeqs, SIZEOF_INT, 1, mJumpOut);

	// skip the reserved space and reference sequence index
	const unsigned int REFSEQ_INDEX_SIZE = (2 * SIZEOF_UINT64 + SIZEOF_INT) * mNumRefSeqs;
	fseek64(mJumpOut, RESERVED_SPACE + REFSEQ_INDEX_SIZE, SEEK_CUR);

	// ================================
	// evaluate each reference sequence
	// ================================

	ReferenceSequenceIndex_t refSeqIndex;
	refSeqIndex.resize(mNumRefSeqs);
	ReferenceSequenceIndex_t::iterator rsiIter = refSeqIndex.begin();

	unsigned int currentRefSeq = 0;
	TempReferenceFilenames_t::const_iterator trfIter;

	CConsole::Heading();
	printf("\n- assembling hash positions for each reference sequence:\n");
	CConsole::Reset();

	//CProgressBar<unsigned int>::StartThread(&currentRefSeq, 0, mNumRefSeqs, "reference sequences");

	for(trfIter = mTempReferenceFilenames.begin(); trfIter != mTempReferenceFilenames.end(); ++trfIter, ++rsiIter, ++currentRefSeq) {

		// store the begin offset
		rsiIter->BeginOffset = ftell64(mJumpOut);

		printf("\nDEBUG: reference index: %u\n", currentRefSeq);

		// ---------------------------------------
		// open all of the temporary sorting files
		// ---------------------------------------

		FileHandles_t tempHandles;
		const unsigned int numTempFiles = (unsigned int)trfIter->size();
		tempHandles.resize(numTempFiles);

		TempFilenames_t::const_iterator tfIter = trfIter->begin();
		FileHandles_t::iterator fhIter;

		for(fhIter = tempHandles.begin(); fhIter != tempHandles.end(); ++fhIter, ++tfIter) {
			*fhIter = fopen(tfIter->c_str(), "rb");
			if(!*fhIter) {
				printf("ERROR: Unable to open temporary file (%s) for reading.\n", tfIter->c_str());
				exit(1);
			}
		}

		// ---------------
		// get the top row
		// ---------------

		HashPositions_t sameHash;
		HashPositions_t topRow;
		topRow.reserve(numTempFiles);

		unsigned int owner = 0;
		for(fhIter = tempHandles.begin(); fhIter != tempHandles.end(); ++fhIter, ++owner) {
			HashPosition hp;
			hp.Owner = owner;
			if(hp.Deserialize(*fhIter)) topRow.push_back(hp);
		}

		sort(topRow.begin(), topRow.end(), SortHashPositionDesc());

		// -----------------
		// process the files
		// -----------------

		while(true) {

			// no more hash positions
			if(topRow.empty()) break;
			HashPosition bestPosition = topRow.back();
			topRow.pop_back();

			if(sameHash.empty()) {

				sameHash.push_back(bestPosition);

			} else {

				if(bestPosition.Hash == sameHash[0].Hash) {
					sameHash.push_back(bestPosition);
				} else {
					++rsiIter->NumHashes;
					StoreHash(sameHash);
					sameHash.clear();
					sameHash.push_back(bestPosition);
				}
			}

			// get the next hash position from the appropriate temp file
			{
				HashPosition hp;
				hp.Owner = bestPosition.Owner;
				if(hp.Deserialize(tempHandles[bestPosition.Owner])) topRow.push_back(hp);
			}

			sort(topRow.begin(), topRow.end(), SortHashPositionDesc());
		}

		// store the last hash
		++rsiIter->NumHashes;
		StoreHash(sameHash);

		// store the end offset
		rsiIter->EndOffset = ftell64(mJumpOut);
	}

	//CProgressBar<unsigned int>::WaitThread();

	// ==================================
	// write the reference sequence index
	// ==================================

	fseek(mJumpOut, REFERENCE_SEQUENCE_INDEX_OFFSET, SEEK_SET);

	printf("\n");
	unsigned int index = 0;
	for(rsiIter = refSeqIndex.begin(); rsiIter != refSeqIndex.end(); ++rsiIter, ++index) {
		const off_type length = rsiIter->EndOffset - rsiIter->BeginOffset;
		printf("INDEX: %u, %llu - %llu: length: %llu, num hashes: %u\n", index, (unsigned long long)rsiIter->BeginOffset, (unsigned long long)rsiIter->EndOffset, (unsigned long long)length, rsiIter->NumHashes);
		fwrite((char*)&rsiIter->BeginOffset, SIZEOF_UINT64, 1, mJumpOut);
		fwrite((char*)&length,               SIZEOF_UINT64, 1, mJumpOut);
		fwrite((char*)&rsiIter->NumHashes,   SIZEOF_INT,    1, mJumpOut);
	}
}

// creates the hash for a supplied fragment
void CJumpCreator::CreateHash(const char* fragment, const unsigned char fragmentLen, uint64_t& key) {

	// set the key to zero
	key = 0;
	const char translation[26] = { 0, 3, 1, 3, -1, -1, 2, 3, -1, -1, 3, -1, 0, 3, -1, -1, -1, 0, 2, 3, -1, 0, 3, 1, 3, -1 };

	if(fragmentLen > 32) {
		cout << "ERROR: This hash table can only handle fragments smaller or equal to 32 bases." << endl;
		exit(1);
	}	

	// convert each nucleotide to its 2-bit representation
	for(unsigned char i = 0; i < fragmentLen; i++) {

		// convert [A,C,G,T] to [0,1,2,3]
		char tValue = translation[fragment[i] - 'A'];

		// catch any unrecognized nucleotides
		if(tValue < 0) {
			cout << "ERROR: Unrecognized nucleotide in hash table: " << fragment[i] << endl;
			cout << "- fragment: ";
			for(unsigned j = 0; j < fragmentLen; j++) cout << fragment[j];
			cout << endl;
			exit(1);
		}

		// shift the key and add the new value
		key = key << 2 | tValue;
	}
}

// hashes the reference sequence and stores the results in sorted temporary files
void CJumpCreator::HashReferenceSequence(const string& bases) {

	// initialize
	const char* pReference = (const char*)bases.data();
	const unsigned int referenceLength = (unsigned int)bases.size();
	const unsigned int maxPositions = referenceLength - mHashSize + 1;

	// hash each base in the reference
	for(unsigned int i = 0; i < maxPositions; ++i, ++pReference) {
		
		bool skipHash = false;
		for(unsigned int j = 0; j < mHashSize; ++j) {
			char refBase = *(pReference + j);
			if((refBase == 'J') || (refBase == 'X') || (refBase == 'N')) {
				skipHash = true;
				break;
			}
		}

		if(skipHash) continue;

		HashPosition hp;
		hp.Position = i;
		CreateHash(pReference, mHashSize, hp.Hash);
		mHashPositions.push_back(hp);

		// dump our hash positions
		if((int)mHashPositions.size() >= mMaxHashPositions) SerializeHashPositions();
	}

	SerializeHashPositions();
}

// hashes the reference sequences in the specified reference sequence archive
void CJumpCreator::HashReferenceSequences(const string& referenceFilename) {

	// open the reference sequence archive
	MosaikReadFormat::CReferenceSequenceReader refReader;
	refReader.Open(referenceFilename);

	// initialize
	vector<ReferenceSequence> refSeqs;
	vector<ReferenceSequence>::const_iterator rsIter;
	refReader.GetReferenceSequences(refSeqs);
	string bases;

	mNumRefSeqs = (unsigned int)refSeqs.size();
	unsigned int currentRefSeq = 0;

	CConsole::Heading();
	printf("- hashing reference sequences:\n");
	CConsole::Reset();

	CProgressBar<unsigned int>::StartThread(&currentRefSeq, 0, mNumRefSeqs, "reference sequences");

	// hash each reference sequence
	for(rsIter = refSeqs.begin(); rsIter != refSeqs.end(); ++rsIter, ++currentRefSeq) {
		TempFilenames_t tempFiles;
		mTempReferenceFilenames.push_back(tempFiles);
		refReader.GetReferenceSequence(rsIter->Name, bases);
		HashReferenceSequence(bases);
	}

	CProgressBar<unsigned int>::WaitThread();
}

// serializes the hash positions to temporary files
void CJumpCreator::SerializeHashPositions(void) {

	// retrieve a temporary filename
	string tempFilename;
	CFileUtilities::GetTempFilename(tempFilename);

	TempReferenceFilenames_t::iterator trfIter = mTempReferenceFilenames.end() - 1;
	trfIter->push_back(tempFilename);
	
	// open the temporary file
	FILE* temp = fopen(tempFilename.c_str(), "wb");

	if(!temp) {
		cout << "ERROR: Unable to open temporary file (" << tempFilename << ") for writing." << endl;
		exit(1);
	}

	// sort the hash positions
	sort(mHashPositions.begin(), mHashPositions.end(), SortHashPositionAsc());

	// serialize
	HashPositions_t::iterator hpIter;
	for(hpIter = mHashPositions.begin(); hpIter != mHashPositions.end(); ++hpIter) hpIter->Serialize(temp);

	// close the temporary file
	fclose(temp);

	// clear the vector
	mHashPositions.clear();
}

// stores the supplied hash positions in the jump database
void CJumpCreator::StoreHash(HashPositions_t& hashPositions) {
	
	int numHashes = (int)hashPositions.size();

	// DEBUG
	//printf("DEBUG: num hashes: %u\n", numHashes);

	if(numHashes == 0) {
		printf("ERROR: Tried to store an empty hash in StoreHash.\n");
		exit(1);
	}

	// limit the number of hashes that will be written
	if(mLimitPositions && (numHashes > mMaxHashPositions)) numHashes = mMaxHashPositions;

	// shuffle the vector
	random_shuffle(hashPositions.begin(), hashPositions.end());

	// localize the hash
	uint64_t hash = hashPositions[0].Hash;

	// write the hash positions
	unsigned int entrySize = (numHashes + 1) * SIZEOF_INT + SIZEOF_UINT64;
	CheckBufferSize(entrySize);
	char* pBuffer = (char*)mBuffer.data();

	unsigned int bufferOffset = 0;

	memcpy(pBuffer + bufferOffset, (char*)&hash, SIZEOF_UINT64);
	bufferOffset += SIZEOF_UINT64;

	memcpy(pBuffer + bufferOffset, (char*)&numHashes, SIZEOF_INT);
	bufferOffset += SIZEOF_INT;

	// DEBUG
	printf("DEBUG: hash: %llu, numHashes: %u\n", (unsigned long long)hash, numHashes);

	for(int i = 0; i < numHashes; i++) {
		memcpy(pBuffer + bufferOffset, (char*)&hashPositions[i].Position, SIZEOF_INT);
		bufferOffset += SIZEOF_INT;
	}

	fwrite(pBuffer, bufferOffset, 1, mJumpOut);
}
