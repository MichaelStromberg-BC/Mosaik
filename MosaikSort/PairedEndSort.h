// ***************************************************************************
// CPairedEndSort - resolves paired-end reads and creates an alignment archive
//                  sorted by reference sequence position.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Str�mberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <iostream>
#include <list>
#include <set>
#include <map>
#include <string>
#include "sqlite3.h"
#include "AlignmentReader.h"
#include "AlignmentStatus.h"
#include "AlignmentWriter.h"
#include "AlignedRead.h"
#include "ConsoleUtilities.h"
#include "ReadGroup.h"
#include "Mosaik.h"
#include "ProgressBar.h"
#include "ProgressCounter.h"
#include "SequencingTechnologies.h"
#include "UnorderedMap.h"

using namespace std;

// here we assume that we'll need space for Sanger length reads (reference bases,
// query bases, query base qualities)
static double DEFAULT_CONFIDENCE_INTERVAL = 0.9973;
#define GAP '-'
#define DEFAULT_SHORTEST_MEDIAN_DISTANCE 1000000000

// define our serialization status flags
#define PE_IS_LONG_READ           1
#define PE_IS_REVERSE_STRAND      2
#define PE_IS_MATE_REVERSE_STRAND 4
#define PE_IS_FIRST_MATE          8
#define PE_IS_RESOLVED_AS_PAIR    16
#define PE_WAS_RESCUED            32

#define SQL_BUFFER_SIZE           10240
#define DUMMY_MODEL               100
#define MODEL_COUNT_THRESHOLD     0.1

// define the results from our multinomial regression
#define UO_INTERCEPT       -11.935959238309490
#define UO_RL_COEFFICIENT    0.359469508790554
#define UO_RL2_COEFFICIENT  -0.002122089928114
#define UO_AQ_COEFFICIENT    1.060903764204119
#define UO_AQ2_COEFFICIENT  -0.006546300635631

#define UU_INTERCEPT        13.378769291414342
#define UU_RL_COEFFICIENT    0.555448038249160
#define UU_RL2_COEFFICIENT  -0.003392797389766
#define UU_AQ_COEFFICIENT    0.527564123421933
#define UU_AQ2_COEFFICIENT  -0.002989447057163

#define UM_INTERCEPT        18.273277656027158
#define UM_RL_COEFFICIENT    0.275353627313707
#define UM_RL2_COEFFICIENT  -0.002483569170412
#define UM_AQ_COEFFICIENT   -0.039141946749387
#define UM_AQ2_COEFFICIENT   0.006464992145754

#define MM_INTERCEPT         3.847409799175512
#define MM_RL_COEFFICIENT    0.247658486576270
#define MM_RL2_COEFFICIENT  -0.002155537232682
#define MM_AQ_COEFFICIENT    0.077801154924765
#define MM_AQ2_COEFFICIENT   0.003795915092381

// define our model count data structure
struct ModelType {
	unsigned char ID;
	unsigned int Count;

	// constructor
	ModelType(void)
		: ID(0)
		, Count(0)
	{}

	// our less-than operator
	bool operator<(const ModelType& mt) const {
		return mt.Count < Count;
	}
};

class CPairedEndSort {
public:
	// constructor
	CPairedEndSort(const unsigned int numCachedReads);
	// destructor
	~CPairedEndSort(void);
	// configures which read pair types should be resolved
	void ConfigureResolution(const bool uo, const bool uu, const bool um, const bool mm);
	// disables fragment alignment quality calculation
	void DisableFragmentAlignmentQuality(void);
	// allows any fragment length when evaluating unique mate-pairs
	void EnableAllUniqueFragmentLengths(void);
	// enables closest mate selection in um read pairs
	void EnableClosestMultipleMateSelection(void);
	// enables consed renaming
	void EnableConsedRenaming(void);
	// enables duplicate read filtering
	void EnableDuplicateFiltering(const string& duplicateDirectory);
	// enables the sampling of all read pairs
	void EnableFullFragmentLengthSampling(void);
	// resolves the paired-end reads found in the specified input file
	void ResolvePairedEndReads(const string& inputFilename, const string& outputFilename);
	// sets the desired confidence interval
	void SetConfidenceInterval(const double& percent);

private:
	// define our sort configuration structure
	struct SortSettings {
		unsigned char AlignmentModel1;
		unsigned char AlignmentModel2;
		string DuplicateDirectory;
		string UnresolvedFilename;
		double ConfidenceInterval;
		unsigned int NumCachedReads;

		SortSettings() 
			: ConfidenceInterval(DEFAULT_CONFIDENCE_INTERVAL)
			, NumCachedReads(0)
		{}
	} mSettings;
	// define our boolean flags structure
	struct FlagData {
		bool AllowAllUniqueFragmentLengths;
		bool FindClosestMultipleMate;
		bool RemoveDuplicates;
		bool RenameMates;
		bool ResolveMM;
		bool ResolveUM;
		bool ResolveUO;
		bool ResolveUU;
		bool SampleAllFragmentLengths;
		bool UseFragmentAlignmentQuality;

		FlagData()
			: AllowAllUniqueFragmentLengths(false)
			, FindClosestMultipleMate(false)
			, RemoveDuplicates(false)
			, RenameMates(false)
			, ResolveMM(false)
			, ResolveUM(false)
			, ResolveUO(false)
			, ResolveUU(false)
			, SampleAllFragmentLengths(false)
			, UseFragmentAlignmentQuality(true)
		{}
	} mFlags;
	// retrieves an alignment from the specified temporary file and adds it to the specified list
	void AddAlignment(FILE* tempFile, const unsigned int owner, list<Alignment>& alignments);
	// exchanges alignment info between two mates
	static inline void ExchangeMateInfo(Alignment* pMate1Al, Alignment* pMate2Al);
	// returns the current alignment model based on the order and orientation of the mates
	static inline unsigned char GetCurrentModel(unsigned int m1Begin, bool m1IsReverseStrand, unsigned int m2Begin, bool m2IsReverseStrand);
	// calculates the fragment alignment quality based on the fragment class
	static inline unsigned char GetFragmentAlignmentQuality(const unsigned char aq, const unsigned int readLength, const bool isUO, const bool isUU, const bool isMM);
	// retrieves an alignment from the specified temporary file
	bool GetAlignment(FILE* tempFile, const unsigned int owner, Alignment& al);
	// records the observed gaps in the specified reference 
	void RecordReferenceGaps(Alignment& al);
	// serializes the specified vector to a temporary file
	uint64_t Serialize(list<Alignment>& alignmentCache, const unsigned int numEntries);
	// our temporary file vector
	vector<string> mTempFiles;
	// our output buffer
	unsigned char* mBuffer;
	unsigned int mBufferLen;
	// our reference gap hash map vector and associated iterator
	vector<unordered_map<unsigned int, unsigned short> > mRefGapVector;
	unordered_map<unsigned int, unsigned short>::iterator mRefGapIter;
};

// returns the current alignment model based on the order and orientation of the mates
inline unsigned char CPairedEndSort::GetCurrentModel(unsigned int m1Begin, bool m1IsReverseStrand, unsigned int m2Begin, bool m2IsReverseStrand) {

	unsigned char currentModel = DUMMY_MODEL;

	if(m1Begin < m2Begin) {

		if(!m1IsReverseStrand && !m2IsReverseStrand) currentModel = 0;
		if(!m1IsReverseStrand && m2IsReverseStrand)  currentModel = 1;
		if(m1IsReverseStrand  && !m2IsReverseStrand) currentModel = 2;
		if(m1IsReverseStrand  && m2IsReverseStrand)  currentModel = 3;

	} else {

		if(!m2IsReverseStrand && !m1IsReverseStrand) currentModel = 4;
		if(!m2IsReverseStrand && m1IsReverseStrand)  currentModel = 5;
		if(m2IsReverseStrand  && !m1IsReverseStrand) currentModel = 6;
		if(m2IsReverseStrand  && m1IsReverseStrand)  currentModel = 7;
	}

	return currentModel;
}

// calculates the fragment alignment quality based on the fragment class
inline unsigned char CPairedEndSort::GetFragmentAlignmentQuality(const unsigned char aq, const unsigned int readLength, const bool isUO, const bool isUU, const bool isMM) {

	// calculate the new fragment alignment quality
	double actualAQ = 0.0;

	if(isUO) {
		actualAQ = UO_INTERCEPT + UO_RL_COEFFICIENT * readLength + UO_RL2_COEFFICIENT * readLength * readLength + UO_AQ_COEFFICIENT * aq + UO_AQ2_COEFFICIENT * aq * aq;
	} else if(isUU) {
		actualAQ = UU_INTERCEPT + UU_RL_COEFFICIENT * readLength + UU_RL2_COEFFICIENT * readLength * readLength + UU_AQ_COEFFICIENT * aq + UU_AQ2_COEFFICIENT * aq * aq;
	} else if(isMM) {
		actualAQ = MM_INTERCEPT + MM_RL_COEFFICIENT * readLength + MM_RL2_COEFFICIENT * readLength * readLength + MM_AQ_COEFFICIENT * aq + MM_AQ2_COEFFICIENT * aq * aq;
	} else {
		actualAQ = UM_INTERCEPT + UM_RL_COEFFICIENT * readLength + UM_RL2_COEFFICIENT * readLength * readLength + UM_AQ_COEFFICIENT * aq + UM_AQ2_COEFFICIENT * aq * aq;
	}

	// adjust the alignment quality
	if(actualAQ < 0.0)  actualAQ = 0.0;
	if(actualAQ > 99.0) actualAQ = 99.0;
	return (unsigned char)actualAQ;
}

// exchanges alignment info between two mates
inline void CPairedEndSort::ExchangeMateInfo(Alignment* pMate1Al, Alignment* pMate2Al) {
	pMate1Al->MateReferenceIndex  = pMate2Al->ReferenceIndex;
	pMate1Al->MateReferenceBegin  = pMate2Al->ReferenceBegin;
	pMate1Al->MateReferenceEnd    = pMate2Al->ReferenceEnd;
	pMate1Al->IsMateReverseStrand = pMate2Al->IsReverseStrand;
	pMate1Al->IsResolvedAsPair    = true;

	pMate2Al->MateReferenceIndex  = pMate1Al->ReferenceIndex;
	pMate2Al->MateReferenceBegin  = pMate1Al->ReferenceBegin;
	pMate2Al->MateReferenceEnd    = pMate1Al->ReferenceEnd;
	pMate2Al->IsMateReverseStrand = pMate1Al->IsReverseStrand;
	pMate2Al->IsResolvedAsPair    = true;
}