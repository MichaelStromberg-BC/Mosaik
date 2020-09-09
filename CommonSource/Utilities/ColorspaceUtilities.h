// ***************************************************************************
// CColorspaceUtilities - conversion to colorspace from basespace & vice versa
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
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "Alignment.h"
#include "MosaikString.h"
#include "UnorderedSet.h"

using namespace std;

#define PACK_SHORT(a,b) (((a) << 8) | (b))
#define PACK_SORTED_SHORT(a,b) ((a) < (b) ? ((a) << 8) | (b) : ((b) << 8) | (a))
#define ARRAY_EXTENSION 10

// define our region data structure
struct RegionT {
	unsigned short Begin;
	unsigned short Length;

	// constructor
	RegionT(unsigned short beg)
		: Begin(beg)
		, Length(0)
	{}
};

// define our reference and query data structure
struct AlignRegionT {
	string Reference;
	string Query;
	char Seed;
};

typedef map<short, char> CS_MAP_t;
typedef map<short, char> BS_MAP_t;
typedef tr1::unordered_set<short> Ambiguity_SET_t;
typedef vector<RegionT> RegionVector;

class CColorspaceUtilities {
public:
	// constructor
	CColorspaceUtilities(void);
	// destructor
	~CColorspaceUtilities(void);
	// converts the supplied alignment from colorspace to basespace
	void ConvertAlignmentToBasespace(Alignment& al);
	// converts the supplied read from basespace to pseudo-colorspace
	void ConvertReadBasespaceToPseudoColorspace(CMosaikString& s);
	// converts the supplied read from colorspace to pseudo-colorspace
	void ConvertReadColorspaceToPseudoColorspace(CMosaikString& s);
	// converts the supplied read from pseudo-colorspace to colorspace
	void ConvertReadPseudoColorspaceToColorspace(CMosaikString& s);
	// sets the reference sequences
	void SetReferenceSequences(char** pBsRefSeqs);

private:
	// converts a colorspace sequence with provided seed base into basespace
	void ConvertColorspaceToBasespace(char seed, const string& colorspaceSeq, string& basespaceSeq);
	// derives the base at a specified location (meant to handle IUPAC ambiguity codes gracefully)
	char DeriveBase(Alignment& al, const unsigned short pos);
	// derives the base for the query at the beginning of an alignment (handles IUPAC ambiguity codes)
	char DeriveBeginQueryBase(Alignment& al);
	// records regions of contiguous identity in the alignment
	void FindIndenticalRegions(char* pReference, char* pQuery, const unsigned short pairwiseLen, RegionVector& rv);
	// retrieves the reference and query bases (basespace) for the specified match region
	void GetMatchRegion(string& bsRef, string& bsQry, const unsigned short regionBegin, const unsigned short regionEnd, Alignment& al);
	// retrieves the reference and query bases (basespace) for the specified mismatch region
	unsigned short GetMismatchRegion(string& bsRef, string& bsQry, const unsigned short regionBegin, const unsigned short regionEnd, const bool isLastRegion, Alignment& al);
	// returns the simplified version of the IUPAC ambiguity code - based on base frequencies in the human genome 36.3
	static inline char GetSimplifiedBase(char c);
	// adds the ambiguity code equivalents
	void InitializeAmbiguitySet(void);
	// adds the colorspace to basespace conversions
	void InitializeBasespaceMap(void);
	// adds the basespace to colorspace conversions
	void InitializeColorspaceMap(void);
	// returns true if the base is a canonical base { A, C, G, T }
	static inline bool IsCanonicalBase(char c);
	// replaces the gaps in the pairwise alignment with the appropriate colorspace transitions
	void PatchColorspaceGaps(Alignment& al);
	// replaces a gap with the appropriate colorspace transitions
	void PatchColorspaceGap(char* pCS, Alignment& al, unsigned short& currentIndex);
	// updates the observed gaps array and returns true if gaps are found in the alignment
	bool UpdateObservedGaps(const char* pReference, const char* pQuery, const unsigned short pairwiseLen);
	// our colorspace to basespace conversion map
	BS_MAP_t mBSMap;
	// our basespace to colorspace conversion map
	CS_MAP_t mCSMap;
	// our ambiguity code equality set
	Ambiguity_SET_t mAmbiguitySet;
	// stores the number of gaps observed at a particular location
	unsigned short* mNumGapsObserved;
	unsigned short mNumGapsObservedLen;
	// our basespace reference sequence
	char** mpBsRefSeqs;
};

// returns the simplified version of the IUPAC ambiguity code - based on base frequencies in the human genome 36.3
inline char CColorspaceUtilities::GetSimplifiedBase(char c) {

	switch(c) {
		case 'M':
		case 'R':
		case 'V':
			c = 'A';
			break;

		case 'S':
			c = 'G';
			break;

		case 'B':
		case 'D':
		case 'H':
		case 'K':
		case 'W':
		case 'Y':
			c = 'T';
			break;

		case 'X':
			c = 'N';
			break;
	}

	return c;
}

// returns true if the base is a canonical base { A, C, G, T }
inline bool CColorspaceUtilities::IsCanonicalBase(char c) {
	if((c == 'A') || (c == 'C') || (c == 'G') || (c == 'T')) return true;
	return false;
}
