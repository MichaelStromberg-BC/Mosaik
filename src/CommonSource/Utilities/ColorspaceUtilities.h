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

using namespace std;

#define PACK_SHORT(a,b) (((a) << 8) | (b))
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

typedef map<short, char> CS_MAP_t;
typedef map<short, char> BS_MAP_t;
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
	// records regions of contiguous identity in the alignment
	void FindIndenticalRegions(char* pReference, char* pQuery, const unsigned short pairwiseLen, RegionVector& rv);
	// returns the simplified version of the IUPAC ambiguity code - based on base frequencies in the human genome 36.3
	static inline char GetSimplifiedBase(char c);
	// adds the colorspace to basespace conversions
	void InitializeBasespaceMap(void);
	// adds the basespace to colorspace conversions
	void InitializeColorspaceMap(void);
	// replaces the gaps in the pairwise alignment with a dibase transition code '4'
	void PatchColorspaceGaps(char* pReference, char* pQuery, const unsigned short pairwiseLen);
	// updates the observed gaps array and returns true if gaps are found in the alignment
	bool UpdateObservedGaps(const char* pReference, const char* pQuery, const unsigned short pairwiseLen);
	// our colorspace to basespace conversion map
	BS_MAP_t mBSMap;
	// our basespace to colorspace conversion map
	CS_MAP_t mCSMap;
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

		case 'N':
			c = 'X';
			break;
	}

	return c;
}
