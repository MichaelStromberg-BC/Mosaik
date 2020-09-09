// ***************************************************************************
// CColorspaceUtilities - conversion to colorspace from basespace & vice versa
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "ColorspaceUtilities.h"

// constructor
CColorspaceUtilities::CColorspaceUtilities(void) 
: mNumGapsObserved(NULL)
, mNumGapsObservedLen(0)
, mpBsRefSeqs(NULL)
{
	// initialize our sets and maps
	InitializeAmbiguitySet();
	InitializeBasespaceMap();
	InitializeColorspaceMap();
}

// destructor
CColorspaceUtilities::~CColorspaceUtilities(void) {
	mNumGapsObserved = 0;
	delete [] mNumGapsObserved;
}

// converts the supplied alignment from colorspace to basespace
void CColorspaceUtilities::ConvertAlignmentToBasespace(Alignment& al) {

	// convert the alignment to character arrays
	const unsigned int pairwiseLen = al.Reference.Length();
	char* pReference = al.Reference.Data();
	char* pQuery     = al.Query.Data();

	char* pBSRef     = mpBsRefSeqs[al.ReferenceIndex] + al.ReferenceBegin;

	string bsRefMismatchRegion, bsQryMismatchRegion, bsRefMatchRegion, bsQryMatchRegion;
	ostringstream refSB, qrySB;

	// update the observed gaps array and patches colorspace gaps if needed
	if(UpdateObservedGaps(pReference, pQuery, pairwiseLen)) PatchColorspaceGaps(al);

	// record all of the regions where the reference and query are identical
	RegionVector rv;
	FindIndenticalRegions(pReference, pQuery, pairwiseLen, rv);

	// adjust the alignment if the entire alignment matches the reference
	if(rv.size() == 1) {

		// initialize
		const unsigned short matchLength = rv[0].Length;
		const unsigned short matchBegin  = rv[0].Begin;
		const unsigned short matchEnd    = matchBegin + matchLength - 1;

		// NOTE: we can have two mismatch regions flanking this match region
		GetMatchRegion(bsRefMatchRegion, bsQryMatchRegion, matchBegin, matchEnd, al);

		// handle the left mismatch region
		string bsRefLeftMismatchRegion, bsQryLeftMismatchRegion;
		bool hasLeftMismatchRegion = false;
		if(matchBegin > 0) {
			GetMismatchRegion(bsRefLeftMismatchRegion, bsQryLeftMismatchRegion, 0, matchBegin - 1, false, al);
			hasLeftMismatchRegion = true;
		}

		// handle the right mismatch region
		string bsRefRightMismatchRegion, bsQryRightMismatchRegion;
		bool hasRightMismatchRegion = false;
		if(matchEnd < (pairwiseLen - 1)) {
			GetMismatchRegion(bsRefRightMismatchRegion, bsQryRightMismatchRegion, matchEnd + 1, pairwiseLen - 1, true, al);
			hasRightMismatchRegion = true;
		}

		if(!hasLeftMismatchRegion && !hasRightMismatchRegion) {

			al.Reference = bsRefMatchRegion.c_str();
			al.Query     = bsQryMatchRegion.c_str();

		} else {

			refSB << bsRefLeftMismatchRegion << bsRefMatchRegion << bsRefRightMismatchRegion;
			qrySB << bsQryLeftMismatchRegion << bsQryMatchRegion << bsQryRightMismatchRegion;

			// assign the reference and query in the alignment structure
			al.Reference = refSB.str().c_str();
			al.Query     = qrySB.str().c_str();
		}

	} else {

		RegionVector::const_iterator prevIter = rv.begin();
		RegionVector::const_iterator regIter  = prevIter + 1;

		for(; regIter != rv.end(); ++regIter, ++prevIter) {

			// retrieve the reference and query for the mismatch region
			const unsigned short numMismatchBases = GetMismatchRegion(bsRefMismatchRegion, bsQryMismatchRegion, 
				prevIter->Begin + prevIter->Length, regIter->Begin - 1, false, al);

			// retrieve the reference and query for the match region
			GetMatchRegion(bsRefMatchRegion, bsQryMatchRegion, prevIter->Begin, 
				prevIter->Begin + prevIter->Length - 1, al);

			refSB << bsRefMatchRegion;
			qrySB << bsQryMatchRegion;

			if(numMismatchBases > 0) {
				refSB << bsRefMismatchRegion;
				qrySB << bsQryMismatchRegion;
			}		
		}

		// retrieve the reference and query for the match region
		GetMatchRegion(bsRefMatchRegion, bsQryMatchRegion, prevIter->Begin, 
			prevIter->Begin + prevIter->Length - 1, al);

		refSB << bsRefMatchRegion;
		qrySB << bsQryMatchRegion;

		// assign the reference and query in the alignment structure
		al.Reference = refSB.str().c_str();
		al.Query     = qrySB.str().c_str();
	}

	// -------------------------------
	// update the number of mismatches
	// -------------------------------

	const unsigned int bsPairwiseLen = al.Reference.Length();

	al.NumMismatches = 0;
	Ambiguity_SET_t::const_iterator amIter;
	for(unsigned short i = 0; i < bsPairwiseLen; ++i) {
		amIter = mAmbiguitySet.find(PACK_SORTED_SHORT(pReference[i], pQuery[i]));
		if(amIter == mAmbiguitySet.end()) al.NumMismatches++;
	}

	// ------------------------------------------------------------------------------------------
	// convert the colorspace transition qualities to base qualities
	// NOTE: this algorithm will simply take the minimum of the two qualities that overlap a base
	// ------------------------------------------------------------------------------------------

	const unsigned short numColorspaceQualities = al.BaseQualities.Length();
	const unsigned short lastCSQIndex = numColorspaceQualities - 1;

	CMosaikString csQualities = al.BaseQualities;
	al.BaseQualities.Reserve(numColorspaceQualities + 1);
	al.BaseQualities.SetLength(numColorspaceQualities + 1);

	const char* pCSQual = csQualities.CData();
	char* pBSQual       = al.BaseQualities.Data();

	// handle the first base quality
	*pBSQual = *pCSQual;
	++pBSQual;

	// handle the internal base qualities
	for(unsigned short i = 1; i < numColorspaceQualities; ++i, ++pBSQual) {
		*pBSQual = min(pCSQual[i - 1], pCSQual[i]);
	}

	// handle the final base quality
	*pBSQual = pCSQual[lastCSQIndex];

	// adjust the endpoints (basespace endpoint = colorspace + 1)
	++al.ReferenceEnd;
	++al.QueryEnd;
	al.QueryLength = al.QueryEnd - al.QueryBegin + 1;
}

// converts a colorspace sequence with provided seed base into basespace
void CColorspaceUtilities::ConvertColorspaceToBasespace(char seed, const string& colorspaceSeq, string& basespaceSeq) {

	// make the basespace sequence just as long as the colorspace sequence
	const unsigned short csLen = colorspaceSeq.size();
	basespaceSeq.resize(csLen);

	// create traversal pointers to our strings
	const char* pCS = colorspaceSeq.data();
	char* pBS       = (char*)basespaceSeq.data();

	// convert each colorspace/seed combo into a basespace nucleotide
	BS_MAP_t::const_iterator bsIter;
	for(unsigned int i = 0; i < csLen; ++i, ++pCS, ++pBS) {

		// find the appropriate seed/colorspace transition combination
		bsIter = mBSMap.find(PACK_SHORT(seed, *pCS));
		if(bsIter == mBSMap.end()) {
			printf("ERROR: Unknown combination found when converting to basespace: [%c] & [%c]\n", seed, *pCS);
			exit(1);
		}

		seed = *pBS = bsIter->second;
	}
}

// converts the supplied read from basespace to pseudo-colorspace
void CColorspaceUtilities::ConvertReadBasespaceToPseudoColorspace(CMosaikString& s) {

	char* pPrev   = s.Data();
	char* pString = pPrev + 1;

	// simplify various ambiguity codes
	*pPrev = GetSimplifiedBase(*pPrev);

	CS_MAP_t::const_iterator csIter;
	for(unsigned int i = 1; i < s.Length(); ++i, ++pString, ++pPrev) {

		// simplify various ambiguity codes
		*pString = GetSimplifiedBase(*pString);

		csIter = mCSMap.find(PACK_SHORT(*pPrev, *pString));
		if(csIter == mCSMap.end()) {
			printf("ERROR: Unknown combination found when converting to colorspace: [%c] & [%c]\n", *pPrev, *pString);
			exit(1);
		}

		*pPrev = csIter->second;
	}

	// adjust the read
	s.TrimEnd(1);
}

// converts the supplied read from colorspace to pseudo-colorspace
void CColorspaceUtilities::ConvertReadColorspaceToPseudoColorspace(CMosaikString& s) {
	char* pBases = s.Data();
	for(unsigned int i = 0; i < s.Length(); ++i, ++pBases) {
		switch(*pBases) {
			case '0':
				*pBases = 'A';
				break;
			case '1':
				*pBases = 'C';
				break;
			case '2':
				*pBases = 'G';
				break;
			case '3':
				*pBases = 'T';
				break;
			case 'X':
				break;
			case '-':
				*pBases = 'N';
				break;
			case '.':
				// here we pick an arbitrary colorspace transition, this will have at
				// least 25 % of being correct as opposed to specifying an 'N'.
				*pBases = 'A';
				break;
			default:
				printf("ERROR: Unrecognized nucleotide (%c) when converting read to pseudo-colorspace.\n", pBases[i]);
				exit(1);
				break;
		}
	}
}

// converts the supplied read from pseudo-colorspace to colorspace
void CColorspaceUtilities::ConvertReadPseudoColorspaceToColorspace(CMosaikString& s) {
	char* pBases = s.Data();
	for(unsigned int i = 0; i < s.Length(); ++i, ++pBases) {
		switch(*pBases) {
			case 'A':
				*pBases = '0';
				break;
			case 'C':
				*pBases = '1';
				break;
			case 'G':
				*pBases = '2';
				break;
			case 'T':
				*pBases = '3';
				break;
			case 'X':
			case 'N':
				break;
			default:
				printf("ERROR: Unrecognized nucleotide (%c) when converting read to colorspace.\n", pBases[i]);
				exit(1);
				break;
		}
	}
}

// derives the base at a specified location (meant to handle IUPAC ambiguity codes gracefully)
char CColorspaceUtilities::DeriveBase(Alignment& al, const unsigned short pos) {

	// initialize
	const unsigned short lastIndex = al.Reference.Length() - 1;

	const char* pCSRef = al.Reference.CData();
	const char* pCSQry = al.Query.CData();

	const char* pBSRef = mpBsRefSeqs[al.ReferenceIndex] + al.ReferenceBegin;

	// ===================================================
	// special case: deriving bases that begin at region 0
	// ===================================================

	if(pos == 0) return pBSRef[0];

	// =====================================================
	// normal case: extract the base at the current position
	// =====================================================

	char canonicalBase = pBSRef[pos - mNumGapsObserved[pos]];
	if(IsCanonicalBase(canonicalBase)) return canonicalBase;

	// =====================================================================
	// case 1: find the canonical base by searching earlier in the alignment
	// =====================================================================

	bool hasCanonicalBase = false;
	unsigned int checkIndex = pos - 1;
	char seed;

	while(true) {

		// check if the current base is canonical
		seed = pBSRef[checkIndex - mNumGapsObserved[checkIndex]];
		if(IsCanonicalBase(seed)) {
			hasCanonicalBase = true;
			break;
		}

		// decrement the check index
		if(checkIndex == 0) break;
		--checkIndex;
	}

	// move forward to the original index
	if(hasCanonicalBase) {
		canonicalBase = seed;
		while(checkIndex < pos) {
			canonicalBase = mBSMap[PACK_SHORT(canonicalBase, GetSimplifiedBase(pCSQry[checkIndex]))];
			++checkIndex;
		}

		return canonicalBase;
	}

	// ===================================================================
	// case 2: find the canonical base by searching later in the alignment
	// ===================================================================

	if(pos < lastIndex) {
		checkIndex = pos + 1;

		while(true) {

			// check if the current base is canonical
			seed = pBSRef[checkIndex - mNumGapsObserved[checkIndex]];
			const char csBase = pCSQry[checkIndex - 1 - mNumGapsObserved[checkIndex - 1]];

			if(IsCanonicalBase(seed) && (csBase != '-')) {
				hasCanonicalBase = true;
				break;
			}

			// increment the check index
			if(checkIndex == lastIndex) break;
			++checkIndex;
		}

		// move backward to the original index
		if(hasCanonicalBase) {
			canonicalBase = seed;
			while(checkIndex > pos) {
				const char csBase = pCSQry[checkIndex - 1 - mNumGapsObserved[checkIndex - 1]];
				if(csBase != '-') canonicalBase = mBSMap[PACK_SHORT(canonicalBase, GetSimplifiedBase(csBase))];
				--checkIndex;
			}

			return canonicalBase;
		}
	}

	// =================================================================
	// final case: if we cannot derive the base, we have a major problem
	// =================================================================

	printf("ERROR: Unable to derive the base at position [%u]\n", pos);
	exit(1);
	return canonicalBase;
}

// derives the base for the query at the beginning of an alignment (handles IUPAC ambiguity codes)
char CColorspaceUtilities::DeriveBeginQueryBase(Alignment& al) {

	// initialize
	const unsigned short lastIndex = al.Reference.Length() - 1;

	const char* pCSRef = al.Reference.CData();
	const char* pCSQry = al.Query.CData();

	const char* pBSRef = mpBsRefSeqs[al.ReferenceIndex] + al.ReferenceBegin;

	// =====================================================
	// normal case: extract the base at the current position
	// =====================================================

	char canonicalBase = pBSRef[0];
	if(IsCanonicalBase(canonicalBase)) return canonicalBase;

	// =========================================================================
	// special case: find the canonical base by searching later in the alignment
	// =========================================================================

	unsigned short checkIndex = 1;
	char seed;
	bool hasCanonicalBase = false;

	while(true) {

		// check if the current base is canonical
		seed = pBSRef[checkIndex - mNumGapsObserved[checkIndex]];
		const char csBase = pCSQry[checkIndex - 1 - mNumGapsObserved[checkIndex - 1]];

		if(IsCanonicalBase(seed) && IsCanonicalBase(csBase)) {
			hasCanonicalBase = true;
			break;
		}

		// increment the check index
		if(checkIndex == lastIndex) break;
		++checkIndex;
	}

	// move backward to the original index
	if(hasCanonicalBase) {
		canonicalBase = seed;
		while(checkIndex > 0) {
			const char csBase = pCSQry[checkIndex - 1 - mNumGapsObserved[checkIndex - 1]];
			if(csBase != '-') canonicalBase = mBSMap[PACK_SHORT(canonicalBase, GetSimplifiedBase(csBase))];
			--checkIndex;
		}

		return canonicalBase;
	}

	// =================================================================
	// final case: if we cannot derive the base, we have a major problem
	// =================================================================

	printf("ERROR: Unable to derive the base at the beginning of the colorspace alignment.\n");
	exit(1);
	return canonicalBase;
}

// records regions of contiguous identity in the alignment
void CColorspaceUtilities::FindIndenticalRegions(char* pReference, char* pQuery, const unsigned short pairwiseLen, RegionVector& rv) {

	for(unsigned short i = 0; i < pairwiseLen; i++) {
		if(pReference[i] == pQuery[i]) {
			RegionT r(i);
			unsigned short end = i;
			while((end < pairwiseLen) && (pReference[end] == pQuery[end])) {
				++r.Length;
				++end;
			}
			rv.push_back(r);
			i = end;
		}
	}
}

// retrieves the reference and query bases (basespace) for the specified match region
void CColorspaceUtilities::GetMatchRegion(string& bsRef, string& bsQry, const unsigned short regionBegin, const unsigned short regionEnd, Alignment& al) {

	// initialize
	char* pReference = mpBsRefSeqs[al.ReferenceIndex] + al.ReferenceBegin;
	const unsigned short regionLength = regionEnd - regionBegin + 1;

	// check for non-canonical bases in our matching region
	bool hasNonCanonicalBase = false;
	for(unsigned int pos = regionBegin; pos <= regionEnd; pos++) {
		if(!IsCanonicalBase(pReference[pos])) {
			hasNonCanonicalBase = true;
			break;
		}
	}

	// copy the basespace reference associated with the matching region
	bsRef.resize(regionLength + 1);
	char* pBsRef = (char*)bsRef.data();
	memcpy(pBsRef, pReference + regionBegin - mNumGapsObserved[regionBegin], regionLength + 1);

	// use the transitions to extract the basespace query associated with the matching region
	if(hasNonCanonicalBase) {

		string bsQryMatchBases;
		string csQryMatchBases = al.Query.Substring(regionBegin, regionLength);
		const char matchSeed = DeriveBase(al, regionBegin);
		ConvertColorspaceToBasespace(matchSeed, csQryMatchBases, bsQryMatchBases);

		const unsigned short numQueryBases = bsQryMatchBases.size();
		bsQry.resize(numQueryBases + 1);
		char* pBsQry = (char*)bsQry.data();
		*pBsQry = matchSeed;
		++pBsQry;
		memcpy(pBsQry, bsQryMatchBases.data(), numQueryBases);

	} else {
		bsQry = bsRef;
	}
}

// retrieves the reference and query bases (basespace) for the specified mismatch region
unsigned short CColorspaceUtilities::GetMismatchRegion(string& bsRef, string& bsQry, const unsigned short regionBegin, const unsigned short regionEnd, const bool isLastRegion, Alignment& al) {

	// initialize
	char* pReference = mpBsRefSeqs[al.ReferenceIndex] + al.ReferenceBegin;
	const unsigned short regionLength = regionEnd - regionBegin + 1;

	AlignRegionT csMismatch;
	csMismatch.Reference = al.Reference.Substring(regionBegin, regionLength);
	csMismatch.Query     = al.Query.Substring(regionBegin,     regionLength);

	// convert the mismatch reference region to basespace
	csMismatch.Seed = DeriveBase(al, regionBegin);
	ConvertColorspaceToBasespace(csMismatch.Seed, csMismatch.Reference, bsRef);
	if(regionBegin == 0) bsRef.insert(0, 1, csMismatch.Seed);

	// convert the mismatch reference region to basespace
	if((regionBegin == 0) && !IsCanonicalBase(csMismatch.Seed)) csMismatch.Seed = DeriveBeginQueryBase(al);
	ConvertColorspaceToBasespace(csMismatch.Seed, csMismatch.Query, bsQry);
	if(regionBegin == 0) bsQry.insert(0, 1, csMismatch.Seed);

	unsigned short numMismatchBases = bsQry.size();

	// shorten the mismatch region by one if we have a matching region that follows
	if(!isLastRegion) {
		--numMismatchBases;
		bsRef.resize(numMismatchBases);
		bsQry.resize(numMismatchBases);
	}

	return numMismatchBases;
}

// adds the ambiguity code equivalents
void CColorspaceUtilities::InitializeAmbiguitySet(void) {

	// add the normal bases
	mAmbiguitySet.insert(PACK_SORTED_SHORT('A','A'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('C','C'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('G','G'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('T','T'));

	// add the biallelic bases
	mAmbiguitySet.insert(PACK_SORTED_SHORT('M','A'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('M','C'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('R','A'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('R','G'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('W','A'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('W','T'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('S','C'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('S','G'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('Y','C'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('Y','T'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('K','G'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('K','T'));

	// add the triallelic bases
	mAmbiguitySet.insert(PACK_SORTED_SHORT('V','A'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('V','C'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('V','G'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('H','A'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('H','C'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('H','T'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('D','A'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('D','G'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('D','T'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('B','C'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('B','G'));
	mAmbiguitySet.insert(PACK_SORTED_SHORT('B','T'));
}

// adds the colorspace to basespace conversions
void CColorspaceUtilities::InitializeBasespaceMap(void) {

	//                bas clr     bas
	mBSMap[PACK_SHORT('A','A')] = 'A';
	mBSMap[PACK_SHORT('A','C')] = 'C';
	mBSMap[PACK_SHORT('A','G')] = 'G';
	mBSMap[PACK_SHORT('A','T')] = 'T';
	mBSMap[PACK_SHORT('A','E')] = 'N';
	mBSMap[PACK_SHORT('A','F')] = '-';

	mBSMap[PACK_SHORT('C','A')] = 'C';
	mBSMap[PACK_SHORT('C','C')] = 'A';
	mBSMap[PACK_SHORT('C','G')] = 'T';
	mBSMap[PACK_SHORT('C','T')] = 'G';
	mBSMap[PACK_SHORT('C','F')] = 'N';
	mBSMap[PACK_SHORT('C','I')] = '-';

	mBSMap[PACK_SHORT('G','A')] = 'G';
	mBSMap[PACK_SHORT('G','C')] = 'T';
	mBSMap[PACK_SHORT('G','G')] = 'A';
	mBSMap[PACK_SHORT('G','T')] = 'C';
	mBSMap[PACK_SHORT('G','I')] = 'N';
	mBSMap[PACK_SHORT('G','L')] = '-';

	mBSMap[PACK_SHORT('T','A')] = 'T';
	mBSMap[PACK_SHORT('T','C')] = 'G';
	mBSMap[PACK_SHORT('T','G')] = 'C';
	mBSMap[PACK_SHORT('T','T')] = 'A';
	mBSMap[PACK_SHORT('T','L')] = 'N';
	mBSMap[PACK_SHORT('T','O')] = '-';

	// here we choose bases that are not part of the IUPAC ambiguity codes so that
	// these bases will always register as a mismatch during an alignment.
	mBSMap[PACK_SHORT('N','E')] = 'A';
	mBSMap[PACK_SHORT('N','F')] = 'C';
	mBSMap[PACK_SHORT('N','I')] = 'G';
	mBSMap[PACK_SHORT('N','L')] = 'T';
	mBSMap[PACK_SHORT('N','O')] = 'N';
	mBSMap[PACK_SHORT('N','P')] = '-';

	mBSMap[PACK_SHORT('-','F')] = 'A';
	mBSMap[PACK_SHORT('-','I')] = 'C';
	mBSMap[PACK_SHORT('-','L')] = 'G';
	mBSMap[PACK_SHORT('-','O')] = 'T';
	mBSMap[PACK_SHORT('-','P')] = 'N';
	mBSMap[PACK_SHORT('-','E')] = '-';
}

// adds the basespace to colorspace conversions
void CColorspaceUtilities::InitializeColorspaceMap(void) {

	//                bas bas     clr
	mCSMap[PACK_SHORT('A','A')] = 'A';
	mCSMap[PACK_SHORT('A','C')] = 'C';
	mCSMap[PACK_SHORT('A','G')] = 'G';
	mCSMap[PACK_SHORT('A','T')] = 'T';
	mCSMap[PACK_SHORT('A','N')] = 'E';
	mCSMap[PACK_SHORT('A','-')] = 'F';

	mCSMap[PACK_SHORT('C','A')] = 'C';
	mCSMap[PACK_SHORT('C','C')] = 'A';
	mCSMap[PACK_SHORT('C','G')] = 'T';
	mCSMap[PACK_SHORT('C','T')] = 'G';
	mCSMap[PACK_SHORT('C','N')] = 'F';
	mCSMap[PACK_SHORT('C','-')] = 'I';

	mCSMap[PACK_SHORT('G','A')] = 'G';
	mCSMap[PACK_SHORT('G','C')] = 'T';
	mCSMap[PACK_SHORT('G','G')] = 'A';
	mCSMap[PACK_SHORT('G','T')] = 'C';
	mCSMap[PACK_SHORT('G','N')] = 'I';
	mCSMap[PACK_SHORT('G','-')] = 'L';

	mCSMap[PACK_SHORT('T','A')] = 'T';
	mCSMap[PACK_SHORT('T','C')] = 'G';
	mCSMap[PACK_SHORT('T','G')] = 'C';
	mCSMap[PACK_SHORT('T','T')] = 'A';
	mCSMap[PACK_SHORT('T','N')] = 'L';
	mCSMap[PACK_SHORT('T','-')] = 'O';

	// here we choose bases that are not part of the IUPAC ambiguity codes so that
	// these bases will always register as a mismatch during an alignment.
	mCSMap[PACK_SHORT('N','A')] = 'E';
	mCSMap[PACK_SHORT('N','C')] = 'F';
	mCSMap[PACK_SHORT('N','G')] = 'I';
	mCSMap[PACK_SHORT('N','T')] = 'L';
	mCSMap[PACK_SHORT('N','N')] = 'O';
	mCSMap[PACK_SHORT('N','-')] = 'P';

	// since gaps are inserted during alignment, we don't need to consider making
	// these transitions pass-through the alignment routines.
	mCSMap[PACK_SHORT('-','A')] = 'F';
	mCSMap[PACK_SHORT('-','C')] = 'I';
	mCSMap[PACK_SHORT('-','G')] = 'L';
	mCSMap[PACK_SHORT('-','T')] = 'O';
	mCSMap[PACK_SHORT('-','N')] = 'P';
	mCSMap[PACK_SHORT('-','-')] = 'E';
}

// replaces a gap with the appropriate colorspace transitions
void CColorspaceUtilities::PatchColorspaceGap(char* pCS, Alignment& al, unsigned short& currentIndex) {

	// initialize
	const unsigned short lastIndex = al.Reference.Length();

	// find the right gap position
	unsigned int leftGapPos, rightGapPos;
	leftGapPos = rightGapPos = currentIndex;
	while((pCS[rightGapPos] == '-') && (rightGapPos < lastIndex)) ++rightGapPos;

	// determine which base occurs to the left of the gap (basespace)
	const char leftGapBase = DeriveBase(al, leftGapPos);

	// determine which base occurs to the right of the gap (basespace)
	const char rightGapBase = mBSMap[PACK_SHORT(leftGapBase, GetSimplifiedBase(pCS[rightGapPos]))];

	// update the gap
	for(unsigned short i = leftGapPos; i <= rightGapPos; i++) {
		if(i == leftGapPos) {
			pCS[i] = mCSMap[PACK_SHORT(leftGapBase, '-')];
		} else if(i == rightGapPos) {
			pCS[i] = mCSMap[PACK_SHORT('-', rightGapBase)];
		} else {
			pCS[i] = mCSMap[PACK_SHORT('-', '-')];
		}
	}

	currentIndex = rightGapPos;
}

// replaces the gaps in the pairwise alignment with the appropriate colorspace transitions
void CColorspaceUtilities::PatchColorspaceGaps(Alignment& al) {

	// initialize
	const unsigned short pairwiseLen = al.Reference.Length();
	char* pReference = al.Reference.Data();
	char* pQuery     = al.Query.Data();

	const unsigned short lastIndex = pairwiseLen - 1;

	// here we take advantage of the fact that gaps should NEVER occur in the
	// reference and query sequence at the same position
	for(unsigned short i = 0; i < pairwiseLen; ++i) {
		if(pReference[i] == '-') PatchColorspaceGap(pReference, al, i);
		if(pQuery[i]     == '-') PatchColorspaceGap(pQuery,     al, i);
	}
}

// sets the reference sequences
void CColorspaceUtilities::SetReferenceSequences(char** pBsRefSeqs) {
	mpBsRefSeqs = pBsRefSeqs;
}

// updates the observed gaps array and returns true if gaps are found in the alignment
bool CColorspaceUtilities::UpdateObservedGaps(const char* pReference, const char* pQuery, const unsigned short pairwiseLen) {

	// redimension the observed gaps array
	if(pairwiseLen > mNumGapsObservedLen) {
		delete [] mNumGapsObserved;
		mNumGapsObservedLen = pairwiseLen + ARRAY_EXTENSION;
		mNumGapsObserved = new unsigned short[mNumGapsObservedLen];
	}

	// populate the observed gaps array
	unsigned short numReferenceGaps = 0;
	unsigned short numGapPositions  = 0;
	for(unsigned short i = 0; i < pairwiseLen; ++i) {
		mNumGapsObserved[i] = numReferenceGaps;
		if(pReference[i] == '-')                         ++numReferenceGaps;
		if((pReference[i] == '-') || (pQuery[i] == '-')) ++numGapPositions;
	}

	// return true if we have any gapped positions
	return (numGapPositions > 0);
}
