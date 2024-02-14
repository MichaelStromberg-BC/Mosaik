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
	// initialize our basespace and colorspace conversion maps
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

	// update the observed gaps array and patches colorspace gaps if needed
	if(UpdateObservedGaps(pReference, pQuery, pairwiseLen)) {
		PatchColorspaceGaps(pReference, pQuery, pairwiseLen);
	}

	// record all of the regions where the reference and query are identical
	RegionVector rv;
	FindIndenticalRegions(pReference, pQuery, pairwiseLen, rv);

	// adjust the alignment if the entire alignment matches the reference
	if(rv.size() == 1) {

		al.Reference.Copy(mpBsRefSeqs[al.ReferenceIndex] + rv[0].Begin + al.ReferenceBegin, rv[0].Length + 1);
		al.Query = al.Reference;

	} else {

		// resolve the mismatched regions
		ostringstream refSB, qrySB;
		RegionVector::const_iterator prevIter = rv.begin();
		RegionVector::const_iterator regIter = prevIter + 1;

		for(; regIter != rv.end(); ++regIter, ++prevIter) {

			// extract the mismatch region
			const unsigned short mmBegin = prevIter->Begin + prevIter->Length;
			const unsigned short mmEnd   = regIter->Begin - 1;

			string bsRefMismatchBases, bsQryMismatchBases;
			string csRefMismatchBases = al.Reference.Substring(mmBegin, mmEnd - mmBegin + 1);
			string csQryMismatchBases = al.Query.Substring(mmBegin,     mmEnd - mmBegin + 1);

			// convert the mismatch region to basespace
			const char lastBase = mpBsRefSeqs[al.ReferenceIndex][al.ReferenceBegin + mmBegin - mNumGapsObserved[mmBegin]];
			ConvertColorspaceToBasespace(lastBase, csQryMismatchBases, bsQryMismatchBases);
			ConvertColorspaceToBasespace(lastBase, csRefMismatchBases, bsRefMismatchBases);
			const unsigned short numMismatchBases = bsQryMismatchBases.size() - 1;

			// extract the preceeding identity region
			CMosaikString matchRegion;
			matchRegion.Copy(mpBsRefSeqs[al.ReferenceIndex] + prevIter->Begin + al.ReferenceBegin, prevIter->Length + 1);

			// add the basespace regions to our string builder
			if(numMismatchBases == 0) {
				refSB << matchRegion;
				qrySB << matchRegion;
			} else {
				refSB << matchRegion << bsRefMismatchBases.substr(0, numMismatchBases);
				qrySB << matchRegion << bsQryMismatchBases.substr(0, numMismatchBases);
			}
		}

		// assign the reference and query sequence
		CMosaikString matchRegion;
		matchRegion.Copy(mpBsRefSeqs[al.ReferenceIndex] + al.ReferenceBegin + prevIter->Begin - mNumGapsObserved[prevIter->Begin], prevIter->Length + 1);

		refSB << matchRegion;
		qrySB << matchRegion;

		al.Reference = refSB.str().c_str();
		al.Query     = qrySB.str().c_str();
	}

	// update the number of mismatches
	// TODO: This should be augmented to support IUPAC ambiguity codes
	const unsigned int bsPairwiseLen = al.Reference.Length();

	al.NumMismatches = 0;
	for(unsigned short i = 0; i < bsPairwiseLen; ++i) {
		if(pReference[i] != pQuery[i]) al.NumMismatches++;
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
	// TODO: check that this is true even for reverse complement reads
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

// adds the colorspace to basespace conversions
void CColorspaceUtilities::InitializeBasespaceMap(void) {

	mBSMap[PACK_SHORT('A','A')] = 'A';
	mBSMap[PACK_SHORT('A','C')] = 'C';
	mBSMap[PACK_SHORT('A','G')] = 'G';
	mBSMap[PACK_SHORT('A','T')] = 'T';
	mBSMap[PACK_SHORT('A','N')] = '-';
	mBSMap[PACK_SHORT('A','X')] = 'X';

	mBSMap[PACK_SHORT('C','A')] = 'C';
	mBSMap[PACK_SHORT('C','C')] = 'A';
	mBSMap[PACK_SHORT('C','G')] = 'T';
	mBSMap[PACK_SHORT('C','T')] = 'G';
	mBSMap[PACK_SHORT('C','N')] = '-';
	mBSMap[PACK_SHORT('C','X')] = 'X';

	mBSMap[PACK_SHORT('G','A')] = 'G';
	mBSMap[PACK_SHORT('G','C')] = 'T';
	mBSMap[PACK_SHORT('G','G')] = 'A';
	mBSMap[PACK_SHORT('G','T')] = 'C';
	mBSMap[PACK_SHORT('G','N')] = '-';
	mBSMap[PACK_SHORT('G','X')] = 'X';

	mBSMap[PACK_SHORT('T','A')] = 'T';
	mBSMap[PACK_SHORT('T','C')] = 'G';
	mBSMap[PACK_SHORT('T','G')] = 'C';
	mBSMap[PACK_SHORT('T','T')] = 'A';
	mBSMap[PACK_SHORT('T','N')] = '-';
	mBSMap[PACK_SHORT('T','X')] = 'X';

	mBSMap[PACK_SHORT('-','A')] = 'A';
	mBSMap[PACK_SHORT('-','C')] = 'C';
	mBSMap[PACK_SHORT('-','G')] = 'G';
	mBSMap[PACK_SHORT('-','T')] = 'T';
	mBSMap[PACK_SHORT('-','N')] = '-';
	mBSMap[PACK_SHORT('-','X')] = 'X';

	mBSMap[PACK_SHORT('X','A')] = 'A';
	mBSMap[PACK_SHORT('X','C')] = 'C';
	mBSMap[PACK_SHORT('X','G')] = 'G';
	mBSMap[PACK_SHORT('X','T')] = 'T';
	mBSMap[PACK_SHORT('X','N')] = '-';
	mBSMap[PACK_SHORT('X','X')] = 'X';
}

// adds the basespace to colorspace conversions
void CColorspaceUtilities::InitializeColorspaceMap(void) {

	mCSMap[PACK_SHORT('A','A')] = 'A';
	mCSMap[PACK_SHORT('A','C')] = 'C';
	mCSMap[PACK_SHORT('A','G')] = 'G';
	mCSMap[PACK_SHORT('A','T')] = 'T';
	mCSMap[PACK_SHORT('A','N')] = 'N';
	mCSMap[PACK_SHORT('A','X')] = 'X';

	mCSMap[PACK_SHORT('C','A')] = 'C';
	mCSMap[PACK_SHORT('C','C')] = 'A';
	mCSMap[PACK_SHORT('C','G')] = 'T';
	mCSMap[PACK_SHORT('C','T')] = 'G';
	mCSMap[PACK_SHORT('C','N')] = 'N';
	mCSMap[PACK_SHORT('C','X')] = 'X';

	mCSMap[PACK_SHORT('G','A')] = 'G';
	mCSMap[PACK_SHORT('G','C')] = 'T';
	mCSMap[PACK_SHORT('G','G')] = 'A';
	mCSMap[PACK_SHORT('G','T')] = 'C';
	mCSMap[PACK_SHORT('G','N')] = 'N';
	mCSMap[PACK_SHORT('G','X')] = 'X';

	mCSMap[PACK_SHORT('T','A')] = 'T';
	mCSMap[PACK_SHORT('T','C')] = 'G';
	mCSMap[PACK_SHORT('T','G')] = 'C';
	mCSMap[PACK_SHORT('T','T')] = 'A';
	mCSMap[PACK_SHORT('T','N')] = 'N';
	mCSMap[PACK_SHORT('T','X')] = 'X';

	mCSMap[PACK_SHORT('-','A')] = 'A';
	mCSMap[PACK_SHORT('-','C')] = 'C';
	mCSMap[PACK_SHORT('-','G')] = 'G';
	mCSMap[PACK_SHORT('-','T')] = 'T';
	mCSMap[PACK_SHORT('-','N')] = 'N';
	mCSMap[PACK_SHORT('-','X')] = 'X';

	mCSMap[PACK_SHORT('X','A')] = 'A';
	mCSMap[PACK_SHORT('X','C')] = 'C';
	mCSMap[PACK_SHORT('X','G')] = 'G';
	mCSMap[PACK_SHORT('X','T')] = 'T';
	mCSMap[PACK_SHORT('X','N')] = 'N';
	mCSMap[PACK_SHORT('X','X')] = 'X';
}

// replaces the gaps in the pairwise alignment with a dibase transition code '4'
// NOTE: we arbitrarily add an extra gap transition - this should be modified so
//       that we create a true transition to the next base.
void CColorspaceUtilities::PatchColorspaceGaps(char* pReference, char* pQuery, const unsigned short pairwiseLen) {

	const unsigned int lastIndex = pairwiseLen - 1;

	// here we take advantage of the fact that gaps should NEVER occur in the
	// reference and query sequence at the same position
	for(unsigned short i = 0; i < pairwiseLen; ++i) {

		// patch the gaps in the reference
		if(pReference[i] == '-') {
			unsigned int index = i;
			while((pReference[index] == '-') && (index < lastIndex)) index++;
			const unsigned int length = index - i + 1;
			memset(pReference + i, 'N', length);
			i = index;
		}

		// patch the gaps in the query
		if(pQuery[i] == '-') {
			unsigned int index = i;
			while((pQuery[index] == '-') && (index < lastIndex)) index++;
			const unsigned int length = index - i + 1;
			memset(pQuery + i, 'N', length);
			i = index;
		}
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
