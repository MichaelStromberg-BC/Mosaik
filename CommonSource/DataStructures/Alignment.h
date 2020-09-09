// ***************************************************************************
// Alignment.h - stores everything related to alignments.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include "MosaikString.h"
#include "AlignmentStatus.h"
#include "Mosaik.h"

using namespace std;

// our alignment structure [all members identified as (temp) are not serialized to disk]
struct Alignment {

	unsigned int MateReferenceBegin;   // required for SAM/BAM
	unsigned int MateReferenceEnd;     // required for SAM/BAM
	unsigned int MateReferenceIndex;   // required for SAM/BAM
	unsigned int ReferenceBegin;
	unsigned int ReferenceEnd;
	unsigned int ReferenceIndex;
	unsigned int Owner;                // the temporary file that contains the alignment
	unsigned int ReadGroupCode;        // the read group code (temp)
	unsigned short QueryLength;        // used during filtering (temp)
	unsigned short NumMismatches;      // number of mismatches
	unsigned short QueryBegin;
	unsigned short QueryEnd;
	unsigned char Quality;             // alignment quality
	bool IsFirstMate;                  // is this alignment from the first mate in a paired-end read
	bool IsMateReverseStrand;          // read orientation for the mate
	bool IsPairedEnd;                  // is the read sequenced as a paired-end read
	bool IsResolvedAsPair;             // is the alignment part of resolved paired-end read
	bool IsReverseStrand;              // read orientation
	bool WasRescued;                   // was the alignment rescued during local alignment search
	char* ReferenceName;               // only filled via CAlignmentReader (temp)
	CMosaikString Reference;
	CMosaikString Query;
	CMosaikString BaseQualities;
	CMosaikString Name;              // the read name

	// constructors
	Alignment(void)
		: MateReferenceBegin(0)
		, MateReferenceEnd(0)
		, MateReferenceIndex(0)
		, ReferenceIndex(0)
		, ReadGroupCode(0)
		, Quality(0)
		, IsFirstMate(false)
		, IsMateReverseStrand(false)
		, IsPairedEnd(false)
		, IsResolvedAsPair(false)
		, IsReverseStrand(false)
		, WasRescued(false)
		, ReferenceName(NULL)
	{}

	// our less-than operator
	bool operator<(const Alignment& al) const {
		if(ReferenceIndex == al.ReferenceIndex) return ReferenceBegin < al.ReferenceBegin;
		return ReferenceIndex < al.ReferenceIndex;
	}

	// returns the number of bytes used in the serialized representation
	unsigned int GetSerializationSize(void) const {
		const unsigned int pairwiseLength = Reference.Length();
		const bool isLongRead = ((QueryEnd > 255) || (pairwiseLength > 255) ? true : false);
		return pairwiseLength + BaseQualities.Length() + SIZEOF_INT * 3 + SIZEOF_SHORT + 2 + (isLongRead ? SIZEOF_SHORT * 3 : 3);
	}

	// serializes the current alignment into the supplied buffer (pre-allocated)
	void Serialize(char* &pBuffer) const {

		// store the reference sequence start position
		memcpy(pBuffer, (char*)&ReferenceBegin, SIZEOF_INT);
		pBuffer += SIZEOF_INT;

		// store the reference sequence end position
		memcpy(pBuffer, (char*)&ReferenceEnd, SIZEOF_INT);
		pBuffer += SIZEOF_INT;

		// store the reference sequence index
		memcpy(pBuffer, (char*)&ReferenceIndex, SIZEOF_INT);
		pBuffer += SIZEOF_INT;

		// store the alignment quality
		*pBuffer = Quality; ++pBuffer;

		// store the alignment status flag
		unsigned char status = AF_UNKNOWN;

		if(IsPairedEnd && IsFirstMate)  status |= AF_IS_FIRST_MATE;
		if(IsPairedEnd && !IsFirstMate) status |= AF_IS_SECOND_MATE;
		if(IsReverseStrand)             status |= AF_IS_REVERSE_STRAND;
		if(WasRescued)                  status |= AF_WAS_RESCUED;

		*pBuffer = status; ++pBuffer;

		// store the number of mismatches
		memcpy(pBuffer, (char*)&NumMismatches, SIZEOF_SHORT);
		pBuffer += SIZEOF_SHORT;

		const unsigned short pairwiseLength = (unsigned short)Reference.Length();
		const bool isLongRead = ((QueryEnd > 255) || (pairwiseLength > 255) ? true : false);

		if(isLongRead) {

			// store the pairwise length
			memcpy(pBuffer, (char*)&pairwiseLength, SIZEOF_SHORT);
			pBuffer += SIZEOF_SHORT;

			// store the query begin
			memcpy(pBuffer, (char*)&QueryBegin, SIZEOF_SHORT);
			pBuffer += SIZEOF_SHORT;

			// store the query end
			memcpy(pBuffer, (char*)&QueryEnd, SIZEOF_SHORT);
			pBuffer += SIZEOF_SHORT;

		} else {

			// store the pairwise length
			*pBuffer = (unsigned char)pairwiseLength; ++pBuffer;

			// store the query begin
			*pBuffer = (unsigned char)QueryBegin; ++pBuffer;

			// store the query end
			*pBuffer = (unsigned char)QueryEnd; ++pBuffer;
		}

		// pack the pairwise-alignment string
		CMosaikString packString(Reference);
		packString.Pack(Query);

		// store the packed pairwise alignment
		memcpy(pBuffer, packString.CData(), pairwiseLength);
		pBuffer += pairwiseLength;

		// store the pairwise query base qualities
		const unsigned int bqLength = BaseQualities.Length();
		memcpy(pBuffer, BaseQualities.CData(), bqLength);
		pBuffer += bqLength;
	}
};
