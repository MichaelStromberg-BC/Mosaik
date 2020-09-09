// ***************************************************************************
// AlignedRead.h - stores a collection of mate 1 and mate 2 alignments.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#pragma once

#include <vector>
#include "Alignment.h"
#include "MosaikString.h"

using namespace std;

namespace Mosaik {
	struct AlignedRead {
		unsigned int ReadGroupCode;
		CMosaikString Name;
		vector<Alignment> Mate1Alignments;
		vector<Alignment> Mate2Alignments;
		bool IsLongRead;
		bool IsPairedEnd;

		// constructor
		AlignedRead()
			: ReadGroupCode(0)
			, IsLongRead(false)
			, IsPairedEnd(false)
		{}

		// serializes the current aligned read into the supplied string
		void Serialize(string& buffer) {

			// initialize
			const unsigned int numMate1Alignments = (unsigned int)Mate1Alignments.size();
			const unsigned int numMate2Alignments = (unsigned int)Mate2Alignments.size();

			const bool haveMate1 = (numMate1Alignments > 0 ? true : false);
			const bool haveMate2 = (numMate2Alignments > 0 ? true : false);

			vector<Alignment>::const_iterator alIter;

			// calculate the entry size
			const unsigned char readNameLen = (unsigned char)Name.Length();
			unsigned int entrySize = readNameLen + SIZEOF_INT + 2;

			if(haveMate1) {
				for(alIter = Mate1Alignments.begin(); alIter != Mate1Alignments.end(); ++alIter) entrySize += alIter->GetSerializationSize();
				entrySize += SIZEOF_INT;
			}

			if(haveMate2) {
				for(alIter = Mate2Alignments.begin(); alIter != Mate2Alignments.end(); ++alIter) entrySize += alIter->GetSerializationSize();
				entrySize += SIZEOF_INT;
			}

			// adjust the buffer size
			buffer.resize(entrySize);
			char* pBuffer = (char*)buffer.data();
			//const char* pBegin = pBuffer;

			// derive our read status
			unsigned char readStatus = RF_UNKNOWN;
			if(haveMate1) readStatus |= RF_HAVE_MATE1;
			if(haveMate2) readStatus |= RF_HAVE_MATE2;

			// store the read name
			*pBuffer = readNameLen; ++pBuffer;

			memcpy(pBuffer, Name.CData(), readNameLen);
			pBuffer += readNameLen;

			// store the read group code
			memcpy(pBuffer, (char*)&ReadGroupCode, SIZEOF_INT);
			pBuffer += SIZEOF_INT;

			// store the read status flag
			*pBuffer = readStatus; ++pBuffer;

			// store the number of mate 1 alignments
			if(haveMate1) {
				memcpy(pBuffer, (char*)&numMate1Alignments, SIZEOF_INT);
				pBuffer += SIZEOF_INT;
			}

			// store the number of mate 2 alignments
			if(haveMate2) {
				memcpy(pBuffer, (char*)&numMate2Alignments, SIZEOF_INT);
				pBuffer += SIZEOF_INT;
			}

			// serialize each mate 1 alignment
			if(haveMate1) {
				for(alIter = Mate1Alignments.begin(); alIter != Mate1Alignments.end(); ++alIter) alIter->Serialize(pBuffer);
			}

			// serialize each mate 2 alignment
			if(haveMate2) {
				for(alIter = Mate2Alignments.begin(); alIter != Mate2Alignments.end(); ++alIter) alIter->Serialize(pBuffer);
			}

			//printf("entry size: %u, buffer size: %u\n", entrySize, (pBuffer - pBegin));
		}
	};
}
