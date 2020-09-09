// ***************************************************************************
// CMemoryUtilities - used in monitoring memory usage and dynamic buffer
//                    reallocation.
// ---------------------------------------------------------------------------
// (c) 2006 - 2009 Michael Strömberg
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Dual licenced under the GNU General Public License 2.0+ license or as
// a commercial license with the Marth Lab.
// ***************************************************************************

#include "MemoryUtilities.h"

// checks if the buffer is large enough to accomodate the requested size
void CMemoryUtilities::CheckBufferSize(char* &pBuffer, unsigned int& bufferLen, const unsigned int requestedBytes) {
	try {
		if(requestedBytes > bufferLen) {
			bufferLen = requestedBytes + 10;
			delete [] pBuffer;
			pBuffer = new char[bufferLen];
		}
	} catch(bad_alloc) {
		cout << "ERROR: Out of memory when allocating " << requestedBytes << " bytes." << endl;
		exit(1);
	}
}

// checks if the buffer is large enough to accomodate the requested size
void CMemoryUtilities::CheckBufferSize(unsigned char* &pBuffer, unsigned int& bufferLen, const unsigned int requestedBytes) {
	try {
		if(requestedBytes > bufferLen) {
			bufferLen = requestedBytes + 10;
			delete [] pBuffer;
			pBuffer = new unsigned char[bufferLen];
		}
	} catch(bad_alloc) {
		cout << "ERROR: Out of memory when allocating " << requestedBytes << " bytes." << endl;
		exit(1);
	}
}

// returns the number of bytes used by the current process
uint64_t CMemoryUtilities::GetMemoryUsage(void) {

	uint64_t numBytesUsed = 0;

#ifdef WIN32

	DWORD aProcesses[1024], cbNeeded;

	if(EnumProcesses(aProcesses, sizeof(aProcesses), &cbNeeded)) {
		DWORD numProcesses = cbNeeded / sizeof(DWORD);
		uint64_t processBytesUsed = 0;
		for(unsigned int i = 0; i < numProcesses; ++i) {
			processBytesUsed = GetMemoryUsedByProcess(aProcesses[i]);
			if(processBytesUsed > numBytesUsed) numBytesUsed = processBytesUsed;
		}
	}

#else

	process_id_t pid = getpid();

	char filename[FILENAME_MAX];
	snprintf(filename, FILENAME_MAX, "/proc/%u/stat", pid);
	FILE* in = fopen(filename, "rb");

	if(!in) {
		cout << "ERROR: Unable to get the memory usage for the current process." << endl;
		exit(1);
	}

	// grab the stat information all at once
	char buffer[MEM_USAGE_BUFFER_LEN];
	size_t numBytesRead = fread(buffer, 1, MEM_USAGE_BUFFER_LEN - 1, in);
	fclose(in);

	buffer[numBytesRead] = 0x0;

	// get the resident program size
	long numRssMemoryPages;

	sscanf(buffer, "%*d %*s %*c %*d %*d %*d %*d %*d %*u %*u %*u %*u %*u %*u "
		"%*u %*d %*d %*d %*d %*d %*d %*u %*u %ld", &numRssMemoryPages);

	numBytesUsed = numRssMemoryPages * BYTES_PER_MEMORY_PAGE;

#endif


	return numBytesUsed;
}

// returns the number of bytes used by the specified process
uint64_t CMemoryUtilities::GetMemoryUsedByProcess(const process_id_t& pid) {
#ifdef WIN32
	PROCESS_MEMORY_COUNTERS pmc;
	HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, pid);

	if(!GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) return 0;
	return (uint64_t)pmc.WorkingSetSize;
#else
	return 0;
#endif
}
