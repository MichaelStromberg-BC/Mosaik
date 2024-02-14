%module mosaik
%include "std_string.i"
%include "std_vector.i"

%{
#include "MosaikAlignment.h"
%}

/* Apply all of the integer typemaps to uint64_t */
%apply int { uint64_t };

/* Let's just grab the original header file here */
%include "MosaikAlignment.h"

namespace std {
   %template(AlignmentVector) vector<Mosaik::Alignment>;
   %template(RefSeqVector) vector<Mosaik::ReferenceSequence>;
   %template(ReadGroupVector) vector<Mosaik::ReadGroup>;
}
