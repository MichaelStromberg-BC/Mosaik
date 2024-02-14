# ================
# compiler options
# ================

# this is the directory that contains perl.h, XSUB.h, and Extern.h
PERL_CORE_DIR=/usr/lib64/perl5/5.10.0/x86_64-linux-thread-multi/CORE
SWIG=swig

all: SWIG
 
SWIG:
	$(SWIG) -c++ -perl MosaikAlignment.i
	$(CC) -fpic -c fastlz.c
	$(CXX) -fpic -c MosaikAlignment.cpp MosaikAlignment_wrap.cxx -I$(PERL_CORE_DIR)
	$(CXX) -shared MosaikAlignment.o MosaikAlignment_wrap.o fastlz.o -o mosaik.so

.PHONY: clean

clean:
	rm -f *.o *.so *~ *.pm *.cxx
