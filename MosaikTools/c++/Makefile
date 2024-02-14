# ==================================
# define our source and object files
# ==================================

SOURCES=MosaikConversionMain.cpp MosaikAlignment.cpp
OBJECTS=$(SOURCES:.cpp=.o)

SOURCES2=MosaikReaderMain.cpp MosaikAlignment.cpp
OBJECTS2=$(SOURCES2:.cpp=.o)

CSOURCES=fastlz.c
COBJECTS=$(CSOURCES:.c=.o)

# ================
# compiler options
# ================

CPPFLAGS=-Wall -O3
LDFLAGS=-Wl

all: MosaikConversion MosaikReaderTest
 
MosaikConversion: $(OBJECTS) $(COBJECTS)
	@echo "- linking MosaikConversion"
	@$(CXX) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LIBS)

MosaikReaderTest: $(OBJECTS2) $(COBJECTS)
	@echo "- linking MosaikReaderTest"
	@$(CXX) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean

clean:
	rm -f *.o MosaikConversion MosaikReaderTest *~
