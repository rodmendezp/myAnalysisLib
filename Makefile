SHELL = /bin/bash

.DELETE_ON_ERROR:

.PHONY: all clean


ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs)
ROOTINCDIR  := $(shell $(ROOTCONFIG) --incdir)

CXX       := g++
CXXFLAGS  += -O2 -Wall -fPIC $(ROOTCFLAGS)
LD        = g++
LDFLAGS   = -O2 $(ROOTLDFLAGS)

INCLUDES  := -I/$(ROOTINCDIR) -I$(CLASTOOL)/include -I$(ANALYSER)/include -I/.
LIBS      := $(ROOTLIBS) -L$(CLASTOOL)/slib/Linux -lClasTool -L$(ANALYSER)/slib -lTIdentificator

SOURCES := myROOTUtils.cpp FitPCorr.cpp pCorrStepan.cpp
OBJECTS := $(SOURCES:.cpp=.o)

##############################################################################
all: libmyROOTLib.so

myROOTUtils.o: myROOTUtils.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@	

FitPCorr.o: FitPCorr.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
		
pCorrStepan.o: pCorrStepan.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

libmyROOTLib.so: $(OBJECTS)
		$(CXX) $(CXXFLAGS) $(INCLUDES) -shared -o $@ $^

clean:
                @rm -rf pCorrStepan.o FitPCorr.o myROOTUtils.o libmyROOTLib.so
