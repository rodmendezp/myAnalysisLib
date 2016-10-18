SHELL = /bin/bash

.DELETE_ON_ERROR:

.PHONY: all clean


ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs)
ROOTINCDIR  := $(shell $(ROOTCONFIG) --incdir)

CXX       := g++
CXXFLAGS  += -std=c++11 -O2 -Wall -fPIC $(ROOTCFLAGS)
LD        = g++
LDFLAGS   = -O2 $(ROOTLDFLAGS)

INCLUDES  := -I/$(ROOTINCDIR) -I$(CLASTOOL)/include -I$(ANALYSER)/include
LIBS      := $(ROOTLIBS) -L$(CLASTOOL)/slib/Linux -lClasTool -L$(ANALYSER)/slib -lTIdentificator

##############################################################################
all: myROOTLib.so

myROOTLib.so: myROOTLib.cpp
        $(CXX) $(CXXFLAGS) $(INCLUDES) -shared $< -o $@

clean:
        @rm -rf myROOTLib.so
