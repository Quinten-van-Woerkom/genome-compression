CC=gcc
CXX=g++-8
RM=rm -f
CPPFLAGS=-Wall -std=c++17 -Iinclude
LDFLAGS=-lstdc++fs
LDLIBS=

SRCS=compress.cpp src/dna.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: clean compress

compress: $(SRCS)
	$(CXX) -o $@ $(SRCS) $(LDLIBS) $(LDFLAGS) $(CPPFLAGS)

clean:
	$(RM) $(subst .cpp, ,$(SRCS))