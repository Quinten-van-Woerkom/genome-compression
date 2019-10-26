CC=gcc
CXX=g++-8
RM=rm -f
CPPFLAGS=-Wall -std=c++17 -Iinclude
LDFLAGS=-lstdc++fs
LDLIBS=

SRCS=compress.cpp src/dna.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: compress

compress: $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDLIBS) $(LDFLAGS) $(CPPFLAGS)

%: %.o
	$(CXX) -o $@ $< $(LDFLAGS)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CPPFLAGS)

%: %.cpp
	$(CXX) -o $@ $< $(LDFLAGS)

clean:
	$(RM) $(OBJS)