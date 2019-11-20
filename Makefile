CC=gcc
CXX=g++-8
RM=rm -f
CPPFLAGS=-Wall -std=c++17 -Iinclude
LDFLAGS=-lstdc++fs
LDLIBS=

MAIN=compress.cpp
TEST=tests/test.cpp
SRCS=src/dna.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: compress test

test: $(SRCS) $(TEST)
	$(CXX) -o $@ $(TEST) $(SRCS) $(LDLIBS) $(LDFLAGS) $(CPPFLAGS)

compress: $(SRCS) $(MAIN)
	$(CXX) -o $@ $(MAIN) $(SRCS) $(LDLIBS) $(LDFLAGS) $(CPPFLAGS)

clean:
	$(RM) $(subst .cpp, ,$(SRCS))
	$(RM) $(subst .cpp, ,$(MAIN))
	$(RM) $(subst .cpp, ,$(TEST))
	$(RM) $(subst .cpp,.o,$(SRCS))
	$(RM) $(subst .cpp,.o,$(MAIN))
	$(RM) $(subst .cpp,.o,$(TEST))