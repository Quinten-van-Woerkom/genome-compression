CC=gcc
CXX=clang++
RM=rm -f
CPPFLAGS=-Wall -std=c++17 -Iinclude
ADDED_CPPFLAGS=
LDFLAGS=-lstdc++fs
LDLIBS=

MAIN=compress.cpp
TEST=tests/test.cpp
JUMP=local_alignment.cpp
SRCS=src/dna.cpp src/file_reader.cpp src/shared_tree.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

release: ADDED_CPPFLAGS=-O3

all release: compress test local_alignment

test: $(SRCS) $(TEST)
	$(CXX) -o $@ $(TEST) $(SRCS) $(LDLIBS) $(LDFLAGS) $(CPPFLAGS) $(ADDED_CPPFLAGS)

compress: $(SRCS) $(MAIN)
	$(CXX) -o $@ $(MAIN) $(SRCS) $(LDLIBS) $(LDFLAGS) $(CPPFLAGS) $(ADDED_CPPFLAGS)

local_alignment: $(SRCS) $(JUMP)
	$(CXX) -o $@ $(JUMP) $(SRCS) $(LDLIBS) $(LDFLAGS) $(CPPFLAGS) $(ADDED_CPPFLAGS)

clean:
	$(RM) $(subst .cpp, ,$(SRCS))
	$(RM) $(subst .cpp, ,$(MAIN))
	$(RM) $(subst .cpp, ,$(JUMP))
	$(RM) test
	$(RM) $(subst .cpp,.o,$(SRCS))
	$(RM) $(subst .cpp,.o,$(MAIN))
	$(RM) $(subst .cpp,.o,$(TEST))
	$(RM) $(subst .cpp,.o,$(JUMP))