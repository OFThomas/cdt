CC=gcc
CPP=g++
RM=rm -f
PROG:=numerics

CPPFLAGS=-g $(shell root-config --cflags)
LDFLAGS=-g $(shell root-config --ldflags)
LDLIBS=-lpython3.6m -lm -I/usr/include/python3.6
#$(shell root-config --libs)

SRCS=numerics.cpp
OBJS=$(subst .cc,.o,$(SRCS))

all: $(PROG)

tool: $(OBJS)
	$(CXX) $(LDFLAGS) -o $(PROG) $(OBJS) $(LDLIBS) 

clean:
	$(RM) $(OBJS)

