CC = gcc
CXX = g++

CFLAGS = -std=c11 -Wall
CXXFLAGS = -std=c++17 -Wall
LDLIBS = -lm

APPNAME = Ne2x
CPPAPP = Ne2x_cpp
LIBNAME = libne2x.a

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    CFLAGS += -static
endif

.PHONY: all clean test asan

all: $(APPNAME) $(LIBNAME) $(CPPAPP)

$(APPNAME): Ne2x.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

Ne2x.o: Ne2x.c Ne2x_api.h
	$(CC) $(CFLAGS) -c -o $@ Ne2x.c

Ne2x_lib.o: Ne2x.c Ne2x_api.h
	$(CC) $(CFLAGS) -DNE2X_NO_MAIN -c -o $@ Ne2x.c

$(LIBNAME): Ne2x_lib.o
	ar rcs $@ $^

main.o: main.cpp Ne2x_api.h
	$(CXX) $(CXXFLAGS) -c -o $@ main.cpp

$(CPPAPP): main.o $(LIBNAME)
	$(CXX) $(CXXFLAGS) -o $@ main.o -L. -lne2x $(LDLIBS)

test: $(APPNAME) $(CPPAPP)
	bash regression_test.sh

asan:
	$(CC) -fsanitize=address -g3 -std=c11 -Wall -o Ne2x_asan Ne2x.c -lm

clean:
	rm -f Ne2x.o Ne2x_lib.o main.o $(APPNAME) $(CPPAPP) $(LIBNAME) Ne2x_asan *.d
