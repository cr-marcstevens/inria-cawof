CFLAGS=-std=c99 -g -O2
CPPFLAGS=-I../include
LDFLAGS=

.PHONY: all clean help

all: build

build:
	@cd src && $(MAKE)
	@cd examples && $(MAKE)

clean:
	@cd src && $(MAKE) clean
	@cd examples && $(MAKE) clean

help:
	@echo -e "Usage"
	@echo -e "make [all]\t Build cawofa and examples"
	@echo -e "make build\t Build a executable file from src/cawofa.c"
	@echo -e "make examples\t Build executable files from examples/"
	@echo -e "make clean\t Return to the initial state of directory before the compilation"
	@echo -e "make help\t Give the principal target of Makefile with a short description"
