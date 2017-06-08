CFLAGS=-std=c99 -Wall -Wextra -g -O2
CPPFLAGS=-I../include
LDFLAGS=-lm -lgsl -lgslcblas

.PHONY: all clean help

all: library bound behaviour

library:
	@cd ../src && $(MAKE)

bound: bound.o ../src/entropy_tools.o ../src/bound_function.o ../src/wf_prange.o ../src/wf_stern.o ../src/wf_dumer.o ../src/wf_mmt.o ../src/wf_bjmm.o ../src/wf_nn.o
	$(CC) $(CFLAGS) -o bound bound.o ../src/entropy_tools.o ../src/bound_function.o ../src/wf_prange.o ../src/wf_stern.o ../src/wf_dumer.o ../src/wf_mmt.o ../src/wf_bjmm.o ../src/wf_nn.o $(LDFLAGS)

bound.o: bound.c ../include/entropy_tools.h ../include/bound_function.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o bound.o bound.c

behaviour: behaviour.o ../src/entropy_tools.o ../src/entropy_tools.o ../src/wf_prange.o ../src/wf_stern.o ../src/wf_dumer.o ../src/wf_mmt.o ../src/wf_bjmm.o ../src/wf_nn.o
	$(CC) $(CFLAGS) -o behaviour behaviour.o ../src/entropy_tools.o ../src/wf_prange.o ../src/wf_stern.o ../src/wf_dumer.o ../src/wf_mmt.o ../src/wf_bjmm.o ../src/wf_nn.o $(LDFLAGS)

behaviour.o: behaviour.c ../include/entropy_tools.h ../include/bound_function.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o behaviour.o behaviour.c

clean: 	
	@rm -f *.o *~ bound /src/*.o *~
	@rm -f *.o *~ behaviour bound /src/*.o *~
	
help: 
	@echo -e "Usage"
	@echo -e "make [all]\t Compile behaviour and bound"
	@echo -e "make bound\t Compile a executable file from bound.c"
	@echo -e "make behaviour\t Compile a executable file from behaviour.c"
	@echo -e "make clean\t Return to the initial state of directory before the compilation\t\t "
	@echo -e "make help\t Give the principal target of Makefile with a short description"
