CFLAGS=-std=c99 -g -O2
CPPFLAGS=-I../include
LDFLAGS=-lm -lgsl -lgslcblas

.PHONY: all clean help

all: cawof

cawof: cawof.o entropy_tools.o wf_prange.o wf_stern.o wf_dumer.o wf_mmt.o wf_bjmm.o wf_nn.o
	$(CC) $(CFLAGS) -o cawof cawof.o entropy_tools.o wf_prange.o wf_stern.o wf_dumer.o wf_mmt.o wf_bjmm.o wf_nn.o $(LDFLAGS)

cawofa.o: cawof.c cawof.h ../include/entropy_tools.h 
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o cawof.o cawof.c

entropy_tools.o: entropy_tools.c ../include/entropy_tools.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o entropy_tools.o entropy_tools.c

wf_prange.o: wf_prange.c ../include/wf_prange.h ../include/entropy_tools.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o wf_prange.o wf_prange.c 

wf_stern.o: wf_stern.c ../include/wf_stern.h ../include/entropy_tools.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o wf_stern.o wf_stern.c 

wf_dumer.o: wf_dumer.c ../include/wf_dumer.h ../include/entropy_tools.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o wf_dumer.o wf_dumer.c 

wf_mmt.o: wf_mmt.c ../include/wf_mmt.h ../include/entropy_tools.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o wf_mmt.o wf_mmt.c 

wf_bjmm.o: wf_bjmm.c ../include/wf_bjmm.h ../include/entropy_tools.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o wf_bjmm.o wf_bjmm.c 

wf_nn.o: wf_nn.c ../include/wf_nn.h ../include/entropy_tools.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o wf_nn.o wf_nn.c 

bound_function.o: bound_function.c ../include/bound_function.h ../include/entropy_tools.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o entropy_tools.o entropy_tools.c

clean: 	
	@rm -f *.o *~ cawof

help: 
	@echo -e "Usage"
	@echo -e "make [all]\t Compile cawofa"
	@echo -e "make CaWoFa\t Compile a executable file from cawofa.c"
	@echo -e "make clean\t Return to the initial state of directory before the compilation\t\t "
	@echo -e "make help\t Give the principal target of Makefile with a short description"
