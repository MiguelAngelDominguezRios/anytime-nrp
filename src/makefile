
CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio1262/cplex
INCLUDE=$(CPLEX_HOME)/include
LINK=$(CPLEX_HOME)/lib/x86-64_linux/static_pic
LFLAGS=-lcplex -lrt -lm -lpthread
CFLAGS=-DIL_STD -I$(INCLUDE) -O0 -g3 -c -fmessage-length=0 -std=c++0x
CC=g++

NRP:	main_NRP.o
	$(CC) -L$(LINK) -o NRP ./main_NRP.o $(LFLAGS)
main_NRP.o: 	main_NRP.cpp algorithms2.h biobjective_algorithms.h
	$(CC) $(CFLAGS) main_NRP.cpp

clean:
	rm NRP main_NRP.o

