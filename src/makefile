CPLEX_HOME=/opt/ibm/ILOG/CPLEX_Studio1262/cplex
INCLUDE=$(CPLEX_HOME)/include
STATIC_LIB=$(CPLEX_HOME)/lib/x86-64_linux/static_pic
DYNAMIC_LIB=$(CPLEX_HOME)/bin/x86-64_linux
CPLEX_LIB_STATIC=cplex
CPLEX_LIB_DYNAMIC=cplex1262
LFLAGS=-lrt -lm -lpthread
CFLAGS=-DIL_STD -I$(INCLUDE) -O0 -g3 -c -fmessage-length=0 -std=c++0x
CC=g++

NRP-static:    main_NRP.o
	$(CC) -L$(STATIC_LIB) -o NRP ./main_NRP.o -l$(CPLEX_LIB_STATIC) $(LFLAGS)
NRP-dynamic:    main_NRP.o
	$(CC) -L$(DYNAMIC_LIB) -o NRP ./main_NRP.o -l$(CPLEX_LIB_DYNAMIC) -ldl $(LFLAGS)
	@echo "Remember when running to do: export LD_LIBRARY_PATH=\$$LD_LIBRARY_PATH:$(DYNAMIC_LIB)";
main_NRP.o:     main_NRP.cpp
	$(CC) $(CFLAGS) main_NRP.cpp

clean:
	rm NRP main_NRP.o main_NRP.d

