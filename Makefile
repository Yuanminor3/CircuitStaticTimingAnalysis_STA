CC=g++
CFLAGS=-Wall -std=c++11

all: sta

sta: main.cpp
	$(CC) $(CFLAGS) main.cpp -o sta

clean:
	rm -f sta
	rm -f ckt_traversal.txt

b15_1: sta
	./sta NLDM_lib_max2Inp b15_1.isc

b18_1: sta
	./sta NLDM_lib_max2Inp b18_1.isc

b19_1: sta
	./sta NLDM_lib_max2Inp b19_1.isc

b22: sta
	./sta NLDM_lib_max2Inp b22.isc

c17: sta
	./sta NLDM_lib_max2Inp c17.isc

c1908_: sta
	./sta NLDM_lib_max2Inp c1908_.isc

c2670: sta
	./sta NLDM_lib_max2Inp c2670.isc

c3540: sta
	./sta NLDM_lib_max2Inp c3540.isc

c5315: sta
	./sta NLDM_lib_max2Inp c5315.isc

c7552: sta
	./sta NLDM_lib_max2Inp c7552.isc


.PHONY: clean run

run: sta
	@echo "Usage: make run ARGS= "lib_file ckt_file"
	./sta $(ARGS)