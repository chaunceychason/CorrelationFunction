CC=icc
CFLAGS=-Wall -O3 -vec_report3 -g -lm -openmp


DM_CF_Parallel : HW3ParallelCorrelationFuncDM3_1.c
	$(CC) -o $@ $< $(CFLAGS)

.PHONY: clean

clean :
	rm -f DM_CF_Parallel





