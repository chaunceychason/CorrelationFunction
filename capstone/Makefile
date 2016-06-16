CC=icc
CFLAGS=-Wall -O3 -vec_report3 -g -lm -openmp

all: DM_CF_Parallel DM_CF_Parallel_Variable DM_CF_Parallel_Research_Variable DM_CF_GRIDED_PARALLEL Serial_DM_CF

DM_CF_Parallel : ParallelCorrelationFuncDM3_2.c
	$(CC) -o $@ $< $(CFLAGS)


DM_CF_Parallel_Variable : ParallelCorrelationFuncDM3_3.c
	$(CC) -o $@ $< $(CFLAGS)

DM_CF_Parallel_Research_Variable : ParallelCorrelationFuncDMResearch.c
	$(CC) -o $@ $< $(CFLAGS)

DM_CF_GRIDED_PARALLEL : ParallelCorrelationFuncDMResearchGRIDED.c
	$(CC) -o $@ $< $(CFLAGS)

Serial_DM_CF : SerialCorrelationFunc_FLOPS.c
	$(CC) -o $@ $< $(CFLAGS)

.PHONY: clean

clean :
	rm -f DM_CF_Parallel
	rm -f DM_CF_Parallel_Variable
	rm -f DM_CF_Parallel_Research_Variable
	rm -f Serial_DM_CF
	rm -f DM_CF_GRIDED_PARALLEL
