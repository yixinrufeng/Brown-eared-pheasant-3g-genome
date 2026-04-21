//Parameters for the coalescence simulation program : fastsimcoal.exe
1 samples to simulate :
//Population effective sizes (number of genes)
NMB
//Haploid samples sizes 
20
//Growth rates: negative growth implies population expansion
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2 historical event
T1	0	0	0	NBOT	0	0
T2	0	0 	0	NANC	0	0	absoluteResize
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 9.1e-9 OUTEXP
