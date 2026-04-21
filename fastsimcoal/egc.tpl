//Parameters for the coalescence simulation program : fastsimcoal.exe
2 samples to simulate :
//Population effective sizes (number of genes)
NMR
NMB
//Haploid samples sizes 
20
20
//Growth rates: negative growth implies population expansion
0
R1
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0	0	  	
0  	0
//Migration matrix 1
0       MIG1
MIG2    0
//Migration matrix 2
0	0   
0       0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
2 historical event
CHANGM 1 1 0 1 0 1
TSP 0 1 1 NANC 0 2 absoluteResize
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 9.1e-9 OUTEXP
