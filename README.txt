       simNGS -- software for simulating next-generation sequencing data
       
** Introduction
        simNGS is software for simulating "perfect" observations from Illumina
sequencing machines using the statistical models behind the AYB base-calling
software. Observations are perfect in the sense that they are generated
according the model and do not incorporate any additional effects that may be
present in real data ("dust", bubbles, merged clusters, sequence-heterogeneous
clusters, etc).
	simNGS takes fasta format sequences and a file describing the covariance
of noise between bases and cycles observed in an actual run of the machine,
randomly generates noisy intensities representing the signals for the sequence 
at each cycle and calculates likelihoods for all possible base calls.


** Compilation
	simNGS is written for BSD-style Unices according to the C99 standard and
should be relatively portable. Make files for both generic Linux and Mac Os-X 
in the the src/ directory.
	The blas and lapack libraries for matrix operations and linear algebra
must be installed and in the standard library path.

Compiling:
	cd src
	make -f Makefile.linux
which should produce a binary "bin/simNGS".

** Usage
	Fasta format sequence are read from stdin and log-likelihoods for the
sequences, after the addition of noise, are written to stdout in a format 
described in "Output format" below. Messages, progress indicators and a 
summary of errors in the generated data are written to stderr.
	simNGS accepts several command-line arguments to change its behaviour,
brief descriptions of which are available by running 'simNGS --help'.
More detailed descriptions are available in doc/parameters.txt

	cat sequences.fa | simNGS runfile > sequences.like


** Example
	Produces likelihoods for sequences in test100.fa and outputs to 
test.like, treating the (single-ended) runfile "s_2_0005.runfile" as 
paired-end.

cat test100.fa | bin/simNGS -p data/s_2_0005.runfile  > test.like

Description of runfile:
> # 76 cycle PhiX data from Sanger, single-ended. Very good run
> Treating single-ended model as paired-end.
> Using seed 1264523667
> Finished generating       100 sequences
> Summary of errors, calling by maximum likelihood
> Cycle  Count  Phred   lower, upper   Count  Phred   lower, upper
>   1:       8  10.97 (  8.24, 13.86)      6  12.22 (  9.04, 15.56)
>   2:       4  13.98 ( 10.07, 18.05)      4  13.98 ( 10.07, 18.05)

The final few lines are a per-cycle description of the raw rate if the bases
are called by maximum likelihood. "Count" is the raw number of errors 
generated, so the proportion of bases at a cycle with an error is 
Count / Number of Sequences. Phred is the Phred-like "quality" score for the 
bases, and lower and upper are an 95% confidence for the Phred score generated
by transforming a Wilson score interval.

Output, trimmed and artificially split over several lines:
> 2	5	1165	560
> 	1.441031e+00 9.830551e-01 2.178242e+00 2.034257e+00
> 	4.003852e+00 4.413710e+00 2.539974e+00 3.105548e+00
> ... more intensities
> 2	5	1668	739
> ... more 

The lane was 2, the tile was 5. The coordinates of the first cluster are 
1165,560 (randomly generated as if from a GAII machine [1], uniformly over the
tile and not thinned to take into account merged clusters). The likelihoods are 
encoded as -log(likelihood), so smaller values are the more likely, and the 
difference between these values (a function of the likelihood-ratio statistic) 
is a measure of the relative confidence between two calls.


** Output format
	The output from simNGS is a simple format based on the "_int.txt" from
the Illumina platform. Informally, the file consists of several lines, one line
per sequence, with tab separated fields. The first four fields contain the tile,
lane number and x & y coordinates of the cluster. There are a number of
additional field, equal to the number of cycles, containing (minus) the 
log-likelihoods for the base-calls at cycle: four space-separated non-negative 
real numbers in C-style "%e" form.

More formally, the grammar in Extended Backus-Naur form is:
	FILE ::= LINE +
	LINE ::= LANE , "\t" , TILE , "\t" , X , "\t" , Y , LOGLIKES + , EOL
	TILE ::= digit +
	LANE ::= digit +
	X    ::= digit +
	Y    ::= digit +
	LOGLIKES ::= "\t" , LOGLIKE , " " , LOGLIKE , " " , LOGLIKE , " " , LOGLIKE
	LOGLIKE  ::= non-zero digit , "." , digit + , EXPO 
	EXPO ::= "e" , ( "+" | "-" ) , digit +

Paired-end runs are concatenated, so appear as a longer single-ended run with
twice as many cycles.


** Runfile format
	The runfile describes how noise and cluster intensites was distributed
in a real run of an Illumina machine, as estimated by AYB.
It consists of:
 1) One line of free-text describing the run (prefixed by a '#').
 2) One line containing the lane number, tile number, and the shape and scale
of the distribution of the brightness of clusters.
 3) One line describing the dimensions of the covariance matrix for the first
end of the run, the number of rows and columns are each equal to four times the
number of cycles.
 4) Several lines, one per row of the covariance matrix, containing the elements
of that row.
	If the run was paired-ended, (3) and (4) are repeated for a second
matrix, which should have the same dimensions as the first.

	The matrices in runfile are covariance matrices for the noise in the run
and so must be symmetric positive-definite.

[1] X \in {0, ... , 1794}
    Y \in {0, ... , 2047}
