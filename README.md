# find_VSG_Ndomains
# Identify the N-terminal Domain Type of a VSG Sequence

The N-terminal domain of a VSG sequence can be classified as one of the following types/subtypes:

Type A:
* A1
* A2
* A3

Type B:
* B1
* B2

## Usage

```
python find_VSG_Nterm.py file path
```

For more details about arguments see [Input](#input).

## Required Software
* [Python](https://www.python.org/) (written using Python 2.7.15)
* [Biopython](https://anaconda.org/anaconda/biopython) (written using version 1.68)
* [pandas](https://anaconda.org/anaconda/pandas) (written using version 0.24.2)
* [HMMER](http://hmmer.org) (written using version 3.1b2)
	- HMMER must be installed in PATH.

## Input
* (1) FASTA file containing one or more VSG protein sequences. - (`file`)
* (1) Path to directory containing HMM profiles to run sequences against. - (`path`)

### Required files for HMMER hmmscan
* VSG N-terminal first pass TypeA and TypeB hmm profile ([Wickstead et al.](https://www.sciencedirect.com/science/article/pii/S0166685114000772)) - (`VSG-N-mergeAB.hmm`)
* VSG N-terminal second pass sequence recovery profile - (`VSG-N-DomainBoundary.hmm`)
* VSG N-terminal TypeSubtype hmm profile - (`VSG-N-TypeSubtype.hmm`)

HMM profile files are located in [HMMprofiles](HMMprofiles).  
The files must be in a single directory (and the path to this directory provided as `path`).  
HMM profile files must first be compressed and indexed using HMMER hmmpress.

For example:
```
hmmpress profile.hmm
```

## Output
* (1) Terminal text containing progress of the executed script.
* (1) CSV file containing a summary of the identified N-terminal domain(s).
* (1) FASTA Directory containing:
	- (1) FASTA file containing the most probable N-terminal domain(s) of the input VSG(s).
	- (1) intermediate FASTA file containing VSG sequences that failed first pass to determine N-domain boundary.
	- (1) FASTA file containing VSG sequences for which N-terminal domains could not be defined.
* (1) HMM Directory containing:
	- (2) HMMER hmmscan table output files containing data to determine the N-terminal domain boundary of each VSG sequence as well as type (TypeA and TypeB).
	- (2) HMMER hmmscan standard output files containing annotation and statistical information about each of the input VSG(s) (TypeA and TypeB).
	- (1) HMMER hmmscan table output file containing subtype data on the N-terminal domain sequence(s) of the input VSG(s) (TypeA1-3/TypeB1-2).
	- (1) HMMER hmmscan standard output file containing annotation and statistical information about the N-terminal domain(s) subtype of the input VSG(s) (TypeA1-3/TypeB1-2).
