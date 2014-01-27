INTRODUCTION
-------

mps-sim is a repository for code used to simulate massively parallel sequencing data. 

License
-------

The code within this repository is available under a 3-clause BSD license. See the License.txt file 
for more information.

Citation
--------

If you use this documentation or the mps-sim software for your own research, please cite:

Harvey, M. G., B. T. Smith, T. C. Glenn, B. C. Faircloth, and R. T. Brumfield. 2013. Sequence capture versus restriction site associated DNA sequencing for phylogeography. arXiv:1312.6439.

Please also provide the link to this software repository in your manuscript:

    https://github.com/mgharvey/mps-sim

Dependencies
--------

This list may be incomplete.

- `ms (Hudson 2002) <http://home.uchicago.edu/~rhudson1/source/mksamples.html>`_
- `seq-gen (Rambaut and Grass 1997) <http://tree.bio.ed.ac.uk/software/seqgen/>`_
- `BioPython (Cock et al. 2009) <http://biopython.org/wiki/Main_Page>`_

Make sure ms and seq-gen are in $PATH. 

Usage
--------

This repository contains two scripts: sim_uce_seqs.py simulates data resembling sequences from 
sequence capture of genomic ultraconserved elements (e.g., Faircloth et al. 2012), and 
sim_gbs_snps.py simulates data resembling genotypes from genotyping by sequencing (e.g., Hohenlohe 
et al. 2010, Elshire et al. 2011). Both take the same arguments. To run the scripts, at the command 
prompt type:

```
python sim_uce_seqs.py config_file.txt
```

or

```
python sim_gbs_snps.py config_file.txt
```

after replacing config_file.txt with the name of your configuration file. The configuration file
contains user-provided information about how to simulate data and which output files to generate. 
The following pieces of information should be provided (see Config_file_example.txt for format):

- number of datasets 

- number of loci per dataset

- length of loci 

- per-site theta value for seq-gen
-- This scales trees before simulating sequence data. See the seq-gen manual for details.

- is an outgroup being included (see "ms command string" below for details)? 

- maximum number of SNPs/locus (sim_gbs_snps.py only)
-- Most pipelines for GBS data set a maximum number of SNPs. Loci above this threshold are removed 
because they are more likely to contain paralogous sequences.

- ms command string 
-- The number of loci in the command must be 1. If you would like to include outgroup information (e.g. 
for rooting trees or polarizing SNPs in the AFS), make sure the LAST POPULATION in the ms command 
comprises two haplotypes and diverges first (ie. coalesces last). If you would like to generate 
unphased diploid input files for G-PhoCS, all populations must have an even number of samples equal 
to 2x the desired number of diploid individuals.

- the output directory

- if you would like to output the following files (default is "Yes" for all):
	- a file with summary statistics
	- separate Nexus alignments of sequences 
	- a Nexus alignment of all SNPs concatenated 
	- a Nexus alignment of one random SNP from each locus concatenated
	- a hapmap genotype file (includes 1 random SNP/locus)
	- input files for dadi (Gutenkunst et al. 2009; includes 1 random SNP/locus)
	- input files for G-PhoCS (Gronau et al. 2011)
