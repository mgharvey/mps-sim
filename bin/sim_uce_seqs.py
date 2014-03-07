#!/usr/bin/env python

"""

Name: sim_uce_seqs.py 
Author: Michael G. Harvey
Date: 15 October 2013

See readme for details.

"""

import os
import sys
import random
import re
import readline
import argparse
import ConfigParser
from subprocess import check_call, call, PIPE
from Bio import AlignIO


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"config_file",
			type=str,
			help="""The configuration file"""
		)
	return parser.parse_args()


def ConfigSectionMap(section, Config):
	dict1 = {}
	options = Config.options(section)
	for option in options:
		try:
			dict1[option] = Config.get(section, option)
			if dict1[option] == -1:
				DebugPrint("skip: %s" % option)
		except:
			print("exception on %s!" % option)
			dict1[option] = None
	return dict1


def interpret_ms(ms, seq_length, outgroup):
	# Parse ms command
	ms_commands = ms.split(' ') 
	if ms_commands[1] != "1":
		print "\nWARNING: Number of loci in ms command not equal to 1. Setting to 1."
		ms_commands[1] = "1"
	pop_sizes = list()
	if "-t" not in ms_commands:
		print "\nWARNING: Theta not specified in ms command. Scaling branch lengths to 1."		
		h = 1
	if "-I" not in ms_commands:
		print "\nWARNING: No '-I' command in ms string. Treating as single population."
		p = 1
	for i, ms_command in enumerate(ms_commands):
		if ms_command == "-t":
			th = float(ms_commands[i+1]) # theta
			h = float(th)/float(seq_length) # theta per site
		if ms_command == "-I":
			p = int(ms_commands[i+1])
			ntax = 0 
			q = 0
			for z in range(p):
				pop_sizes.append(int(ms_commands[i+(2+q)])) # pop sizes
				ntax += pop_sizes[q] # total number of taxa/samples
				q += 1
	if outgroup == "True":
		if pop_sizes[-1] != 2:
			print "\nWARNING: When using outgroup, last population size should be \"2\""			
	new_command_list = list()
	for i in range(len(ms_commands)):
		new_command_list.append(ms_commands[i])
		new_command_list.append(" ")
	new_command_string = ''.join(new_command_list)
	return new_command_string, pop_sizes, h, ntax, p


def create_output(out_dir):
	# Create directories to which output files will be written
	outdir = ("{0}".format(out_dir)) # output directory
	FNULL = open(os.devnull, 'w')
	call("mkdir -p {0}".format(outdir), shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
	os.chdir(outdir)
	call("mkdir {0}tmp".format(outdir), shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
	call("mkdir {0}nex".format(outdir), shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
	FNULL.close()


def add_zeros(out_dir):
	# Add zeros to beginning of single-digit samples in seq-gen output	
	os.chdir("{0}".format(out_dir))
	f1 = open("{0}tmp/seq-gen_output".format(out_dir), 'r')
	f2 = open("{0}tmp/seq-gen_output_reformatted".format(out_dir), 'wb')
	for line in f1:
		parts = line.split()
		if len(parts[0]) == 1:
			new = "000{0}".format(parts[0])
		elif len(parts[0]) == 2:
			new = "00{0}".format(parts[0])
		elif len(parts[0]) == 3:
			new = "0{0}".format(parts[0])
		elif len(parts[0]) == 4:
			new = "{0}".format(parts[0])
		new_line = "{0} {1}\n".format(new, parts[1])	
		f2.write(new_line)
	f1.close()
	f2.close()


def reorder_seqgen(outdir):
	# Re-order seq-gen output according to sample number
	f1 = open("{0}tmp/seq-gen_output_reformatted".format(outdir), 'r')
	f2 = open("{0}tmp/seq-gen_output_reordered".format(outdir), 'wb')
	next(f1)
	lines = [line for line in f1 if line.strip()]
	lines.sort()
	f2.writelines(lines)
	f1.close()
	f2.close()


def write_nexus(outdir, l, ntax, seq_length, p, pop_sizes):				
	# Export Nexus alignments of all loci	
	nf = open("{0}/nex/locus{1}.nex".format(outdir, l+1), 'wb')				
	nf.write("#NEXUS\n")
	nf.write("begin taxa;\n")
	nf.write("\tdimensions ntax={0};\n".format(ntax/2))
	nf.write("\ttaxlabels\n")
	for i in range (0, ((ntax)/2)): # write taxa
		nf.write("\tInd_{0}\n".format(i+1))
	nf.write(";\n")
	nf.write("end;\n")
	nf.write("\n")
	nf.write("begin characters;\n")
	nf.write("\tdimensions nchar={0};\n".format(seq_length))
	nf.write("\tformat datatype=dna missing=? gap=-;\n")
	nf.write("\tmatrix\n")
	i = 0
	j = 1	
	ntaxdone = 0
	while j < (p+1): # for each population
		count = range(ntaxdone+1, ntaxdone+pop_sizes[j-1]+1) # Generate list of random integers	
		ntaxdone += pop_sizes[j-1]
		rrs = random.sample(count, pop_sizes[j-1]) 
		k = 0
		while k < pop_sizes[j-1]: # for each taxon/sample in this population
			f1 = open("{0}/tmp/seq-gen_output_reordered".format(outdir, j), 'r')
			lines = f1.readlines()
			rr1 = rrs[k] # use random numbers to select first haplotype
			l1 = (lines[rr1-1])
			rr2 = rrs[k+1] # use random numbers to select second haplotype
			l2 = (lines[rr2-1]) 
			i += 1
			k += 2							
			o = 5
			nf.write("\tInd_{0}\t".format(i)) # write taxon name and sequence
			# Start loop to unphase haplotypes
			while o < (int(seq_length)+5): # for each base
				if l1[o] == "A":
					if l2[o] == "A":
						nf.write("A")
					elif l2[o] == "G":
						nf.write("R")
					elif l2[o] == "C":
						nf.write("M")
					elif l2[o] == "T":
						nf.write("W")
				elif l1[o] == "G":
					if l2[o] == "A":
						nf.write("R")
					elif l2[o] == "G":
						nf.write("G")
					elif l2[o] == "C":
						nf.write("S")
					elif l2[o] == "T":
						nf.write("K")
				elif l1[o] == "C":
					if l2[o] == "A":
						nf.write("M")	
					elif l2[o] == "G":
						nf.write("S")
					elif l2[o] == "C":
						nf.write("C")
					elif l2[o] == "T":
						nf.write("Y")
				elif l1[o] == "T":
					if l2[o] == "A":
						nf.write("W")
					elif l2[o] == "G":
						nf.write("K")
					elif l2[o] == "C":
						nf.write("Y")
					elif l2[o] == "T":
						nf.write("T")
				else:
					nf.write("N")
				f1.close()			
				o += 1
			nf.write("\n")	
		j += 1					
	nf.write(";\n")
	nf.write("end;\n")
	nf.close()


def write_gphocs(dataset_dir):
	dir = "{0}nex/".format(dataset_dir)
	files = list()
	prefiles = os.listdir(dir)
	for prefile in prefiles: # Remove hidden files
		if not prefile.startswith('.'):
			files.append(prefile)
	gf = open("{0}gphocs_input.txt".format(dataset_dir), 'wb')
	gf.write("{0}\n".format(len(files)))
	for i, file in enumerate(files):
		alignment = AlignIO.read("{0}{1}".format(dir, file), "nexus")
		gf.write("locus{0} {1} {2}\n".format(i+1, len(alignment), alignment.get_alignment_length()))
		for record in alignment:
			gf.write("{0} {1}\n".format(record.id, record.seq))	
		gf.flush()
	gf.close()


def write_nexus_concat(samples, alignment, dataset_dir, output):
	ff = open("{0}{1}".format(dataset_dir, output), 'wb')
	ff.write("#NEXUS\n\n")
	ff.write("Begin data;\n")
	ff.write("\tDimensions ntax={0} nchar={1};\n".format(len(samples), len(alignment)))
	ff.write("\tFormat datatype=dna gap=-;\n")
	ff.write("\tMatrix\n")
	for i, sample in enumerate(samples):
		ff.write("{0}\t".format(sample))
		for k in range(len(alignment)):
			ff.write(alignment[k][i])
		ff.write("\n")
	ff.write(";\n")
	ff.write("End;\n")
	ff.close()			


def extract_SNPs(dataset_dir, outgroup):			
	cat_align = list()
	cat_all_align = list()
	dir = "{0}nex/".format(dataset_dir)
	files = os.listdir(dir)
	all_samples = list()
	for m, file in enumerate(files): # Loop to get all sample names
		alignment = AlignIO.read("{0}nex/{1}".format(dataset_dir, file), "nexus")
		samples = list()
		for record in alignment:
			samples.append(record.id)
			if record.id not in all_samples:
				all_samples.append(record.id)
	poly = 0
	spoly = 0
	for m, file in enumerate(files): # Loop to extract polymorphic sites		
		alignment = AlignIO.read("{0}nex/{1}".format(dataset_dir, file), "nexus")
		alignment_polymorphic = list()	
		for w in xrange(alignment.get_alignment_length()):
			bases = list(alignment[:,w])
			if outgroup == "True":
				obases = bases[:-1]
			else:
				obases = bases
			uniq = set()
			for i, obase in enumerate(obases, 1):
				if obase not in ["-", "?", "N", "n"]:
					uniq.add(obase)			
			if len(uniq) > 1: # Is it polymorphic? 
				alignment_polymorphic.append(bases)	
				poly += 1
		for n in range(len(alignment_polymorphic)):
			cat_all_align.append(alignment_polymorphic[n])
		if len(alignment_polymorphic) > 0:
			# Randomly extract one polymorphic site from each alignment					
			alignment_single = list()
			rn = random.randint(0, len(alignment_polymorphic)-1)
			snp_seq = alignment_polymorphic[int(rn)]
			cat_align.append(snp_seq)
			spoly += 1	
	return all_samples, cat_align, cat_all_align, poly, spoly


def write_summary(dataset_dir, ms, genes, seq_length, outgroup, poly, spoly, biall):
	sf = open("{0}summary.txt".format(dataset_dir), 'wb')
	sf.write("ms command: {0}\n".format(ms))
	sf.write("number of loci: {0}\n".format(genes))
	sf.write("locus length: {0}\n".format(seq_length))
	sf.write("outgroup: {0}\n".format(outgroup))
	if outgroup == "True":
		sf.write("number of polymorphic sites (not including outgroup): {0}\n".format(poly))
		sf.write("number of loci with polymorphic sites (not including outgroup): {0}\n".format(spoly))
	elif outgroup == "False":
		sf.write("number of polymorphic sites: {0}\n".format(poly))
		sf.write("number of loci with polymorphic sites: {0}\n".format(spoly))
	sf.write("number of biallelic sites included in hapmap/dadi files: {0}".format(biall))
	sf.close()


def write_hapmap(dataset_dir, samples, alignment):
	biall = 0
	hf = open("{0}hapmap.hmp.txt".format(dataset_dir), 'wb')
	hf.write("rs\talleles\tchrom\tpos\tstrand\tassembly\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode")
	for sample in samples:
		hf.write("\t{0}".format(sample))
	hf.write("\n")
	for i in range(len(alignment)):
		bases = alignment[i]
		uniq = set()
		new_bases = list()
		for base in bases:	
			if base in ["A", "a"]:
				new_bases.append("A")
				new_bases.append("A")
			elif base in ["C", "c"]:
				new_bases.append("C")
				new_bases.append("C")
			elif base in ["G", "g"]:
				new_bases.append("G")
				new_bases.append("G")
			elif base in ["T", "t"]:
				new_bases.append("T")
				new_bases.append("T")
			elif base in ["M", "m"]:
				new_bases.append("A")
				new_bases.append("C")
			elif base in ["R", "r"]:
				new_bases.append("A")
				new_bases.append("G")
			elif base in ["W", "w"]:
				new_bases.append("A")
				new_bases.append("T")
			elif base in ["S", "s"]:
				new_bases.append("C")
				new_bases.append("G")
			elif base in ["Y", "y"]:
				new_bases.append("C")
				new_bases.append("T")
			elif base in ["K", "k"]:
				new_bases.append("G")
				new_bases.append("T")
			elif base in ["N", "n", "-", "?"]:
				new_bases.append("N")
				new_bases.append("N")
		uniqs = list(set(new_bases).intersection('ACGT'))
		if len(uniqs) == 2: # only include biallelic SNPs
			biall += 1
			hf.write("loc{0}\t{1}/{2}\t0\t{3}\t+\tNA\tSIM\tSIM\tNoRef\tCustom\tQC+".format(i+1, uniqs[0], uniqs[1], i+1))
			for j, sample in enumerate(samples):
				hf.write("\t{0}".format(alignment[i][j]))
			hf.write("\n")
	return biall


def make_pop_dict(p, pop_sizes, samples, outgroup):
	populations = list()
	projection = list()
	new_samples = list()
	for sample in samples:
		new_samples.append("{0}a".format(sample))
		new_samples.append("{0}b".format(sample))
	for i in range(int(p)):		
		if outgroup == "True":
			projection.append(min(pop_sizes[:-1]))
		elif outgroup == "False":
			projection.append(min(pop_sizes))		
		for j in range(pop_sizes[i]):
			populations.append("p{0}".format(i+1))
	if outgroup == "True":
		populations[-2] = 'outgroup'
		populations = populations[:-1]
		new_samples = new_samples[:-1]
		projection = projection[:-1]
	pop_assignments = dict(zip(new_samples, populations))
	return pop_assignments, projection


def data_dict_from_nex(nexus, pop_assignments, outgroup):
	alignment = AlignIO.read("{0}".format(nexus), "nexus")
	ntax = len(alignment)
	nchar = alignment.get_alignment_length()
	taxa = list()
	for record in alignment:
		taxa.append("{0}a".format(record.id))
		taxa.append("{0}b".format(record.id))
	if outgroup == "True":
		taxa = taxa[:-1]
	populations = set(pop_assignments.values())
	populations.discard('outgroup')
	for taxon,pop in pop_assignments.items():
		if pop=='outgroup':
			outgroup_ii = taxa.index(taxon)
			break
	else:
		outgroup_ii = None	
	dd = {}
	for i in range(nchar):
		site_name = 'site_{0}'.format(i+1)		
		bases = list(alignment[:,i])	
		new_bases = list()
		for base in bases:	
			if base in ["A", "a"]:
				new_bases.append("A")
				new_bases.append("A")
			elif base in ["C", "c"]:
				new_bases.append("C")
				new_bases.append("C")
			elif base in ["G", "g"]:
				new_bases.append("G")
				new_bases.append("G")
			elif base in ["T", "t"]:
				new_bases.append("T")
				new_bases.append("T")
			elif base in ["M", "m"]:
				new_bases.append("A")
				new_bases.append("C")
			elif base in ["R", "r"]:
				new_bases.append("A")
				new_bases.append("G")
			elif base in ["W", "w"]:
				new_bases.append("A")
				new_bases.append("T")
			elif base in ["S", "s"]:
				new_bases.append("C")
				new_bases.append("G")
			elif base in ["Y", "y"]:
				new_bases.append("C")
				new_bases.append("T")
			elif base in ["K", "k"]:
				new_bases.append("G")
				new_bases.append("T")
			elif base in ["N", "n", "-", "?"]:
				new_bases.append("N")
				new_bases.append("N")
		uniqs = set(new_bases).intersection('ACGT')
		if len(uniqs) == 2: # only include biallelic SNPs
			entry = {}
			dd[site_name] = entry
			segregating = tuple(uniqs)
			entry['segregating'] = segregating
			if outgroup_ii is not None:
				entry['outgroup_allele'] = new_bases[outgroup_ii]
			call_dict = {}
			entry['calls'] = call_dict
			for pop in populations:
				call_dict[pop] = [0,0]
			for taxon, new_base in zip(taxa, new_bases):
				pop = pop_assignments[taxon]
				# Ignore the outgroup here
				if pop is 'outgroup':
					continue
				if new_base in segregating:
					if new_base == segregating[0]:
						aa = 0
					elif new_base == segregating[1]:
						aa = 1
					call_dict[pop][aa] += 1
	return dd


def write_dadi(dd, dataset_dir, outgroup, pop_sizes, projection):
	import dadi
	# Generate frequency spectrum from data dictionary
	if outgroup == "False":
		if len(pop_sizes) == 1:
			fs = dadi.Spectrum.from_data_dict(dd, ['p1'], projection, polarized=False)
		elif len(pop_sizes) == 2:
			fs = dadi.Spectrum.from_data_dict(dd, ['p1', 'p2'], projection, polarized=False)
		elif len(pop_sizes) == 3:
			fs = dadi.Spectrum.from_data_dict(dd, ['p1', 'p2', 'p3'], projection, polarized=False)
		else:
			print "WARNING: dadi only accepts frequency spectra of 1 to 3 dimensions. You have {0} populations.".format(len(pop_sizes))
	elif outgroup == "True":
		if len(pop_sizes) == 2:
			fs = dadi.Spectrum.from_data_dict(dd, ['p1'], projection, polarized=True)
		elif len(pop_sizes) == 3:
			fs = dadi.Spectrum.from_data_dict(dd, ['p1', 'p2'], projection, polarized=True)
		elif len(pop_sizes) == 4:
			fs = dadi.Spectrum.from_data_dict(dd, ['p1', 'p2', 'p3'], projection, polarized=True)
		else:
			print "WARNING: dadi only accepts frequency spectra of 1 to 3 dimensions. You have {0} populations.".format(len(pop_sizes))
	# Print fs to file
	fs.tofile('dadi_input.fs')


def main():
	print "\n--------------------------------------------------------------------------------"
	# parse configuration file
	args = get_args()
	Config = ConfigParser.ConfigParser()
	Config.read("{0}".format(args.config_file))
	datasets = int(ConfigSectionMap("Arguments", Config)['datasets'])
	print "\nSimulating {0} dataset(s)\n".format(datasets)
	genes = int(ConfigSectionMap("Arguments", Config)['genes'])
	print "Simulating {0} genes/dataset\n".format(genes)
	seq_length = int(ConfigSectionMap("Arguments", Config)['seq_length'])
	print "Simulating alignments of {0} bp for each locus\n".format(seq_length)
	th = float(ConfigSectionMap("Arguments", Config)['theta_seqgen'])
	print "Per-site theta for seq-gen set at {0}\n".format(th)
	outgroup = str(Config.getboolean("Arguments", "outgroup"))	
	print "Outgroup: {0}\n".format(outgroup)	
	ms = str(ConfigSectionMap("Arguments", Config)['ms'])
	print "ms command: {0}\n".format(ms)
	out_dir = str(ConfigSectionMap("Output", Config)['out_dir'])	
	print "Output directory: {0}\n".format(out_dir)	
	print "Files to output:\n"
	summary = str(Config.getboolean("Output", "summary"))
	if summary == "True":
		print "\tSummary file"
	nexus = str(Config.getboolean("Output", "nexus"))
	if nexus == "True":
		print "\tNexus alignment for each locus"
	nexus_all = str(Config.getboolean("Output", "nexus_all"))
	if nexus_all == "True":
		print "\tNexus alignment of all SNPs concatenated"
	nexus_one = str(Config.getboolean("Output", "nexus_one"))
	if nexus_one == "True":
		print "\tNexus alignment of one random SNP per locus concatenated"
	hapmap = str(Config.getboolean("Output", "hapmap"))
	if hapmap == "True":
		print "\tHapMap genotype file containing one random SNP per locus"
	dadi = str(Config.getboolean("Output", "dadi"))
	if dadi == "True":
		print "\tFrequency spectrum for dadi"
	gphocs = str(Config.getboolean("Output", "gphocs"))
	if gphocs == "True":
		print "\tInput file for G-PhoCS\n"
	print "--------------------------------------------------------------------------------"
	# extract information from ms command
	ms_info = interpret_ms(ms, seq_length, outgroup)
	ms = ms_info[0]
	pop_sizes = ms_info[1]
	h = ms_info[2]
	#th = h # to set theta to ms theta/locus length
	ntax = ms_info[3]
	p = ms_info[4]
	# loop for datasets
	for i in range(datasets): 
		print "\nSimulating dataset {0} of {1} and generating Nexus alignments...\n \
			".format(i+1, datasets)
		dataset_dir = "{0}dataset_{1}/".format(out_dir, i+1)
		create_output(dataset_dir)
		# loop for genes
		for j in range(genes): 
			FNULL = open(os.devnull, 'w')
			# ms command
			check_call("ms {0} | tail +4 | grep -v // >{1}/tmp/ms_treefile ".format(ms, \
				dataset_dir), shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE) 
			# seq-gen command
			check_call("seq-gen -mHKY -l {0} -s {1} -i 0.5 -a 0.5 <{2}/tmp/ms_treefile \
				>{3}/tmp/seq-gen_output".format(seq_length, th, dataset_dir, dataset_dir), \
				shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE) 
			# reformat seq-gen output
			add_zeros(dataset_dir)
			reorder_seqgen(dataset_dir)
			# write nexus alignments
			if nexus == "True": 
				write_nexus(dataset_dir, j, ntax, seq_length, p, pop_sizes)
			sys.stdout.write("\rUCE loci simulated: {0}".format((int(j+1))))
			sys.stdout.flush()		
			FNULL.close()
		if gphocs == "True": 
			print "\n\nWriting G-PhoCS input file..."
			write_gphocs(dataset_dir)
		if nexus or summary or hapmap or dadi == "True":
			es = extract_SNPs(dataset_dir, outgroup)
			samples = es[0]
			cat_align = es[1]
			cat_all_align = es[2]
			poly = es[3]
			spoly = es[4]
			if nexus_all == "True": 
				print "\nWriting concatenated Nexus alignment of all SNPs..."
				write_nexus_concat(samples, cat_all_align, dataset_dir, "all_SNPs_concatenated.nex")
			if nexus_one == "True": 
				print "\nWriting concatenated Nexus alignment of one SNP per locus..."
				write_nexus_concat(samples, cat_align, dataset_dir, "one_SNP_per_locus_concatenated.nex")				
			if hapmap == "True": 
				print "\nWriting HapMap file..."
				biall = write_hapmap(dataset_dir, samples, cat_align) # cf. cat_all_align
				# default is to write hapmap with only single SNP per locus
			if summary == "True": 
				print "\nWriting summary statistics..."
				write_summary(dataset_dir, ms, genes, seq_length, outgroup, poly, spoly, biall)
			if dadi == "True": 
				print "\nWriting dadi allele frequency spectrum..."
				# Generate data dictionary from NeXus alignment
				mpd = make_pop_dict(p, pop_sizes, samples, outgroup)
				pop_assignments = mpd[0]
				projection = mpd[1]
				nexus = "{0}one_SNP_per_locus_concatenated.nex".format(dataset_dir)	
				# default is to write fs with only single SNP per locus
				dd = data_dict_from_nex(nexus, pop_assignments, outgroup)
				write_dadi(dd, dataset_dir, outgroup, pop_sizes, projection)    
		FNULL = open(os.devnull, 'w')
		call("rm -rf {0}tmp".format(dataset_dir), shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
		call("rm {0}seedms".format(dataset_dir), shell=True, stdout=PIPE, stderr=PIPE, stdin=PIPE)
		FNULL.close()
	print "\nSimulations complete"		
	print "\n--------------------------------------------------------------------------------\n"

if __name__ == '__main__':
    main()