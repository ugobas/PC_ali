# PC_ali
Multiple sequence and structure aligner PC_ali
Program PC_ali
Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa
Email: <ubastolla@cbm.csic.es>

(It includes a modified version of the needlemanwunsch aligner programmed by Dr. Andrew C. R. Martin in the Profit suite of programs, (c) SciTech Software 1993-2007)

PC_ali performs hybrid multiple structure and sequence alignments based on the structure+sequence similarity score PC_sim. Besides the MSA, it outputs pairwise simila
rity scores and divergence scores and neighbor-joining phylogenetic tree obtained with the hybrid evolutionary divergence measure based on PC_sim, and it outputs a pd
b file with multiply superimposed structures. Optionally, it computes violations of the molecular clock for each pair of proteins.

PC_ali takes as input either not aligned sequences (option -seq) or an MSA (option -ali). PDB file names must be specified as sequence names.

Usage:
PC_ali   -seq <sequences in FASTA format, with names of PDB files>
	 -ali <MSA file in FASTA format, with names of PDB files>
	# The pdb code is optionally followed by the chain index
	# Ex: >1opd.pdb A or >1opdA or >1opd_A

	 -pbdir <directory of pdb files>  (default: current directory)
	 -pdbext <extension of pdb files>  (default: none)
#### Optional parameters:
	 -out <Name of output files> (default: alignment file)
	 -ali_tm     ! Perform pairwise alignments that target TM score
	 -ali_co     ! Perform pairwise alignments that target Contact Overlap
	 -ali_ss     ! Perform alignments that target sec.structure
	 -ss_mult    ! alignments that target sec.structure are MSA
	 -shift_max <Maximum shift for targeting sec.str.>
	 -print_pdb    ! Print structure superimposition in PDB format
	 -print_sim    ! Print similarity measures
	 -print_div    ! Print divergence measures
	 -print_cv     ! Print clock violations
	 -func <file with function similarity for pairs of proteins>

Computed similarity measures:
(1) Aligned fraction ali,
(2) Sequence identity SI,
(3) Contact overlap CO,
(4) TM-score TM (Zhang & Skolnick Proteins 2004 57:702)
(5) PC_sim, based on the main Principal Component of the four above similarity scores
They are printed in <>.prot.sim for all pairs of protein sequences, and also for multiple conformations of the same sequence (if present) if required with -print_sim

Computed divergence measures:
(1) Tajima-Nei divergence TN=-log((SI-S0)/(1-S0)) with S0=0.06 (Tajima F & Nei 1984, Mol Biol Evol 1:269),
(2) Contact_divergence CD=-log((q-q0(L))/(1-q0(L)) (Pascual-Garcia et al Proteins 2010 78:181-96),
(3) TM_divergence=-log((TM-TM0)/(1-TM0)), TM0=0.167.
(4) PC_divergence=-log((PC-PC0)/(1-PC0)), PC0 linear combination of S0, TM0, CO(L) and nali0=0.5.
They are printed in <>.prot.div for all pairs of protein sequences, and also for multiple conformations of the same sequence (if present) if required with -print_div

Flux of the program:
(1) In the modality -ali, the program starts from the pairwise alignments obtained from the input MSA. In the modality -seq the starting pairwise alignments are built
 internally.
(2) The program then modifies the pairwise alignments by targeting PC_sim. The similarity matrix is constructed recursively, using the input pairwise alignment for co
mputing the shared contacts and the distance after optimal superimposition (maximizing the TM score) for all pairs of residues and obtaining a new alignments. Two ite
rations are usually enough for getting good results. Optionally, for the sake of comparison, the program can target the TM score (-ali_tm), the Contact Overlap (-ali_
co) and the secondary structure superposition (-ali_ss).
(3) Then, the program builds the multiple alignment based on the maximal cliques of the pairwise alignments. This computation does not require neither a guide tree no
r gap penalty parameters and in most cases it is faster than the progressive multiple alignment.
(4) Finally, the program runs iteratively progressive multiple alignments using as guide tree the average linkage tree obtained with the PC_Div divergence measure of 
the previous step and using as starting alignment the previous multiple alignment. The best MSA is selected as the one with the maximum value of the average PC simila
rity score.
(5) The program prints the optimal MSA and the Neighbor Joining tree obtained from the corresponding PC_Div divergence measure.
(6) Optionally, if -print_pdb is set, the program prints the multiple superimposition obtained by maximizing the TM score.
(7) Furthermore, if -print_cv is set, the program computes and prints for all four divergence measures the violations of the molecular clock averaged over all possibl
e outgroups identified with the Neighbor-Joining criterion, and the corresponding significance score.

In the first pairwise phase, the program computes similarity and divergence scores for all pairs of protein structures. It then clusters all conformations of the same
 protein and computes the structural similarity (divergence) between two proteins as the maximum (minimum) across all the examined conformations.
The similarity and divergence scores are computed for the starting alignment, for the modified pairwise alignments that target different similarity scores (TM score, 
contact overlap and PC_sim) and for the best multiple alignment.

COMPILE:
>unzip PC_ali.zip
>make
>cp PC_ali ~/bin/ (or whatever path directory you like)

RUN:
>PC_ali -seq <sequence file> -pdbdir <path to PDB files>

EXAMPLE: PC_ali -seq 50044_Mammoth.aln -pdbdir <PDBPATH>
(all PDB files named in 50044_Mammoth.aln must be in <PDBPATH>)

OUTPUT:
MSA (PCAli.fas),
NJ tree (PCAli.tree), 
structure similarity (.sim) and structure divergence scores (.div) for each protein pair,
correlations between different types of sequence and structure identity (.id),
MSA of secondary structure (_ss.msa)
