# AVARDA

As described in Monaco et al. 2021

Here we introduce the AntiViral Antibody Response Deconvolution Algorithm (AVARDA), a systematic framework for probabilistic analysis of highly multiplexed antiviral antibody epitope reactivity data. AVARDA integrates three key modules to generate a conservative and probabilistic assessment of antibody responses. The first module uses sequence alignment to define each reactive peptide’s relationship to a comprehensive database of all human viral genomes, which have been translated in all six reading frames. This permits a conservative elaboration of all peptide reactivities that could be associated with each potential viral infection. The second module constructs a sequence homology-based network graph for each virus’s reactive peptides, so as to define the minimum number of independent specificities (response breadth) required to produce the graph. The third module iteratively assigns each peptide to its most likely associated viral infection(s), according to a null model that considers the overall representation of each virus in the VirScan library. We permit individual peptides to be associated with multiple distinct infections, provided there is sufficient evidence for each viral infection on its own (i.e. in the absence of the shared peptides). AVARDA also indicates uncertainty in peptide-virus assignments when there is insufficient evidence to discriminate between infections by related viruses, but when there is sufficient evidence to conclude that an infection has indeed occurred. Linking these modules, the final output of AVARDA provides adjusted p-values for multiple hypothesis testing for exposure to each virus, along with the associated breadths of the antibody responses and the relationships between the reactive peptides.

## AVARDA Script

### Instructions
To run avarda download AVARDA.R and bin2 folder and place in desired directory location.

If using names other than Larman Lab names (i.e. VirScan_Lar_XXXXX) replace the first column of avarda_names.csv with other peptide names.

Set Inputs to run specific values


### Inputs 
Input1: path to file of data to be analyzed, rows are peptides, columns are samples, values are generally binary values (indicating hits)
Input2: thesholdiag value for input 1; set to 1 if using binary matrix
Input3: path to "dictionary" indicating all inter-peptide alignments (doesn't change between runs)
Input4: path to precomp table listing each viruses abundance in VirScan (does not change if using Virscan)
Input5: path to precomp table listing each viruses relative representation compared with other viruses (doesn't change) 
Input6: path to precomp matrix of Virscan peptide - virus alignments
Input7: Output path
Input8: Name of study 
Input9: Table converting between dataset names and avarda names
