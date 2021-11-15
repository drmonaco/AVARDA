# AVARDA

#AVARDA Script


### Instructions
To run avarda place directory in preferred location and set input arguments as absolute paths to various above listed files.

If using names other than Larman Lab names replace the first column of avarda_names.csv with other peptide names.


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
