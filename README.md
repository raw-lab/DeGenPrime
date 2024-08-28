# DeGenPrime

[Introduction](#introduction) <br />
[Citation](#citation) <br />
[Installation](#installation)  <br />
[Options](#options)  <br />
[Input](#input)  <br />
[Output](#output)  <br />
[Examples](#examples)  <br />
[Copyright](#copyright)  <br />

### Introduction
DeGenPrime selects the top PCR primer pairs for one or more phylogenetically similar DNA sequences which are aligned or not aligned on the basis of minimizing melting temperature difference for forward and reverse primers which pass the following filter checks: <br />
- Low Degeneracy
- Few Deletions
- GC content within the 40-60% range.
- Low repetition
- Non-complementary ends
- Minimal risk of hairpins or self and cross-dimerization
- Melting temperature within specified range.
- The range of melting temperatures for PCR primers can be specified by the user but must be within the absolute range of 50.0 – 65.0 degrees Celsius.
DeGenPrime runs off hard filters with no exceptions.
If no primers are found that can pass all of these filters the program will warn the user that no suitable primers were found.

### Installation

#### Global install
```console
git clone https://github.com/raw-lab/DeGenPrime.git
cd DeGenPrime/
mkdir build
cd build
cmake ../src
make -j4
sudo make install
```
#### Local install
```console
git clone https://github.com/raw-lab/DeGenPrime.git
cd DeGenPrime/
mkdir -p ~/bin
cd bin
cmake ../src
make -j4
make install DESTDIR=~/bin
```
Note: Do not forget to export your path:
```console
export PATH=$PATH:/home/DeGenPrime/build/degenprime #for example home directory
```


DeGenPrime is designed to run as a stand-alone console application on any platform capable of running C++ applications.  The program does try to align a file that is misaligned by calling MAFFT within the program, so you must either have this installed or manually align your sequences to use DeGenPrime.

### Options  

```console
./degenprime [--tags] <filename>
```
The filename must always be the last argument, or this program will throw a segmentation error.

#### Command-line arguments
Valid tags include:
```console
--amplicon:           <int>, Set the minimum amplicon length. 
--begin:              <int>, Set the beginning nucleotide. If the user inputs an integer < 0, the program will default this value to zero.
--degenerate          this tag directs the program to use hard filters and some limited degeneracy instead of basing primer design off the consensus sequence and minimizing penalty.  This method runs significantly slower than the default approach.
--delta_g:            <int>, Sets the minimum allowed gibbs free energy for the repetition filter.  Default value is -4.0 kcal/mol.
--end:                <int>, Set the ending nucleotide. If the user inputs an integer > the number of base pairs in the entire sequence, then the program will default this to the last nucleotide in the sequence.
--global or --g,       for lists of sequences that are misaligned, this tag specifies that the file should run MAFFT for global alignment.
--help or --h,         prints this help menu.
--input_file:<file>   this tag species the input file for the program.  It must be an aligned file in fasta or clustal form.  This tag is required for program operation.
--local or --l,        for lists of sequences that are misaligned, this tag specifies that the file should run MAFFT for local alignment.
--max_primer_len:  <int>, Sets the maximum length of the desired primer.  This has a default value of 22 and cannot be larger than 25.
--min_primer_len:  <int>, Sets the minimum length of the desired primer.  This has a default value of 20 and cannot be less than 18.
--min_temp:           <int>, Sets the minimum primer melting temperature. This has a minimum value of 50.0 (degrees Celsius) and must be smaller than --max_temp.
--max_temp:           <int>, Sets the maximum primer melting temperature. This has a maximum value of 65.0 (degrees Celsius) and must be larger than --min_temp.
--primer_conc:        <int>, Sets the concentration of the PCR primer in nM. This has a minimum value of 50.0 nM, and this program will raise any value smaller to this value.
--protein            This tag will cause the program to interpret the input file as a sequence of amino acids and translate the amino acids into nucleotides then save the output file as <filename>_protein.faa.
--output_file:<file> This tag specifies the output csv file where the program data will be saved.  If this tag is not included, the filename will write output to the same filename as the input filename but will replace the file extension with '.csv'
--salt_conc:          <int>, Sets the concentration of monovalent ions in mM. This has a minimum value of 50.0 mM, and this program will raise any value smaller to this value.
--search_fwd:         <string>, The string represents a forward primer.  Searches the collected list of forward primers to see if the argument primer is within them.  If this primer is included, it gives its relative position on the ordered list by penalty.  This tag is not compatible with --test.
--search_rev:         <string>, This tag is similar to --search_fwd, but is for the reverse primer list.  It is also not compatible with --test.
--max_primers:        <int>, Sets the maximum number of output primers. This has a maximum value of 10 and this program will reduce any value larger to this value.
--test:               <string>, The string represents a single primer.  Runs the primer through all filters.  Returns the thermodynamic values of this primer as well as any filters this primer would not pass and its calculated penalty.  This tag is incompatible with --search tags.  Any primer smaller or larger than the size limits will show primer outside size range.
```

### Input
Input files to DeGenPrime can use a variety of formats including single sequences to alignment files.<br /> 
- Nucleotide fasta formats (.faa, .fasta, .ffn, .fna)<br />
- Nucleotide alignments (.clust)<br />

### Output
DeGenPrime will output a few progress messages as the program runs, the recommended primers and their details, and the program runtime at the end to the console.  DeGenPrime also records the operation details and the primers to an output file.  The output filename is selected based on the specified input file.<br />  

#### Standard outputs
- filename.dgp (a text file containing the list of primers, consensus sequence, and details) <br />
- Any errors that occur during runtime will also be written on the console.<br />

### Examples
Your find modifications on chromosome 21 in the p arms above the centromere.  You find out the length of this region from the tips of the 3’ end of the telomeres to the centromere is about 12 million base pairs.  The average length of telomeres on chromosome 21 for a human zygote is 10 thousand base pairs.  Your team collects data from the mutated zygotes’ chromosome 21 and other genetic information collected from NCBI on chromosome 21 into a fasta file called zygote_21.faa.  You want to find 10 good primer pairs on this chromosome to amplify in a PCR reaction.<br /> You would use:<br />

```console
./DeGenPrime --begin:10000 --end:12000000 --max_primers:10 zygote_21.faa
```

You are an immunologist who wants to identify mutations that might have occurred in the genes of a local strain of the influenza A virus so you can produce a new vaccine for the upcoming flu season.  Microbiologists have reported changes to the surface proteins of infected cells and this leads you to hypothesize that a mutation has occurred in segment 4 or segment 6 of the virus.  You obtain the genetic data of these segments from last season’s influenza A as well as data from the current strain into two files called influenzaA_4.faa and influenzaA_6.faa.  You know that the surface proteins normally contain 500 +/- 50 amino acids (which implies 1500 +/- 150 bps in the coding sequence) and you want the top 5 primers for a PCR reaction that will cover at least 40% of the coding region.  Your files are not aligned. <br />  
You would use:<br />

```console
./DeGenPrime --amplicon:660 --local influenzaA_4.faa or influenzaA_6.faa
```

You are an evolutionary biologist who is trying to find evidence of an evolutionary link between a newly discovered archaea from the dead sea region.  The archaea is a halophile and thermophile.  Your theory is that this archaea evolved from the bacterium when it acquired its salt and temperature resistance which enabled it to occupy new niches and evolve through adaptive radiation.  You have aggregated the genetic data from these species into a file called microbe_genes.faa (not aligned) and want to get a PCR reaction with salt and temperature conditions similar to those found in the archaean’s natural habitat.  The Sea of Salt is about 10 times saltier than regular ocean water and has a consistent temperature of 60 +/- 1 degrees Celsius because it is heated from geothermal activity.  You determine the concentration of salt in the Sea of Salt is about 6.3 mM and you want to use a primer concentration of 100 nM to be certain your primer will bond.

```console
./DeGenPrime --salt_conc:6.3 --primer_conc:100 --global --min_temp:59 --max_temp:61 microbe_genes.faa
```

## Copyright  
University of North Carolina at Charlotte, Bryan Fulghum, Sophie Tanker, and Richard Allen White III.  All rights reserved.  DeGenPrime is a bioinformatic tool that can be distributed freely.  
The software is provided “as is” and the copyright owners or contributors are not liable for any direct, indirect, incidental, special, or consequential damages including but not limited to, procurement of goods or services, loss of use, data or profits arising in any way out of the use of this software.<br />

## Citing DeGenPrime

If you are publishing results obtained using DeGenPrime, please cite: <br />

Fulghum BW, Tanker S, White III RA. 2023.  <br />
DeGenPrime provides robust primer design and optimization unlocking the biosphere. [bioRxiv](https://www.biorxiv.org/content/10.1101/2023.08.11.553048v1) <br />

## CONTACT

The informatics point-of-contact for this project is [Dr. Richard Allen White III](https://github.com/raw-lab).  
If you have any questions or feedback, please feel free to get in touch by email.  
[Dr. Richard Allen White III](mailto:rwhit101@uncc.edu)<br /> 
[Bryan Fulghum](mailto:bfulghu2@charlotte.edu
) <br />
Or [open an issue](https://github.com/raw-lab/degenprime/issues).  

[back to top](https://github.com/raw-lab/DeGenPrime/)
