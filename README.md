# DeGenPrime

[Copyright](###Copyright) <br />
[Introduction](###Introduction) <br />
[Citation](###Citation) <br />
[Installation](###Installation) <br />
[Options](###Options)<br />
[Input](###Input)<br />
[Output](###Output) <br />
[Examples](###Examples) <br />

### Copyright  <br />
University of North Carolina at Charlotte, Bryan Fulghum, Sophie Tanker, and Richard Allen White III.  All rights reserved.  DeGenPrime is a bioinformatic tool that can be distributed freely.  
The software is provided “as is” and the copyright owners or contributors are not liable for any direct, indirect, incidental, special, or consequential damages including but not limited to, procurement of goods or services, loss of use, data or profits arising in any way out of the use of this software.<br />

### Introduction <br />
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

### Citation  <br />
DeGenPrime is free software to use. <br /> 
If you use it however, we ask that you please cite the software in publications with:<br />
Fulghum B, Tanker S, White RA III. DeGenPrime – Robust Degenerate Primer Design for Unlocking the Diversity of the Biosphere.  <br />

### Installation   <br />
DeGenPrime is designed to run as a stand-alone console application on any platform capable of running C++ applications.  The program does try to align a file that is misaligned by calling MAFFT within the program, so you must either have this installed or manually align your sequences to use DeGenPrime.

### Options <br />

```
./DeGenPrime [--tags] <filename> (the filename must always be the last argument, or this program will throw a segmentation error.)

Command-line arguments
Valid tags include:
--amplicon:           <int>, Set the minimum amplicon length. <br />
--begin:              <int>, Set the beginning nucleotide.  If the user inputs an integer < 0, the program will default this value to zero.<br />
--end:                <int>, Set the ending nucleotide.  If the user inputs an integer > the number of base pairs in the entire sequence, then the program will default this to the last nucleotide in the sequence.<br />
--global or --g,       for lists of sequences that are misaligned, this tag specifies that the file should run MAFFT for global alignment.<br />
--help or --h,         prints this help menu.<br />
--local or --l,        for lists of sequences that are misaligned, this tag specifies that the file should run MAFFT for local alignment.<br />
--min_temp:           <int>, Sets the minimum primer melting temperature.  This has a minimum value of 50.0 (degrees Celsius) and must be smaller than --max_temp.<br />
--max_temp:           <int>, Sets the maximum primer melting temperature.  This has a maximum value of 65.0 (degrees Celsius) and must be larger than --min_temp.<br />
--primer_conc:        <int>, Sets the concentration of the PCR primer in nM.  This has a minimum value of 50.0 nM, and this program will raise any value smaller to this value.<br />
--protein,             Tells the program that the input sequence is a protein sequence, and the program should unwrap the protein sequence into its base nucleotides instead of trying to find a PCR.  This will produce degenerate nucleotide codes whenever there is any ambiguity.<br />
--salt_conc:          <int>, Sets the concentration of monovalent ions in mM.  This has a minimum value of 50.0 mM, and this program will raise any value smaller to this value.<br />
--max_primers:        <int>, Sets the maximum number of output primers.  This has a maximum value of 10 and this program will reduce any value larger to this value.<br />
```

### Input
Input files to DeGenPrime can use a variety of formats including single sequences to alignment files.<br /> 
- Nucleotide fasta formats (.fasta, .fna, .ffn)<br />
- Protein fasta formats (.faa)<br />
- Protein or nucleotide alignments (.clust)<br />

### Output
DeGenPrime will output a few progress messages as the program runs, the recommended primers and their details, and the program runtime at the end to the console.  DeGenPrime also records the operation details and the primers to an output file.  The output filename is selected based on the specified input file.<br />  

#### Standard outputs
- primers_filename.txt (a list of primers) <br />
- Any errors that occur during runtime will also be written on the console or to stout.<br />

### Examples
You are a federal prosecutor trying to build a case against Dr. Scientist who stands accused of performing illegal genetic engineering experiments on human zygotes.  Your team believes that they made genetic modifications on chromosome 21 in the p arms above the centromere.  You find out the length of this region from the tips of the 3’ end of the telomeres to the centromere is about 12 million base pairs.  The average length of telomeres on chromosome 21 for a human zygote is 10 thousand base pairs.  Your team collects data from the allegedly modified zygotes’ chromosome 21 and other genetic information collected from NCBI on chromosome 21 into a fasta file called zygote_21.faa.  You want to find 10 good primer pairs on this chromosome to amplify in a PCR reaction.<br /> You would use:<br />

```
./DeGenPrime --begin:10000 --end:12000000 --max_primers:10 zygote_21.faa
```

You are an immunologist who wants to identify mutations that might have occurred in the genes of a local strain of the influenza A virus so you can produce a new vaccine for the upcoming flu season.  Microbiologists have reported changes to the geometry of the surface proteins of infected cells and this leads you to hypothesize that a mutation has occurred in segment 4 or segment 6 of the virus.  You obtain the genetic data of these segments from last season’s influenza A as well as data from the current strain into two files called influenzaA_4.faa and influenzaA_6.faa.  You know that the surface proteins normally contain 500 +/- 50 amino acids (which implies 1500 +/- 150 bps in the coding sequence) and you want the top 5 primers for a PCR reaction that will cover at least 40% of the coding region.  Your files are not aligned. <br />  
You would use:<br />

```
./DeGenPrime --amplicon:660 --local influenzaA_4.faa or influenzaA_6.faa
```

You are an evolutionary biologist who is trying to find evidence of an evolutionary link between a newly discovered archaea from the Sea of Salt and primordial bacterium from that region.  The archaea is a halophile and thermophile.  Your theory is that this archaea evolved from the bacterium when it acquired its salt and temperature resistance which enabled it to occupy new niches and evolve through adaptive radiation.  You have aggregated the genetic data from these species into a file called microbe_genes.faa (not aligned) and want to get a PCR reaction with salt and temperature conditions similar to those found in the archaean’s natural habitat.  The Sea of Salt is about 10 times saltier than regular ocean water and has a consistent temperature of 60 +/- 1 degrees Celsius because it is heated from geothermal activity.  You determine the concentration of salt in the Sea of Salt is about 6.3 mM and you want to use a primer concentration of 100 nM to be certain your primer will bond.

```
./DeGenPrime --salt_conc:6.3 --primer_conc:100 --global --min_temp:59 --max_temp:61 microbe_genes.faa
```

[back to top](https://github.com/raw-lab/DeGenPrime/edit/main/README.md)
