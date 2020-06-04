# PBC_plot

**Beta submission! Includes py3 script for per base coverage plot generation. Will be shortly followed by nextflow pipeline for mapping reads to reference, processing and obtaining the input files for the .py plotting script.**

Per base coverage statistics pipeline to assess the expected versus factual coverage of the reference sequence region by the newly sequenced reads (arrival-like statistics). Plots the per base coverage, helping visualize the possible sequencing\mapping pitfalls and repetative region collapses.



_______________________________


- Arrival Statistics is calculated for each position in the reference genome (can be closely related species, if absent). For this the per base coverage file is generated with bedtools genecoverage, based on the reads mapped to reference (using minimap2, e.g.).
- To calculate **mean coverage for short reads** the following formula is used: 
```
c=nreads*lreads/lgen
```
where **nreads** - n of mapped to ref reads, 
      **lreads** - length of short reads(hard-coded in beta to 151 in script),
      **lgene** - length of reference sequence, e.g.genome
      
- To calculate **mean coverage for long reads** the following formula is used:
```
c=sum of lengths of all mapped reads/length of reference sequence
```

___________________________

The script was first designed to indirectly detect and visualize potential NUMTs-derived contaminants in mitochondrial reads. While the arrival statistics is mostly used to detect the collapsed repetative regions, in case of mitochondrial genomes the spike in coverage of particular loci in mitogenome may be evident of the nuclear pseudogenic derived reads, that may  impact the adequacy of mitochondrial genome assembly.
