# fuscia
Detection of chimeric transcripts (aka fusions) from barcoded single cell RNA-seq

## Usage

The main idea of the software is to find transcripts with reads that map to different locations. If reads from a single transcript (defined as having the same cell barcode and molecular barcode) map to different locations, we say that transcript is a "chimeric transcript". Chimeric transcripts may correspond to real fusion events!

### detect_discordant_reads.py
```{python}
python detect_discordant_reads.py star-fusion_output_file scRNA-seq.bam output_dir output_prefix region_plusminus min_mapq
```
Uses STAR-Fusion output from bulk RNA-seq as starting point for where to look for fusion transcripts in scRNA-seq.

#### Positional command line arguments:

* star-fusion_output_file (STAR-Fusion output file from your bulk sample -- e.g. star-fusion.fusion_predictions.tsv)
* scRNA-seq.bam (single cell RNA-seq bam from your sample -- note: corresponding bam.bai must also exist in same directory as .bam file)
* output_dir (full or relative path to output directory)
* output_prefix (output filename prefix to distinguish results from other samples)
* region_plusminus (distance up/downstream from fusion breakpoint to look for discordant reads -- for a 2x base pair window, set region_plusminus to be x base pairs)
* min_mapq (minimum mapping quality required for reads to be analyzed)

### discover_discordant_reads.py
```{python}
python discover_discordant_reads.py scRNA-seq.bam chrA:pos-pos chrB:pos-pos output_dir output_prefix min_mapq
```
You define regions for discovering fusion transcripts in scRNA-seq. Not based on STAR-Fusion results.

#### Positional command line arguments:

* scRNA-seq.bam (single cell RNA-seq bam from your sample -- note: corresponding bam.bai must also exist in same directory as .bam file))
* chrA:pos-pos (chromosome A base pair range where program will look for discordant reads)
* chrB:pos-pos (chromosome B base pair range where program will look for discordant reads)
* output_dir (full or relative path to output directory)
* output_prefix (output filename prefix to distinguish results from other samples)
* min_mapq (minimum mapping quality required for reads to be analyzed)

## Output
The output file has one line for each read coming from a chimeric transcript.

For detect_discordant_reads.py, the output file is  `{output_dir}/{output_prefix}.discordant_reads.tsv`.
For discover_discordant_reads.py, the output file is `{output_dir}/{output_prefix}.discovered_discordant_reads.tsv`.

### Output file column names

* cell_barcode (cell the read originated from)
* molecular_barcode (molecule the read originated from)
* chrom (chromosome of read)
* start (start of read base pair mapping position)
* end (end of read base pair mapping position)
* fusion (detect_discordant_reads.py only -- corresponding bulk RNA-seq fusion called by STAR-Fusion)

## Imports
os
pysam
sys

## Author information
Steven Foltz ([github: envest](https://github.com/envest))
