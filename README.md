# fuscia
Detection of chimeric transcripts (aka fusions) from barcoded single cell RNA-seq

## Usage

### detect_discordant_reads.py
```{python}
python detect_discordant_reads.py star-fusion_output_file scRNA-seq.bam output_dir output_prefix region_plusminus
```
Uses STAR-Fusion output from bulk RNA-seq as starting point for where to look for fusion transcripts in scRNA-seq.

Command line arguments:

* star-fusion_output_file (STAR-Fusion output file from your bulk sample)
* scRNA-seq.bam (single cell RNA-seq bam from your sample)
* output_dir (full or relative path to output directory)
* output_prefix (output filename prefix to distinguish results from other samples)
* region_plusminus (distance up/downstream from fusion breakpoint to look for discordant reads -- for a 2x base pair window, set region_plusminus to be x base pairs)
* min_mapq (minimum mapping quality required for reads to be analyzed)

### discover_discordant_reads.py
```{python}
python discover_discordant_reads.py scRNA-seq.bam chrA:pos-pos chrB:pos-pos output_dir output_prefix
```
You define regions for discovering fusion transcripts in scRNA-seq. Not based on STAR-Fusion results.

Command line arguments:

* scRNA-seq.bam (single cell RNA-seq bam from your sample)
* chrA:pos-pos (chromosome A base pair range)
* chrB:pos-pos (chromosome B base pair range)
* output_dir (full or relative path to output directory)
* output_prefix (output filename prefix to distinguish results from other samples)
* min_mapq (minimum mapping quality required for reads to be analyzed)

## Output
For each read coming from a chimeric transcript, output is a tab-separated file with columns

* cell_barcode (cell the read originated from)
* molecular_barcode (molecule the read originated from)
* chrom (chromosome of read)
* start (start of read)
* end (end of read)
* fusion (detect_discordant_reads.py only -- corresponding bulk RNA-seq fusion called by STAR-Fusion)

For detect_discordant_reads.py, the output file ends with `.discordant_reads.tsv`.
For discover_discordant_reads.py, the output file ends with `.discovered_discordant_reads.tsv`.

## Imports
os
pysam
sys
