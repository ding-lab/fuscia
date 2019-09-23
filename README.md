# fuscia
Detection of chimeric transcripts (aka fusions) from barcoded single cell RNA-seq

## Usage

### detect_discordant_reads.py
```{python}
python detect_discordant_reads.py star-fusion_output_file scRNA-seq.bam output_dir output_prefix region_plusminus
```
Uses STAR-Fusion output from bulk RNA-seq as starting point for where to look for fusion transcripts in scRNA-seq.

### discover_discordant_reads.py
```{python}
python discover_discordant_reads.py scRNA-seq.bam output_dir chrA:pos-pos chrB:pos-pos output_dir output_prefix
```
You define regions for discovering fusion transcripts in scRNA-seq. Not based on STAR-Fusion results.

## Output
For each read coming from a chimeric transcript, output is a tab-separated file with columns

* cell_barcode (cell the read originated from)
* molecular_barcode (molecule the read originated from)
* chrom (chromosome of read)
* start (start of read)
* end (end of read)

## Imports
os
pysam
sys
