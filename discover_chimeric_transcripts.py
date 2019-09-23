# Script to discover transcripts with reads mapping to different chromosomes from 10X single cell RNA-seq data
# Usages: python detect_discordant_reads.py scRNA-seq.bam output_dir chrA:pos-pos chrB:pos-pos output_dir output_prefix

import os
import pysam
import sys

def extract_read_info(read, min_mapq):
  # Given a read from region A, extract relevant information if mate maps to region B, else None
  # Relevant information includes: cell barcode (CB), molecular barcode (UB), and sample index read (BC)
  if read.is_duplicate or read.is_qcfail or read.is_secondary or read.mapping_quality < min_mapq: #or not read.is_proper_pair:
    return(None)
  elif read.has_tag("CB") and read.has_tag("UB"): #and read.has_tag("BC"):
    return(return_read_info(read))
  else: 
    return(None)

def return_read_info(read):
  # If a read passes quality filters and has CB, UB, and BC tags:
  # Return a dictionary of relevant read information, including
  # Mapping position of read and mate, all barcode information
  return_dict = {}  
  return_dict["this_chromosome"] = str(read.reference_name)
  return_dict["this_start_bp"] = int(read.reference_start)
  return_dict["this_end_bp"] = int(read.reference_end)
  tags_dict = {x:y for (x,y) in read.get_tags()}
  return_dict["CB"] = tags_dict["CB"]
  return_dict["UB"] = tags_dict["UB"]
  return(return_dict)

# read arguments
bam_filename = sys.argv[1]
chrA_range = sys.argv[2]
chrB_range = sys.argv[3]
output_dir = sys.argv[4]
output_prefix = sys.argv[5]
min_mapq = int(sys.argv[6])

chrA = str(chrA_range.split(":")[0])
chrA_minpos = int(chrA_range.split(":")[1].split("-")[0])
chrA_maxpos = int(chrA_range.split(":")[1].split("-")[1])
chrB = str(chrB_range.split(":")[0])
chrB_minpos = int(chrB_range.split(":")[1].split("-")[0])
chrB_maxpos = int(chrB_range.split(":")[1].split("-")[1])

chrA_reads = {}
chrB_reads = {}

# iterate over chrA
samfile = pysam.AlignmentFile(bam_filename, "rb")
for read in samfile.fetch(chrA, chrA_minpos, chrA_maxpos):
  read_info = extract_read_info(read, min_mapq)
  if read_info == None:
    next
  else:
    CB = read_info["CB"] 
    UB = read_info["UB"]
    if CB+":"+UB in chrA_reads:
      if read_info["this_chromosome"] in chrA_reads[CB+":"+UB]:
        chrA_reads[CB+":"+UB][read_info["this_chromosome"]] += 1
      else:
        chrA_reads[CB+":"+UB][read_info["this_chromosome"]] = 1
    else:
      chrA_reads[CB+":"+UB] = {read_info["this_chromosome"]: 1}  
samfile.close()

print("through chrA first time")

# iterate over chrB
samfile = pysam.AlignmentFile(bam_filename, "rb")
for read in samfile.fetch(chrB, chrB_minpos, chrB_maxpos):
  read_info = extract_read_info(read, min_mapq)
  if read_info == None:
    next
  else:
    CB = read_info["CB"]
    UB = read_info["UB"]
    if CB+":"+UB in chrB_reads:
      if read_info["this_chromosome"] in chrB_reads[CB+":"+UB]:
        chrB_reads[CB+":"+UB][read_info["this_chromosome"]] += 1
      else:
        chrB_reads[CB+":"+UB][read_info["this_chromosome"]] = 1
    else:
      chrB_reads[CB+":"+UB] = {read_info["this_chromosome"]: 1}
samfile.close()

print("through chrB first time")

# check if any overlap
n_overlapping = 0
overlapping_barcodes = {}
for key in chrA_reads:
  if key in chrB_reads:
    n_overlapping += 1
    overlapping_barcodes[key] = None

print("found "+str(n_overlapping)+" overlapping barcodes")

discordant_reads_list = []
    
# get reads from chimeric transcript
# iterate over chrA (again)
samfile = pysam.AlignmentFile(bam_filename, "rb")
for read in samfile.fetch(chrA, chrA_minpos, chrA_maxpos):
  read_info = extract_read_info(read, min_mapq)
  if read_info == None:
    next
  else:
    CB = read_info["CB"]
    UB = read_info["UB"]
    if CB+":"+UB in overlapping_barcodes:
      discordant_reads_list.append([str(x) for x in [read_info["CB"], read_info["UB"], read_info["this_chromosome"], read_info["this_start_bp"], read_info["this_end_bp"]]])
samfile.close()

print("through chrA second time")

# iterate over chrB
samfile = pysam.AlignmentFile(bam_filename, "rb")
for read in samfile.fetch(chrB, chrB_minpos, chrB_maxpos):
  read_info = extract_read_info(read, min_mapq)
  if read_info == None:
    next
  else:
    CB = read_info["CB"]
    UB = read_info["UB"]
    if CB+":"+UB in overlapping_barcodes:
      discordant_reads_list.append([str(x) for x in [read_info["CB"], read_info["UB"], read_info["this_chromosome"], read_info["this_start_bp"], read_info["this_end_bp"]]])
samfile.close()

print("through chrB second time")

unique_discordant_reads_list = []
for x in discordant_reads_list:
  if x not in unique_discordant_reads_list:
    unique_discordant_reads_list.append(x)

os.makedirs(output_dir, exist_ok = True)
output_file_path = os.path.join(output_dir, output_prefix + ".discovered_discordant_reads.tsv")
output_file = open(output_file_path, "w")
output_file.write("\t".join(["cell_barcode", "molecular_barcode", "chrom", "start", "end"]) + "\n") # list of column headers
for dis_read in sorted(unique_discordant_reads_list):
  output_file.write("\t".join(dis_read) + "\n") # write info for each discordant read
output_file.close()
