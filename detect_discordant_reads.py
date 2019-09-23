# Script to detect discordant reads from 10X single cell RNA-seq data
# Usages: python detect_discordant_reads.py star-fusion_output_file scRNA-seq.bam output_dir output_prefix region_plusminus

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
star_fusion_file = sys.argv[1]
bam_filename = sys.argv[2]
output_dir = sys.argv[3]
output_prefix = sys.argv[4]
plusminus = int(sys.argv[5])
min_mapq = int(sys.argv[6])

star_fusion = open(star_fusion_file, "r")
discordant_reads_list = []
for line in star_fusion:
  if line.startswith("#"):
    next
  else:
    FusionName, JunctionReadCount, SpanningFragCount, SpliceType, LeftGene, LeftBreakpoint, RightGene, RightBreakpoint = line.strip().split()[0:8]
    chromA = LeftBreakpoint.split(":")[0][3:]
    min_posA = int(LeftBreakpoint.split(":")[1]) - plusminus
    max_posA = int(LeftBreakpoint.split(":")[1]) + plusminus
    chromB = RightBreakpoint.split(":")[0][3:]
    min_posB = int(RightBreakpoint.split(":")[1]) - plusminus
    max_posB = int(RightBreakpoint.split(":")[1]) + plusminus
    
    # main business 
    rangeA_reads = {}
    rangeB_reads = {}

    # iterate over rangeA
    samfile = pysam.AlignmentFile(bam_filename, "rb")
    for read in samfile.fetch(chromA, min_posA, max_posA):
      read_info = extract_read_info(read, min_mapq)
      if read_info == None:
        next
      else:
        CB = read_info["CB"] 
        UB = read_info["UB"]
        if CB+":"+UB in rangeA_reads:
          rangeA_reads[CB+":"+UB].append(read_info)
        else:
          rangeA_reads[CB+":"+UB] = [read_info]   
    samfile.close()

    # iterate over rangeB
    samfile = pysam.AlignmentFile(bam_filename, "rb")
    for read in samfile.fetch(chromB, min_posB, max_posB):
      read_info = extract_read_info(read, min_mapq)
      if read_info == None:
        next
      else:
        CB = read_info["CB"]
        UB = read_info["UB"]
        if CB+":"+UB in rangeB_reads:
          rangeB_reads[CB+":"+UB].append(read_info)
        else:
          rangeB_reads[CB+":"+UB] = [read_info]
    samfile.close()

    # check if any overlap
    for key in rangeA_reads:
      if key in rangeB_reads:
        combined_reads = rangeA_reads[key]
        combined_reads.extend(rangeB_reads[key])
        for read_info in combined_reads:
          discordant_reads_list.append([str(x) for x in [read_info["CB"], read_info["UB"], read_info["this_chromosome"], read_info["this_start_bp"], read_info["this_end_bp"], FusionName]])

unique_discordant_reads_list = []
for x in discordant_reads_list:
  if x not in unique_discordant_reads_list:
    unique_discordant_reads_list.append(x)

os.makedirs(output_dir, exist_ok = True)
output_file_path = os.path.join(output_dir, output_prefix + ".discordant_reads.tsv")
output_file = open(output_file_path, "w")
output_file.write("\t".join(["cell_barcode", "molecular_barcode", "chrom", "start", "end", "fusion"]) + "\n") # list of column headers
for dis_read in sorted(unique_discordant_reads_list):
  output_file.write("\t".join(dis_read) + "\n") # write info for each discordant read
output_file.close() 
