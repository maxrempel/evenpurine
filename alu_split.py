from Bio import SeqIO
import os
import argparse
import subprocess
import pandas as pd
import time
import warnings
from Bio import BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=BiopythonWarning)
    from Bio import pairwise2
from Bio.Seq import Seq
import re
from io import StringIO

start_time = time.time()

parser = argparse.ArgumentParser(description='This script has been designed to convert two sequences to Purine code, and align them to each other')
parser._action_groups.pop()
required = parser.add_argument_group('Required arguments')
required.add_argument('--genome', help='Genome in fasta', required=True, type=str, metavar = "")
required.add_argument('--te', help='TE sequence in one-line fasta format', required=True, type=str, metavar = "")

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('--bin_size', help='Bin size bp to split TE sequence. Default = 60', required=False, type=str, default=60, metavar = "")
optional.add_argument('--overlap', help='Overlap  between bins. Default = 50', required=False, type=str, default=50, metavar = "")
optional.add_argument('--perhomol', help='Homology between bins  and genome sequence. Default = 70', required=False, type=str, default=70, metavar = "")
optional.add_argument('--match_length', help='Minimum bin length to match with genome. Default = 60', required=False, type=str, default=60, metavar = "")
parser.parse_args()
args = parser.parse_args()

def rm_file(file):
    if os.path.exists(str(file)):
        os.remove(str(file)) #If the seq Pcode already exists, delete it
def Pcoding(fasta_IN, fasta_OUT): #Transform sequences to Pcode
    rm_file(str(fasta_OUT))

    dictio = str.maketrans({'A': 'G', 'T': 'C'}) # define the purine code
    with open(str(fasta_IN), 'r') as file:
        sequence = file.read().splitlines()

    with open(str(fasta_OUT), 'a') as f:
        for line in sequence:
            if line.startswith(">"):
                print(str(line), file=f)
            if not line.startswith(">"):
                ptrans = line.translate(dictio)
                print(str(ptrans), file=f)

print("Converting genome and TE to purine-code...")
Pcoding(str(args.genome), str('hg38_chr20_Pcode.fa'))
Pcoding(str(args.te), str('aluY_Pcode.fa'))
print("Done!")

def split_odd_even(fasta_IN, fasta_OUT):
    with open(str(fasta_IN), 'r') as file:
        seq = file.read().splitlines()
        for line in seq:
         if not line.startswith('>'):
            odd_chars = line[::2]  # Extract characters at odd positions
            with open(str(fasta_OUT + '_epod.fa'), 'w') as f:
                print(str(">" + fasta_OUT + "_epod" + "\n" + odd_chars), file=f)

            even_chars = line[1::2]  # Extract characters at even positions
            with open(str(fasta_OUT + '_epev.fa'), 'w') as f:
                print(str(">" + fasta_OUT + "_epev" + "\n" + even_chars), file=f)

print("Splitting into odd and even purine-code sequences...")
split_odd_even(str('hg38_chr20_Pcode.fa'), str('genome'))
split_odd_even(str('aluY_Pcode.fa'), str('TE'))
print("Done!")

def split_sequence(sequence, bin_size):
    overlap_size = int(bin_size * (int(args.overlap) / 100))
    bins = [sequence[i:i + bin_size] for i in range(0, len(sequence) - bin_size + 1, bin_size - overlap_size)]
    return bins

def create_bins(fasta_file, bin_size):
    bins = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        bins.extend(split_sequence(sequence, bin_size))
    return bins

perc_overlap = int(args.overlap) / 100
bin_overlap_length = int(perc_overlap * int(args.bin_size))
perc_notoverlap = 1 - perc_overlap
bin_NOToverlap_length = int(perc_notoverlap * int(args.bin_size))

print("Splitting TE sequece into bins...")
bins = create_bins(str('TE_epev.fa'), int(args.bin_size))
rm_file(str('TE_epev_bins.fa'))
with open('TE_epev_bins.fa', 'a') as f:
    for i, bin_sequence in enumerate(bins):
        if i == 0:
            print(str(">Bin_epev_range") + str(args.bin_size), file=f)
            print(bin_sequence, file=f)
        else:
            print(f">Bin_epev_range{(bin_NOToverlap_length * i) + int(args.bin_size)}", file=f)
            print(bin_sequence, file=f)

bins = create_bins(str('TE_epod.fa'), int(args.bin_size))
rm_file(str('TE_epod_bins.fa'))
with open('TE_epod_bins.fa', 'a') as f:
    for i, bin_sequence in enumerate(bins):
        if i == 0:
            print(str(">Bin_epod_range") + str(args.bin_size), file=f)
            print(bin_sequence, file=f)
        else:
            print(f">Bin_epod_range{(bin_NOToverlap_length * i) + int(args.bin_size)}", file=f)
            print(bin_sequence, file=f)
print("Done!")

##### Align bins to the genome with blastn-short
def bin_alignment(phase):
    with open(str(phase + '_blast.tsv'), 'w') as table:
        subprocess.call(['blastn', '-task', 'blastn-short', '-query', str('TE_' + phase + '_bins.fa'),
                         '-subject', 'genome_epod.fa', '-outfmt', '6', '-word_size', '7', '-perc_identity', str(args.perhomol),
                         '-qcov_hsp_perc', str(args.match_length), '-ungapped'], stdout=table)
    table.close

print("Performing epev and opod bins alignment...")
bin_alignment(str('epev'))
bin_alignment(str('epod'))
print("Done!")

def extract_sequence(fasta_file, start_pos, end_pos):
    with open(fasta_file, "r") as file:
        # Parse the FASTA file
        records = list(SeqIO.parse(file, "fasta"))

    # Assuming there is only one record in the file
    if len(records) == 1:
        sequence = records[0].seq
        # Extract the desired subsequence
        subsequence = sequence[start_pos - 1:end_pos]
        return subsequence
    else:
        raise ValueError("Input FASTA file should contain exactly one sequence.")

def regex_1stcolumn(blast_table):
    data = pd.read_csv(str(blast_table), sep = '\t', header = None)
    data.iloc[:, 0] = data.iloc[:, 0].replace(to_replace=r'.*range', value='', regex=True)
    return data

table_epev = regex_1stcolumn(str('epev_blast.tsv'))
table_epod = regex_1stcolumn(str('epod_blast.tsv'))

# Create ranges for masked regions
with open(str(args.genome), 'r') as f:
    sequence = f.read().splitlines()
    for line in sequence:
        if not line.startswith(">"):
            genome_seq = str(line)
seq_object = Seq(genome_seq)

# Find masked regions
masked_regions = [(match.start() + 1, match.end()) for match in re.finditer(r'N+', str(seq_object))]
# Print the masked regions
with open('masked_regions.txt', 'w') as f:
    for start, end in masked_regions:
        if start > 10:
            print(f"{start - 10}\t{end + 10}", file=f)
        else:
            print(f"{start}\t{end + 10}", file=f)
f.close()



#Masked regions to bed
masked_ranges = pd.read_csv(str('masked_regions.txt'), usecols=[0,1], names=['mask_start','mask_end'], sep='\t', header=None)
masked_ranges['chr'] = 'chr20'
masked_ranges = masked_ranges[['chr', 'mask_start', 'mask_end']]
masked_ranges.to_csv('masked_regions.bed', sep='\t', header=False, index=False)

#epev and epod intervals to bed
epev_intervals = pd.read_csv(str('epev_blast.tsv'), usecols=[8,9], names=['match_start','match_end'], sep='\t', header=None)
epev_intervals[['match_start', 'match_end']] = epev_intervals.apply(lambda row: sorted(row), axis=1, result_type='expand')
epev_intervals *= 2
epev_intervals['chr'] = 'chr20'
epev_intervals = epev_intervals[['chr', 'match_start', 'match_end']]
epev_intervals.to_csv('epev_matches.bed', sep='\t', header=False, index=False)

#epev and epod intervals to bed
epev_intervals = pd.read_csv(str('epod_blast.tsv'), usecols=[8,9], names=['match_start','match_end'], sep='\t', header=None)
epev_intervals[['match_start', 'match_end']] = epev_intervals.apply(lambda row: sorted(row), axis=1, result_type='expand')
epev_intervals *= 2
epev_intervals['chr'] = 'chr20'
epev_intervals = epev_intervals[['chr', 'match_start', 'match_end']]
epev_intervals.to_csv('epod_matches.bed', sep='\t', header=False, index=False)

epev_masked_region = subprocess.check_output(['bedtools', 'intersect', '-a', 'masked_regions.bed', '-b', 'epev_matches.bed', '-wa', '-wb'],
                                   universal_newlines=True)
epev_masked_region = pd.read_csv(StringIO(epev_masked_region), usecols=[3,4,5],sep='\t') # 3 column pandas df: chr start end

epod_masked_region = subprocess.check_output(['bedtools', 'intersect', '-a', 'masked_regions.bed', '-b', 'epod_matches.bed', '-wa', '-wb'],
                                   universal_newlines=True)
epod_masked_region = pd.read_csv(StringIO(epod_masked_region), usecols=[3,4,5],sep='\t') # 3 column pandas df: chr start end

def split_odd_even_propairs(seq_in): #, fasta_OUT):
    odd_chars = seq_in[::2]  # Extract characters at odd positions
    #with open(str(fasta_OUT + '_epod.fa'), 'w') as f:
        #print(str(">" + fasta_OUT + "_epod" + "\n" + odd_chars), file=f)

    even_chars = seq_in[1::2]  # Extract characters at even positions
    #with open(str(fasta_OUT + '_epev.fa'), 'w') as f:
        #print(str(">" + fasta_OUT + "_epev" + "\n" + even_chars), file=f)
    return odd_chars, even_chars

def calculate_homology(sequence1, sequence2):
    # Create Seq objects from the input sequences
    seq1 = Seq(sequence1)
    seq2 = Seq(sequence2)

    # Perform pairwise alignment
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)

    # Extract the aligned sequences
    aligned_seq1 = alignments[0][0]
    aligned_seq2 = alignments[0][1]

    # Calculate homology as a percentage
    homology_percentage = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1) * 100

    return homology_percentage

def alignment_protopairs(pandas_blast):
    with open(str("Protopairs.tsv"), "a") as output_table:
        for row in pandas_blast.itertuples():
            bin_end = row[1]
            bin_start = int(bin_end) - int(args.bin_size)

            ### GENOME SEQUENCE
            genome_start = row[9] * 2
            genome_end = row[10] * 2

            if str(pandas_blast) == str(table_epev):
                alu_mask_check = epev_masked_region[(epev_masked_region.iloc[:, 1] == genome_start) & (epev_masked_region.iloc[:, 2] == genome_end)]
            else:
                alu_mask_check = epod_masked_region[(epod_masked_region.iloc[:, 1] == genome_start) & (epod_masked_region.iloc[:, 2] == genome_end)]

            if alu_mask_check.empty:
                if genome_start > genome_end:
                    genome_start, genome_end = genome_end, genome_start

                genomic_sequence = extract_sequence(str(args.genome), genome_start, genome_end)

                genomic_ptrans = str(genomic_sequence).translate(dictio)
                genome_epod, genome_epev = split_odd_even_propairs(str(genomic_ptrans))


                ### ALU SEQUENCE
                match_start = row[7]
                match_end = row[8]
                #start_mask_check = match_start * 2
                #end_mask_check = match_end * 2

                alu_phased_start = match_start + bin_start
                alu_phased_end = match_end + bin_start

                alu_original_start = alu_phased_start * 2
                alu_original_end = alu_phased_end * 2
                alu_sequence = extract_sequence(str(args.te), alu_original_start, alu_original_end)

                alu_ptrans = str(alu_sequence).translate(dictio)
                TE_epod, TE_epev = split_odd_even_propairs(str(alu_ptrans))

                homology_epod = calculate_homology(str(TE_epod), str(genome_epod))
                homology_epev = calculate_homology(str(TE_epev), str(genome_epev))
                print(f"{homology_epod:.2f} \t {homology_epev:.2f}", file=output_table)
            else:
                print(alu_mask_check)
    output_table.close()

            ###COMPARING PROTOPAIRS
            # epod_res = subprocess.check_output(['blastn', '-task', 'blastn-short', '-query', 'tmp_alu_region_epod.fa',
            #                                     '-subject', 'tmp_genomic_region_epod.fa', '-outfmt', '6 length', '-word_size', '7', '-ungapped','-max_target_seqs', '1', '-max_hsps', '1'],
            #                                    universal_newlines=True)
            # match_length_epod = (int(epod_res) * 100) / int(args.bin_size)
            #
            # epev_res = subprocess.check_output(['blastn', '-task', 'blastn-short', '-query', 'tmp_alu_region_epev.fa',
            #                                     '-subject', 'tmp_genomic_region_epev.fa', '-outfmt', '6 length', '-word_size', '7',
            #                                     '-ungapped', '-max_target_seqs', '1', '-max_hsps', '1'],
            #                                    universal_newlines=True)
            # match_length_epev = (int(epev_res) * 100) / int(args.bin_size)
            # print(str(match_length_epev) + "\t" + str(match_length_epod), file=output_table)

print("Calculating homology...")
rm_file(str('Protopairs.tsv'))
dictio = str.maketrans({'A': 'G', 'T': 'C'})  # define the purine code
alignment_protopairs(table_epev)
alignment_protopairs(table_epod)
print("Done!")
#
# end_time = time.time()
# elapsed_time = end_time - start_time
# print(f"Processing time: {elapsed_time:.2f} seconds")