# evenpurine
DNA resonance - Even Purine homology in trans positions - connecting DNA from distant locations

https://docs.google.com/document/d/e/2PACX-1vSXpUeaF5sF33RoXQa6BaUE-WMDrPC2zxTqrxYOkg7KZe_A01Rao-IOB5i3YcODB14XSXAU__oKgyJF/pub

The scripts have essentially the with the same principles: They align purine-code Alu bins against purine-code genome sequences. Each match will correspond to a protopair, which will be used to compute perhomol values % between the two sequences.

The `alu_split.py` script will produce align alu bins that are splitted between odd and even phases.

The `alu_unsplit.py` script will produce alignments with alu bins that are not splitted in alu and even. The conversion will be done after we have the protopairs.

The `alu_unsplit_RANDOM.py` script will calculate perhomol based in a random pattern that do not correspond to the odd/even intervals. Each time you run, you have a new random pattern being used.

The three scripts have the same usage:

````
Required arguments:
  --genome         Genome in fasta
  --te             TE sequence in one-line fasta format

Optional arguments:
  --bin_size       Bin size bp to split TE sequence. Default = 60
  --overlap        Overlap between bins. Default = 50
  --perhomol       Homology between bins and genome sequence. Default = 70
  --match_length   Minimum bin length to match with genome. Default = 60

````

Therefore, to run your analysis, just make:

````
python alu_unsplit_RANDOM.py --genome hg38_chr20.fa --te aluY_seq.fa
````

If you want to change bin size, perhomol or any other parameters, you must to include them in your commandline.
