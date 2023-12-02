# evenpurine
DNA resonance - Even Purine homology in trans positions - connecting DNA from distant locations

https://docs.google.com/document/d/e/2PACX-1vSXpUeaF5sF33RoXQa6BaUE-WMDrPC2zxTqrxYOkg7KZe_A01Rao-IOB5i3YcODB14XSXAU__oKgyJF/pub

The scripts have essentially the with the same principles: They align purine-code Alu bins against purine-code genome sequences. Each match will correspond to a protopair, which will be used to compute perhomol values % between the two sequences.

The first script will produce align alu bins that are splitted between odd and even phases.

The second script will produce alignments with alu bins that are not splitted in alu and even. The conversion will be done after we have the protopairs.

The third script will calculate perhomol based in a random pattern that do not correspond to the odd/even intervals. Each time you run, you have a new random pattern being used.
