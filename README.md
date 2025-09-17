# CpG Island Prediction Program
The program `dinucleotide.py` calculates the frequency of each dinucleotide, compares each dinucleotide’s observed and expected frequencies, and prints out dinucleotides with a significant difference between their observed and expected frequencies. It also models CpG regions via Markov chains, and predicts CpG islands with evaluation metrics (true positive, false positive, false negative).

The general command format is as follows:
```
python dinucleotide.py <fasta_input> -f -hmm <intervals_input>
```
The toggle options are as follows:

| Option        | Description   | 
| ------------- |:-------------| 
| `-f`          | Print dinucleotide frequency summary, a comparison table between observed and expected frequencies, and a list of significant nucleotides. | 
| `-hmm <interval_input>`            | Prints statistics of predicted CpG islands based on the CpG potential and if each predicted island is a true positive, false negative, or false positive. Results are outputted to `results.txt.`      | 

This program deems a dinucleotide to be “significantly different” if the difference between its observed and expected frequencies is greater than 0.02.

The included test files `chrA.fasta` and `chrA.islands` can be used for `<fasta_input>` and `<intervals_input>` respectively, to run `dinucleotide.py`.
