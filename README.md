# Bacaro's Beta
Beta Diversity calculation of 16S read samples based on Bacaro, Ricotta, &amp; Mazzoleni (2007).

## About
### Approach overview
A simple implementation of the Beta diversity metric proposed in [Measuring beta-diversity from taxonomic similarity](https://doi.org/10.1111/j.1654-1103.2007.tb02595.x) by Bacaro, Ricotta, &amp; Mazzoleni (2007).
Originally proposed for plant biology purposes, we re-purpose the formulae for use with microbial community data, specifically amplicon sequencing data.

The approach is simple:
- Construct a pairwise score between samples
- Calculate the Beta Diversity using the pairwise score.

Bacaro _et al_ propose two taxonomic distance measures. One is based on tree parsimony, the other based on differences in Linnean classification of taxa. For simplicity's sake (namely to avoid tree-building), we implement only the latter.

See the [Methodology doc](docs/Methodology.ipynb) for more info of the method.

## Requirements
Requires python >= 3.6.0. <br>
Packages listed in `requirements.txt.` <br>
- pandas
- numpy
- tabulate

## Usage
This is a small-scope software that expects you to have data from [Qiime2](https://qiime2.org/)'s [q2-feature-classifier](https://github.com/qiime2/q2-feature-classifier).
Each sample should have an associated `taxonomy.tsv` file. You may have to unzip Qiime artifacts to find these files.
We recommend that you give the `taxonomy.tsv` files meaningful names, e.g. `sample1_taxonomy.tsv`. 

This software accepts three arguments:
```
--input: A textfile where each line is the filepath for a taxonomy.tsv file
--l: L value used in Bacaro's beta. Specifies the depth of Linnean taxonomy to compare taxa. E.g. 7 is species, 6 is genus, etc. Recommend 7
--output: A directory that exists. Your output files will be written here.
```

To use Bacaros_beta:
- Download the code using git clone or otherwise
- Install dependencies, preferably in the virtual env of your choice (e.g. pip, conda)
- Create a directory for output
- Create a a textfile for input. Each line in the textfile should contain the filepath to a `taxonomy.tsv` file and nothing else.
- invoke the program:

`
python run_beta.py --i my_data.txt --l 7 --o my_results.txt
`

Your results will be printed to the console **and** saved in the directory you specified. 

##### Test
We have included a small test set to demonstrate the functionality of our software.<br>
To run it (from the project directory):
```bash
mkdir results
python run_beta.py --i data/test_input.txt --l 7 --o results
```
