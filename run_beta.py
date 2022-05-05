from beta import beta
import argparse
import sys
import os

parser = argparse.ArgumentParser(prog="Bacaro's Beta",
                                usage="python run_beta.py -i infile -m ['s','t'] -l int -o outdir",
                                description="Calculate pairwise distances between classified OTU/ASV samples and beta diversity among samples."
)
parser.add_argument('--input', action='store', type=str, required=True)
parser.add_argument('--metric',
                    action='store',
                    type=str,
                    required=True
                    )
parser.add_argument('--l', action='store', type=int, required=True)
parser.add_argument('--output', action='store', type=str, required=True)

if __name__ == '__main__':
    args = parser.parse_args()
    with open(args.input, 'r') as f:
        files = [line.strip('\n') for line in f.readlines()]

    L = args.l
    metric = args.metric.lower()
    try:
        assert metric in ['s', 't']
    except:
        raise ValueError("metric must be one of 's' or 't'")

    samples = [{'name': ''.join(os.path.basename(file).split('.')[:-1]),
                'taxa': beta.qiime2_tsv_to_taxa_list(file)} for file in files]
    deltas, b = beta.calculate_beta(samples, L, metric)
    print("Beta Diversity for {0} samples: {1}".format(len(samples), b))
    #print("Distance matrix:")
    #print(deltas.to_markdown())

    outfile_csv = '{0}/{1}.csv'.format(args.output, ''.join(os.path.basename(args.input).split('.')[:-1]))
    outfile_txt = '{0}/{1}.txt'.format(args.output, ''.join(os.path.basename(args.input).split('.')[:-1]))
    deltas.to_csv(outfile_csv)
    with open(outfile_txt, 'w') as f:
        f.write("The calculated beta diversity for {0} was: {1}".format(args.input, b))
