from beta import beta
import sys

if __name__ == '__main__':
    try:
        assert len(sys.argv) == 3
    except:
        print("The program takes exactly two argument; a list of tsv file paths and the L value (number of ranks to calculate against). Recommend L=7")
        exit(1)

    with open(sys.argv[1], 'r') as f:
        files = [line.strip('\n') for line in f.readlines()]

    L = int(sys.argv[2])
    samples = [beta.qiime2_tsv_to_taxa_list(file) for file in files]
    b = beta.calculate_beta(samples, L)
    print("Beta Diversity for {0} samples: {1}".format(len(samples), b))
