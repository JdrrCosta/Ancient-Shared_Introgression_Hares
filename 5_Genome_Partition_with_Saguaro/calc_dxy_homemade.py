import csv
import gzip
import argparse
from collections import defaultdict

def parse_geno(geno):
    return geno.split("/")

def count_differences(alleles1, alleles2):
    return sum(a1 != a2 for a1, a2 in zip(alleles1, alleles2)) / len(alleles1)

def load_pop_file(pop_file):
    pop_dict = defaultdict(list)
    with open(pop_file) as pf:
        for line in pf:
            sample, pop = line.strip().split()
            pop_dict[pop].append(sample)
    return pop_dict

def main():
    parser = argparse.ArgumentParser(description="Calculate overall Dxy between two populations.")
    parser.add_argument("--geno", required=True, help="Input gzipped .geno file")
    parser.add_argument("--popfile", required=True, help="Population file (2 columns: sample pop)")
    parser.add_argument("--pop1", required=True, help="First population name")
    parser.add_argument("--pop2", required=True, help="Second population name")
    args = parser.parse_args()

    pop_dict = load_pop_file(args.popfile)

    if args.pop1 not in pop_dict or args.pop2 not in pop_dict:
        raise ValueError(f"One or both population names ({args.pop1}, {args.pop2}) not found in {args.popfile}")

    pop1_samples = pop_dict[args.pop1]
    pop2_samples = pop_dict[args.pop2]

    total_diffs = 0
    total_comparisons = 0

    with gzip.open(args.geno, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        sample_names = header[2:]
        sample_idx = {name: i + 2 for i, name in enumerate(sample_names)}

        # Sanity check: are all samples in the file?
        missing = [s for s in pop1_samples + pop2_samples if s not in sample_idx]
        if missing:
            raise ValueError(f"Samples not found in file header: {missing}")

        for row in reader:
            for s1 in pop1_samples:
                for s2 in pop2_samples:
                    g1 = row[sample_idx[s1]]
                    g2 = row[sample_idx[s2]]
                    if g1 != "N/N" and g2 != "N/N":
                        a1 = parse_geno(g1)
                        a2 = parse_geno(g2)
                        if len(a1) == len(a2):
                            total_diffs += count_differences(a1, a2)
                            total_comparisons += 1

    if total_comparisons > 0:
        dxy = total_diffs / total_comparisons
        print(f"Overall Dxy between {args.pop1} and {args.pop2}: {dxy:.6f}")
    else:
        print("No valid comparisons found.")

if __name__ == "__main__":
    main()

