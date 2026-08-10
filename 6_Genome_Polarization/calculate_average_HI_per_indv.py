import argparse

def main(input_file, output_file):

    individuals = [
        "Lcas1", "Lcas6", "Lcas7", "Lcas8", "Lcas9",
        "Lcor10", "Lcor6", "Lcor7", "Lcor8", "Lcor9",
        "Ltim1", "Ltim2", "Ltim3", "Ltim4", "Ltim6",
        "Lgra10", "Lgra1", "Lgra4", "Lgra5", "Lgra6",
        "Leur10", "Leur2", "Leur3", "Leur5", "Leur6"
    ]

    # Running sums and counts
    sums = {ind: 0 for ind in individuals}
    counts = {ind: 0 for ind in individuals}

    with open(input_file) as infile:
        infile.readline()  # skip header

        for line in infile:
            parts = line.strip().split()

            polarity_string = parts[3][1:]  # remove leading "S"
            polarity_flag = int(parts[4])

            for i, x in enumerate(polarity_string):

                if x == "_":
                    continue

                val = int(x)

                # Flip if polarity = 1
                if polarity_flag == 1:
                    if val == 0:
                        val = 2
                    elif val == 2:
                        val = 0

                sums[individuals[i]] += val
                counts[individuals[i]] += 1

    with open(output_file, "w") as outfile:
        outfile.write("Individual\tMean_HI\tNSNPs\n")

        for ind in individuals:
            if counts[ind] == 0:
                mean_hi = "NA"
            else:
                mean_hi = sums[ind] / counts[ind]

            outfile.write(f"{ind}\t{mean_hi}\t{counts[ind]}\n")

    print("Done!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute mean HI per individual"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file"
    )

    args = parser.parse_args()

    main(args.input, args.output)
