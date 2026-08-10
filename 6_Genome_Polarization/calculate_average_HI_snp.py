input_file = "Include_Lgra_all_chr_diem.vcf.gz.diem_snp_polarity.txt"
output_file = "ltim_averages.txt"

with open(input_file) as infile, open(output_file, "w") as outfile:
    infile.readline()  # skip header
    outfile.write("Chrom\tStart\tEnd\tLtim_avg\n")
    
    for line in infile:
        parts = line.strip().split()
        
        chrom, start, end = parts[0], parts[1], parts[2]
        polarity_string = parts[3][1:]  # remove leading "S"
        polarity_flag = int(parts[4])   # 0 or 1
        
        # Extract Ltim (positions 11â€“15 â†’ index 10:15)
        ltim = polarity_string[10:15]
        
        values = []
        for x in ltim:
            if x == "_":
                continue
            
            val = int(x)
            
            # Flip if polarity = 1
            if polarity_flag == 1:
                if val == 0:
                    val = 2
                elif val == 2:
                    val = 0
                # 1 stays 1
            
            values.append(val)
        
        avg = sum(values) / len(values) if values else "NA"
        
        outfile.write(f"{chrom}\t{start}\t{end}\t{avg}\n")
