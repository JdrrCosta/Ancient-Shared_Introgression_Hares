#Merge genes in both cactus
cat genes_in_cactus4.txt genes_in_cactus12_2.txt \
| sort -u | grep -v "^LOC" > genes_in_both_cactus_no_loc.txt

#after creating the dataframe add "SYMBOL:" before all gene names or else it wont recognize the gene name
awk -F " " '{print "SYMBOL:" $1,$2}' gene_scores_matrix.txt \
> gene_scores_matrix_4_SYMBOL.txt 
