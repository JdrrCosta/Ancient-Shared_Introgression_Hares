#make concatenated aligments 
python conacatenate.py

#concatenate protein coding
cat *_allgenes.fas > hares_mt_allgenes_final.fasta

#concatenate rrna genes
cat *_rrn.fas > hares_rrn.fasta

#align protein coding
mafft hares_mt_allgenes_final.fasta > hares_mt_allgenes_final_aligned.fasta
#align rrna
mafft hares_rrn.fasta > hares_rrn_aligned.fasta
