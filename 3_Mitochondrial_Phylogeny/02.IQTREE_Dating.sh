#make ML tree and least-squares dating with ancient samples
#input: hares_mt_allgenes_final_aligned_withAncient.fasta (76 contemporary + 3 ancient, codon-aware alignment)
#-keep-ident is needed or IQ-TREE silently drops identical sequences from the tree

./iqtree3 -s hares_mt_allgenes_final_aligned_withAncient.fasta -m GTR+G -pre ml_tree -nt AUTO -keep-ident

#build the LSD2 date file from the alignment headers
python make_dates_file.py hares_mt_allgenes_final_aligned_withAncient.fasta > dates.txt

#least-squares dating (LSD2 via IQ-TREE)
#node calibration: snowshoe-brown MRCA, 3.99 Mya (Ferreira et al. 2021)
#tip dates: 3 ancient samples, -2677 (775-541 BC, Best et al. 2022, converted to years before present)
#see dates.txt for full calibration

./iqtree3 -s hares_mt_allgenes_final_aligned_withAncient.fasta -te ml_tree.treefile --date dates.txt --date-ci 100 -o "Lame1,Lame2,Lame3,Lame4_pruned,Lame5,Lame6,Lame7,Lame8" -pre dated_output -keep-ident

#dated tree: dated_output.timetree.nex (open in FigTree, Node Labels > Display > date)
#rate + calibrated node age + CI: dated_output.timetree.lsd
