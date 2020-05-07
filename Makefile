# Driver script
# Fanli Zhou
# 
# This driver script completes the condon usage analysis. 
# This script takes no arguments.
# 
# usage: make all

# run all analysis
all: results/SP_top.png results/SP_bot.png results/LP_top.png results/LP_bot.png results/AAs_ChiSquare_bot.png results/Codons_ttest_bot.png results/AAs_ChiSquare_top.png results/Codons_ttest_top.png

# pre-processing/clean data
data/clean_data/SP_log_mean.txt data/clean_data/LP_log_mean.txt: src/read_data.py data/raw_data/microarray_data.csv
	python src/read_data.py "data/raw_data/microarray_data.csv" 

# find the top 10% genes for SP 
results/SP_top_gene.txt: src/find_gene.py data/clean_data/SP_log_mean.txt
	python src/find_gene.py "data/clean_data/SP_log_mean.txt"

# find the bottom 10% genes for SP 
results/SP_bot_gene.txt: src/find_gene.py data/clean_data/SP_log_mean.txt
	python src/find_gene.py --top 0 "data/clean_data/SP_log_mean.txt"

# find the top 10% genes for LP 
results/LP_top_gene.txt: src/find_gene.py data/clean_data/LP_log_mean.txt
	python src/find_gene.py "data/clean_data/LP_log_mean.txt"

# find the bottom 10% genes for LP 
results/LP_bot_gene.txt: src/find_gene.py data/clean_data/LP_log_mean.txt
	python src/find_gene.py --top 0 "data/clean_data/LP_log_mean.txt"

# data plot SP top genes
results/SP_top.png: src/data_plot.py data/clean_data/SP_log_mean.txt results/SP_top_gene.txt
	python src/data_plot.py "data/clean_data/SP_log_mean.txt" "results/SP_top_gene.txt"

# data plot SP bottom genes
results/SP_bot.png: src/data_plot.py data/clean_data/SP_log_mean.txt results/SP_bot_gene.txt
	python src/data_plot.py --label "bot" "data/clean_data/SP_log_mean.txt" "results/SP_bot_gene.txt"

# data plot LP top genes
results/LP_top.png: src/data_plot.py data/clean_data/LP_log_mean.txt results/LP_top_gene.txt
	python src/data_plot.py "data/clean_data/LP_log_mean.txt" "results/LP_top_gene.txt"

# data plot LP bottom genes
results/LP_bot.png: src/data_plot.py data/clean_data/LP_log_mean.txt results/LP_bot_gene.txt
	python src/data_plot.py --label "bot" "data/clean_data/LP_log_mean.txt" "results/LP_bot_gene.txt"

# Compare SP and LP top genes to remove intersections
results/SP_top_gene_dist.txt results/LP_top_gene_dist.txt: src/remove_intersections.py results/SP_top_gene.txt results/LP_top_gene.txt
	python src/remove_intersections.py --label "bot" "results/SP_top_gene.txt" "results/LP_top_gene.txt"

# Compare SP and LP bottom genes to remove intersections
results/SP_bot_gene_dist.txt results/LP_bot_gene_dist.txt: src/remove_intersections.py results/SP_bot_gene.txt results/LP_bot_gene.txt
	python src/remove_intersections.py --label "bot" "results/SP_bot_gene.txt" "results/LP_bot_gene.txt"

# Get SP and LP top gene coding sequences 
results/SP_top_seq.txt results/LP_top_seq.txt: src/get_gene_seq.py results/SP_top_gene_dist.txt results/LP_top_gene_dist.txt
	python src/get_gene_seq.py "results/SP_top_gene_dist.txt" "results/LP_top_gene_dist.txt"

# Get SP and LP bottom gene coding sequences 
results/SP_bot_seq.txt results/LP_bot_seq.txt: src/get_gene_seq.py results/SP_bot_gene_dist.txt results/LP_bot_gene_dist.txt
	python src/get_gene_seq.py --label "bot" "results/SP_bot_gene_dist.txt" "results/LP_bot_gene_dist.txt"

# Analyze codon usage in SP and LP top gene coding sequences 
results/AAs_ChiSquare_top.png results/Codons_ttest_top.png: src/codon_usage.py results/SP_top_seq.txt results/LP_top_seq.txt
	python src/codon_usage.py "results/SP_top_seq.txt" "results/LP_top_seq.txt"

# Analyze codon usage in SP and LP bottom gene coding sequences 
results/AAs_ChiSquare_bot.png results/Codons_ttest_bot.png: src/codon_usage.py results/SP_bot_seq.txt results/LP_bot_seq.txt
	python src/codon_usage.py --label "bot" "results/SP_bot_seq.txt" "results/LP_bot_seq.txt"

# Clean up intermediate and results files
clean: 
	rm -rf data/clean_data/*.txt
	rm -rf results/*.txt results/*.png
			