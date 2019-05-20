
This project analyzes codon usage of two groups of genes selected based on microarray data

1. Read the microarray data to get SP and LP group data.

	python read_data.py microarray_data.csv

	Input: 'microarray_data.csv'

	Output: 'SP_log_mean.txt' and 'LP_log_mean.txt'


2. SP or LP expression data were sorted to find the top or bottom 10% (default) genes for each sample.

Only genes appear in all sample were saved. 

	python find_genes.py [-h] [--top TOP] [--cut_off CUT_OFF] data_file

	Input: 'SP_log_mean.txt' or 'LP_log_mean.txt'

	Output: 'SP_(label)_gene.txt' and 'LP_(label)_gene.txt'

Visualization of the results 

example output: .\Plotting_output_example\LP_data_plot_output_example.png

	python data_plot.py [-h] [--name NAME] [--label LABEL] data_files [data_files ...]

	Input: 'SP_log_mean.txt' and 'SP_(label)_gene.txt'

	or 'LP_log_mean.txt' and 'LP_(label)_gene.txt'


3. Compare SP and LP to remove intersections.


	python remove_intersections.py [-h] [--label LABEL] sp_file lp_file

	Input: 'SP_(label)_gene.txt' and 'LP_(label)_gene.txt'

	Output: 'SP_(label)_gene_dist.txt' and 'LP_(label)_gene_dist.txt'


4. Get gene coding sequences.

	python get_gene_seq.py [-h] [--label LABEL] sp_file lp_file

	Input: 'SP_(label)_gene_dist.txt' and 'LP_(label)_gene_dist.txt'

	Output: 'SP_(label)_seq.txt' and 'LP_(label)_seq.txt'


5. Analyze codon usage. 

	python codon_usage.py [-h] [--label LABEL] sp_file lp_file

	Input: 'SP_(label)_seq.txt' and 'LP_(label)_seq.txt'

Output:

	a. ChiSquare test p_value of each amino acid 

	example output: .\Plotting_output_example\AAs_ChiSquare_output_example.png
	
	b. Student's t test p_value of each codon 
		
	example output: .\Plotting_output_example\Codon_ttest_output_example.png
	
	c. Codon usage plot of each amino acid
	
	example output: .\Plotting_output_example\AA_plot_output_example.png


