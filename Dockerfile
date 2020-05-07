# Docker file for codon_usage_analysis
# Authors: Fanli Zhou

# Use continuumio/anaconda3 as base image
FROM continuumio/anaconda3 

# install python packages   
RUN conda install numpy 
RUN conda install scipy 
RUN conda install matplotlib

# Put Anaconda Python in PATH
ENV PATH="/opt/conda/bin:${PATH}"

CMD ["/bin/bash"]

