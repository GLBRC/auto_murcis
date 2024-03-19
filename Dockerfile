#  base image
FROM continuumio/miniconda3

# Set your working directory
WORKDIR /var/auto_murcis/

#RUN mkdir -p /var/auto_murcis/auto_murcis_file/
RUN git clone https://github.com/GLBRC/auto_murcis /var/auto_murcis/

# Create Conda environment for auto_murcis scripts
RUN conda env create -f /var/auto_murcis/auto_murcis_env.yaml

# Activate Conda environment
# reference: https://medium.com/@chadlagore/conda-environments-with-docker-82cdc9d25754
RUN echo "conda activate auto_murcis_env" >> ~/.bashrc
ENV PATH /opt/conda/envs/auto_murcis_env/bin:$PATH