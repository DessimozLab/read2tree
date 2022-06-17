# docker run --rm -t -i -v $PWD/tests:/input -v $PWD/tests/:/reads -v $PWD/out:/app   read2tree_k   --tree --standalone_path /input/marker_genes --reads /reads/sample_1.fastq 

# base image
FROM continuumio/miniconda3

LABEL software="read2tree"

WORKDIR /app

# Create the environment:
COPY environment.yml .

RUN apt-get update && apt-get install -y wget && conda init bash

RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "read2tree_env", "/bin/bash", "-c"]

# Make sure the environment is activated:
RUN echo "Make sure numpy is installed:"

RUN python -c "import numpy"

RUN python -c "import ete3"

RUN python -c "import yaml"


RUN python -m pip install pysam
RUN python -m pip install pyham



ENV PATH /opt/conda/envs/read2tree_env/bin:$PATH


RUN git clone https://github.com/DessimozLab/read2tree.git
WORKDIR read2tree


RUN python setup.py install


ENV PATH  /app/read2tree/bin:$PATH


ENTRYPOINT ["read2tree"]

CMD ["-h"] 


