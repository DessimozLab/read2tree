# base image
FROM continuumio/miniconda3

LABEL software="read2tree"


WORKDIR /app

# Create the environment:
COPY environment.yml .

RUN apt-get -qq update \
    && apt-get install -y --no-install-recommends \
        wget \
    && rm -rf /var/lib/apt/lists/*

RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "read2tree_env", "/bin/bash", "-c"]

# Make sure the environment is activated:
RUN echo "Make sure numpy is installed:" \
    && python -c "import numpy" \
    && python -c "import ete3" \
    && python -c "import pysam"

COPY . .
RUN python setup.py install

ENV PATH /app/read2tree/bin:/opt/conda/envs/read2tree_env/bin:$PATH

WORKDIR /run

ENTRYPOINT ["read2tree"]

CMD ["-h"] 


