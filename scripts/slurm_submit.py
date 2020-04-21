import os
import re
import sys
import glob
import time
import subprocess
import typing as t


class Species(object):
    def __init__(self, species: t.Union[os.PathLike, str, bytes]):
        self.name = os.path.basename(species)
        self.read_files = glob.glob(species+"*")


class JobProducer(object):
    def __init__(self, base_path, species_list):
        self.base_path = base_path
        if isinstance(species_list, (str, bytes)):
            self.species = [Species(species_list)]
        else:
            self.species = [Species(sp) for sp in species_list]
        self.nr_references = len(os.listdir(os.path.join(base_path, "02_ref_dna")))
        self.references = ["RG{:03d}".format(x+1) for x in range(self.nr_references)]
        self.job_dir = "jobs/"
        self.job_log = "logs/"

    def set_job_dir(self, job_dir):
        self.job_dir = job_dir

    def set_job_log_dir(self, log_dir):
        self.job_log = log_dir

    def produce_jobs(self):
        for d in (self.job_dir, self.job_log):
            if not os.path.isabs(d):
                d = os.path.join(self.base_path, d)
            os.makedirs(d, exist_ok=True)
        mjobs = self.mapping_jobs()

    def mapping_jobs(self):
        jobs = []
        for species in self.species:
            for ref in self.references:
                with open(os.path.join(self.base_path, self.job_dir, f"r2t_{ref}-{species.name}"), 'wt') as jobfile:
                    jobfile.write(self.get_r2t_string(species, ref))
                    jobs.append(jobfile.name)
        return jobs

    def get_r2t_string(self, species, reference):
        return f"""#!/bin/bash
#SBATCH --output={self.job_log}/r2t_{species.name}-{reference}.out
#SBATCH --partition={os.getenv('CLUSTER')}
#SBATCH --account=cdessim2_default
#SBATCH --time=0:20:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --job-name r2t_{reference}-{species.name}
#SBATCH --mem=5GB
#SBATCH --export=None

source /scratch/wally/FAC/FBM/DBC/cdessim2/default/aaltenho/miniconda3/etc/profile.d/conda.sh
conda activate r2t

cd {os.path.abspath(self.base_path)}

read2tree --reads {" ".join(species.read_files)} --output_path ./ --single_mapping {reference} --threads 4 --species_name {species.name}"""


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", nargs="+", help="species (reads). If several read files, just provide the "
                                                     "common prefix of them. Basename is species name")
    parser.add_argument("base_path", help="base path of analysis")
    conf = parser.parse_args()

    producer = JobProducer(conf.base_path, conf.species)
    producer.produce_jobs()
