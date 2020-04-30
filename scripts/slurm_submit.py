import math
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
        if len(self.read_files) == 0:
            raise ValueError(f"no read files found for {self.name}")
        self.reads_size = sum(os.stat(f).st_size for f in self.read_files)


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

    def estimated_runtime(self, species):
        if species.reads_size < 200e6:
            return "00:20:00"
        elif species.reads_size < 2e9:
            return "01:00:00"
        elif species.reads_size < 8e9:
            return "04:00:00"
        else:
            return "12:00:00"

    def estimated_memory(self, species):
        gigs = max(5, math.ceil(species.reads_size/(2**30)))
        return f"{gigs}GB"

    def submit_jobs(self, jobs):
        job_ids = []
        for job in jobs:
            res = subprocess.run(["sbatch", job], cwd=self.base_path, stdout=subprocess.PIPE)
            res.check_returncode()
            job_id = res.stdout.split()[-1].decode()
            job_ids.append(job_id)
        return job_ids

    def produce_jobs(self):
        all_jobs = []
        for d in (self.job_dir, self.job_log):
            if not os.path.isabs(d):
                d = os.path.join(self.base_path, d)
            os.makedirs(d, exist_ok=True)
        mjobs = self.mapping_jobs()
        mjob_ids = self.submit_jobs(mjobs)
        all_jobs.extend(mjob_ids)
        merge_job = self.merge_results(mjob_ids)
        all_jobs.append(self.submit_jobs(merge_job))
        return all_jobs

    def mapping_jobs(self):
        jobs = []
        for species in self.species:
            for ref in self.references:
                with open(os.path.join(self.base_path, self.job_dir, f"r2t_map_{ref}-{species.name}"), 'wt') as jobfile:
                    jobfile.write(self.map_job_string(species, ref))
                    jobs.append(jobfile.name)
        return jobs

    def merge_results(self, map_job_ids):
        with open(os.path.join(self.base_path, self.job_dir, f"r2t_merge"), 'wt') as jobfile:
            jobfile.write(self.merge_job_string(map_job_ids))
        return jobfile.name

    def map_job_string(self, species, reference):
        return f"""#!/bin/bash
#SBATCH --output={self.job_log}/r2t_map_{reference}-{species.name}.out
#SBATCH --partition={os.getenv('CLUSTER')}
#SBATCH --account=cdessim2_default
#SBATCH --time={self.estimated_runtime(species)}
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --job-name r2t_map_{reference}-{species.name}
#SBATCH --mem={self.estimated_memory(species)}
#SBATCH --export=None

source /scratch/wally/FAC/FBM/DBC/cdessim2/default/aaltenho/miniconda3/etc/profile.d/conda.sh
conda activate r2t

cd {os.path.abspath(self.base_path)}

read2tree --reads {" ".join(species.read_files)} --output_path ./ --single_mapping {reference} --threads 4 --species_name {species.name}"""

    def merge_job_string(self, map_job_ids):
        condition = ""
        if len(map_job_ids) > 0:
            condition = f"#SBATCH --dependency afterok:{':'.join(map_job_ids)}"

        return f"""#!/bin/bash
#SBATCH --output={self.job_log}/r2t_merge.out
#SBATCH --partition={os.getenv('CLUSTER')}
#SBATCH --account=cdessim2_default
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --job-name r2t_merge
#SBATCH --mem=20GB
#SBATCH --export=None
{condition}

source /scratch/wally/FAC/FBM/DBC/cdessim2/default/aaltenho/miniconda3/etc/profile.d/conda.sh
conda activate r2t

cd {os.path.abspath(self.base_path)}

read2tree --output_path ./ --merge_all_mappings
"""


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", nargs="+", help="species (reads). If several read files, just provide the "
                                                     "common prefix of them. Basename is species name")
    parser.add_argument("base_path", help="base path of analysis")
    conf = parser.parse_args()

    producer = JobProducer(conf.base_path, conf.species)
    producer.produce_jobs()
