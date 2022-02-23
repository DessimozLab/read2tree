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
        self.is_longread = '/ONT/' in self.read_files[0] 


class JobProducer(object):
    def __init__(self, base_path, species_list, submit):
        self.base_path = base_path
        self.submit = submit
        if isinstance(species_list, (str, bytes)):
            self.species = [Species(species_list)]
        else:
            self.species = [Species(sp) for sp in species_list]
        self.references = [z.split('_')[0] for z in os.listdir(os.path.join(base_path, "02_ref_dna"))]
        self.nr_references = len(self.references)
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
        if self.submit:
            mjob_ids = self.submit_jobs(mjobs)
        else:
            mjob_ids = []
        all_jobs.extend(mjob_ids)
        merge_job = self.merge_results(mjob_ids)
        if self.submit:
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
        read_type = '' if not species.is_longread else '--read_type long'
        return f"""#!/bin/bash
#SBATCH --output={self.job_log}/r2t_map_{reference}-{species.name}.out
#SBATCH --account=cdessim2_read2tree
#SBATCH --time={self.estimated_runtime(species)}
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --job-name r2t_map_{reference}-{species.name}
#SBATCH --mem={self.estimated_memory(species)}
#SBATCH --export=None

#source /scratch/wally/FAC/FBM/DBC/cdessim2/default/aaltenho/miniconda3/etc/profile.d/conda.sh
source ~/.conda-init
conda activate r2t

cd {os.path.abspath(self.base_path)}

read2tree {read_type} --reads {" ".join(species.read_files)} --output_path ./ --single_mapping {reference} --threads 4 --species_name {species.name}"""

    def merge_job_string(self, map_job_ids):
        condition = ""
        if len(map_job_ids) > 0:
            condition = f"#SBATCH --dependency afterok:{':'.join(map_job_ids)}"

        return f"""#!/bin/bash
#SBATCH --output={self.job_log}/r2t_merge.out
#SBATCH --account=cdessim2_read2tree
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --job-name r2t_merge
#SBATCH --mem=20GB
#SBATCH --export=None
{condition}

source ~/.conda-init
conda activate r2t

cd {os.path.abspath(self.base_path)}

read2tree --output_path ./ --merge_all_mappings
"""


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", nargs="+", help="species (reads). If several read files, just provide the "
                                                     "common prefix of them. Basename is species name")
    parser.add_argument("--only-jobfiles", default=False, action="store_true",
                        help="do not submit slurm jobs, but rather just create the job files. Note that "
                             "the final merge step can only run after all other jobs finished. User "
                             "needs to make sure this is fulfilled if using this flag.")
    parser.add_argument("base_path", help="base path of analysis")
    conf = parser.parse_args()

    producer = JobProducer(conf.base_path, conf.species, submit=not conf.only_jobfiles)
    producer.produce_jobs()
