import ast
import pandas as pd
import glob
import time
import os

from tqdm import tqdm as tqdm
import subprocess
import json
import numpy as np

from optparse import OptionParser


def _wait_for_running_and_pending_lsf_cluster_jobs(
    maximum_number_of_jobs: int=0,
    sleep_timer: int=30
):
    """
    Wait for the number of running and pending jobs on an LSF cluster to drop below a specified maximum.
    Parameters
    ----------
    maximum_number_of_jobs : int, optional
        The maximum number of running and pending jobs to wait for. If this number is exceeded, the function will wait until the number drops below this value. Default is 0.
    sleep_timer : int, optional
        The number of seconds to wait between checks for the number of running and pending jobs. Default is 30.

    Returns
    -------
    None
    """
    assert type(maximum_number_of_jobs) == int
    assert type(sleep_timer) == int

    while True:
        # Get the number of running and pending jobs.
        num_running_and_pending_jobs = int(os.popen("bjobs -rp | wc -l").read())

        if num_running_and_pending_jobs <= maximum_number_of_jobs:
            break
        
        print(f"Waiting for all jobs to complete... Number of jobs: {num_running_and_pending_jobs}")
        time.sleep(sleep_timer)

    return

def _submit_bsub_job_and_get_job_id(
    bsub_pid: str="test",
    bsub_queue: str="short",
    bsub_allowed_time: str="0:10",
    bsub_memory_gb: str="1",
    command: str="ls",
):
    command = f'bsub -J {bsub_pid} -q {bsub_queue} -W {bsub_allowed_time} -M {bsub_memory_gb}GB -R rusage[mem={bsub_memory_gb}GB] "{command}"'
    job_id = int(subprocess.check_output(command, shell=True).decode("utf-8").strip().split()[1].replace("<", "").replace(">", ""))
    return job_id

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option(
        "--path-to-vcf",
        action="store",
        type="str",
        dest="vcf",
        help="Path to the .vcf file to run deep pileup on. \n",
    )
    
    parser.add_option(
        "--path-to-metadata",
        action="store",
        type="str",
        dest="metadata",
        help="Path to the metadata file. \n",
    )

    parser.add_option(
        "--output-path-to-repository",
        action="store",
        type="str",
        dest="output_path_to_repository",
        help="Output path to the repository \n",
    )
    
    parser.add_option(
        "--path-to-dp-single-position-script",
        action="store",
        type="str",
        dest="path_to_dp_single_position_script",
        help="Absolute to the _deep_pileup_single_position.py file. \n",
    )
    
    parser.add_option(
        "--dkfz-cluster",
        action="store",
        type="str",
        default="True",
        dest="dkfz_cluster",
        help="Boolean if you're within the DKFZ cluster that uses the LSF job system. \n",
    )
    
    parser.add_option(
        "--starting-index",
        action="store",
        type="str",
        default="1",
        dest="starting_index",
        help="Starting index of the .vcf file of what you want to run DP on. \n",
    )
    
    parser.add_option(
        "--ending-index",
        action="store",
        type="str",
        default="-1",
        dest="ending_index",
        help="Ending index of the .vcf file of what you want to run DP on. \n",
    )
    
    (options, args) = parser.parse_args()
    
    path_to_vcf = str(options.vcf)
    assert os.path.exists(path_to_vcf), ".vcf file does not exist."
    
    path_to_metadata = str(options.metadata)
    assert os.path.exists(path_to_metadata), "metadata file does not exist."
    
    output_path_to_repository = str(options.output_path_to_repository)
    
    if not os.path.exists(output_path_to_repository):
        os.mkdir(output_path_to_repository)
    
    run_in_dkfz_cluster = ast.literal_eval(str(options.dkfz_cluster))
    path_to_dp_single_position_script = str(options.path_to_dp_single_position_script)
    
    starting_index = int(options.starting_index)
    ending_index = int(options.ending_index)
    
    # Read in the dataframe.
    data = pd.read_csv(
        path_to_vcf,
        delimiter="\t"
    )
    
    # Specify which indices to run Deep Pileup on.
    data = data.iloc[starting_index:ending_index]
    
    # Delete repeats within the dataframe so Deep Pileup isn't ran on the same gene/position.
    data = data[["GENE", "#CHROM", "POS"]].drop_duplicates()
    data.reset_index(inplace=True, drop=True)
    
    # Iterate through the dataframe and run on each chromosome, gene and position.
    for _, row in tqdm(data.iterrows(), total=data.shape[0], desc="Running DP:"):
        chromosome = row["#CHROM"]
        gene = row["GENE"]
        position = row["POS"]
        
        command = f"python {path_to_dp_single_position_script} --path-to-metadata {path_to_metadata} --output-path-to-repository {output_path_to_repository} --gene {gene} --chromosome {chromosome} --position {position}"
        
        if run_in_dkfz_cluster:            
            _wait_for_running_and_pending_lsf_cluster_jobs(
                maximum_number_of_jobs = 100,
                sleep_timer = 10
            )
            
            job_id = _submit_bsub_job_and_get_job_id(
                bsub_pid = f"{gene}_chr{chromosome}_{position}",
                bsub_queue = "medium",
                bsub_allowed_time = "1:00",
                bsub_memory_gb = "10",
                command = command,
            )
        else:
            p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
            (output, err) = p.communicate()
            p_status = p.wait()
            
    if run_in_dkfz_cluster:
        _wait_for_running_and_pending_lsf_cluster_jobs(
            maximum_number_of_jobs = 0,
            sleep_timer = 10
        )
            
    
    
    
    
    
    