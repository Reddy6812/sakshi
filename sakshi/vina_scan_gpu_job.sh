#!/bin/bash

#SBATCH --job-name=vina_scan_job        # Job name
#SBATCH --output=results_vina_scan.out  # Output file
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<vi.gade@ufl.edu>
#SBATCH --nodes=1
#SBATCH --error=results_vina_scan.err   # Error file
#SBATCH --time=14:00:00                 # Time limit (14 hours)
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=100              # Number of CPU cores per task (all available CPUs)
#SBATCH --mem=15GB                      # Total memory (15 GB, adjusted to the total available)
#SBATCH --gpus=a100:1 
#SBATCH --partition=gpu                 # Partition to submit to (using GPU partition)
#SBATCH --gres=gpu:a100:1                    # Request 1 GPU (adjust if more are needed)
#SBATCH --account=yanjun.li             # Account to charge the resources to
#SBATCH --mail-type=ALL                 # Email notifications (BEGIN, END, FAIL)
#SBATCH --mail-user=your-email@domain.com  # Replace with your email address

# Change to your working directory
cd /blue/yanjun.li/vi.gade1/seabra-li

# Load necessary modules
module load conda

# Activate the 'gscreen' conda environment
conda activate gscreen

# Run the gscreen command to start vina_scan.py
gscreen -s 1 -p vina


