from pathlib import Path
import multiprocessing
import argparse
from src.applications.autodock.parallel_docking.utils import find_docking_inputs, BatchProcessor

def main():
    parser = argparse.ArgumentParser(description='Run parallel docking jobs')
    parser.add_argument('--input_dir', type=str, required=True, 
                      help='Head directory containing ligand/enzyme subdirectories')
    parser.add_argument('--batch_size', type=int, default=2,
                      help='Initial number of jobs per batch')
    parser.add_argument('--cpu_threshold', type=int, default=40,
                      help='CPU usage threshold percentage')
    parser.add_argument('--memory_threshold', type=int, default=40,
                      help='Memory usage threshold percentage')
    
    args = parser.parse_args()
    
    # Find all docking jobs
    print(f"Searching for docking jobs in: {args.input_dir}")
    jobs = find_docking_inputs(args.input_dir)
    print(f"Found {len(jobs)} docking jobs")
    
    if jobs:
        # Initialize batch processor
        processor = BatchProcessor(
            initial_batch_size=args.batch_size,
            min_batch_size=1,
            cpu_threshold=args.cpu_threshold,
            memory_threshold=args.memory_threshold,
            cooldown_time=60  # 1 minute cooldown
        )
        
        # Process all jobs
        processor.process_jobs(jobs)
    else:
        print("No valid docking jobs found")

if __name__ == '__main__':
    main()