
from queue import Queue
import threading
from pathlib import Path
from src.applications.autodock.run_docking import run_vina_docking


class DockingQueue:
    def __init__(self, n_workers=4, cores_per_job=2):
        """
        Initialize docking queue
        
        Parameters:
        -----------
        n_workers : int
            Number of parallel workers
        cores_per_job : int
            Number of CPU cores to allocate per docking job
        """
        self.queue = Queue()
        self.n_workers = n_workers
        self.cores_per_job = cores_per_job
        self.workers = []
        self.failed_jobs = []
        
    def worker(self):
        """Worker function that processes jobs from the queue"""
        while True:
            job = self.queue.get()
            if job is None:  # Poison pill to stop worker
                break
                
            try:
                print(f"Starting job: {Path(job['input_dir']).name}")
                run_vina_docking(**job)
                print(f"Completed job: {Path(job['input_dir']).name}")
            except Exception as e:
                print(f"Failed job {Path(job['input_dir']).name}: {str(e)}")
                self.failed_jobs.append((job, str(e)))
            finally:
                self.queue.task_done()
                
    def submit_jobs(self, job_list):
        """
        Submit and process a list of docking jobs
        
        Parameters:
        -----------
        job_list : list
            List of dictionaries containing job parameters
        """
        print(f"Processing {len(job_list)} jobs with {self.n_workers} workers")
        
        # Start workers
        for _ in range(self.n_workers):
            t = threading.Thread(target=self.worker)
            t.daemon = True  # Allow program to exit if threads are still running
            t.start()
            self.workers.append(t)
            
        # Submit jobs to queue
        for job in job_list:
            self.queue.put(job)
            
        # Add poison pills to stop workers
        for _ in range(self.n_workers):
            self.queue.put(None)
            
        # Wait for all jobs to complete
        for t in self.workers:
            t.join()
            
        # Report any failures
        if self.failed_jobs:
            print("\nFailed jobs:")
            for job, error in self.failed_jobs:
                print(f"- {Path(job['input_dir']).name}: {error}")
