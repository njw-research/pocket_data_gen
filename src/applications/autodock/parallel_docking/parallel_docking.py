from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from src.applications.autodock.run_docking import run_vina_docking

def run_parallel_docking(base_dir: str, job_list: list, max_workers: int = None):
    """
    Run multiple docking jobs in parallel
    """
    if max_workers is None:
        cores_per_docking = 2
        total_cores = multiprocessing.cpu_count()
        max_workers = max(1, total_cores // cores_per_docking)
    
    print(f"Running with {max_workers} parallel workers")
    
    # Run jobs in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for job in job_list:
            future = executor.submit(run_vina_docking, **job)
            futures.append(future)
        
        # Wait for all jobs to complete
        for i, future in enumerate(futures):
            try:
                future.result()  # This will raise any exceptions that occurred
                print(f"Completed job {i}")
            except Exception as e:
                print(f"Job {i} failed with error: {str(e)}")

    print("\nCompleted all jobs")
    