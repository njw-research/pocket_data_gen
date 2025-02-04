from pathlib import Path
from typing import List, Dict
import multiprocessing
import psutil
import time
import json

from src.applications.autodock.parallel_docking.docking_queue import DockingQueue

class ResourceMonitor:
    @staticmethod
    def check_resources(cpu_threshold=90, memory_threshold=90):
        """Check if system resources are available"""
        cpu_percent = psutil.cpu_percent()
        memory_percent = psutil.virtual_memory().percent
        
        resources_ok = (cpu_percent < cpu_threshold and 
                       memory_percent < memory_threshold)
        
        print(f"System Status - CPU: {cpu_percent}%, Memory: {memory_percent}%")
        return resources_ok

class BatchProcessor:
    def __init__(self, 
                 initial_batch_size=4,
                 min_batch_size=1,
                 cpu_threshold=90,
                 memory_threshold=90,
                 cooldown_time=60):
        self.batch_size = initial_batch_size
        self.min_batch_size = min_batch_size
        self.cpu_threshold = cpu_threshold
        self.memory_threshold = memory_threshold
        self.cooldown_time = cooldown_time
        self.monitor = ResourceMonitor()
        
    def process_jobs(self, jobs: List[Dict]):
        """Process jobs in batches with resource monitoring"""
        remaining_jobs = jobs.copy()
        completed_batches = []
        batch_number = 1
        
        while remaining_jobs:
            print(f"\nProcessing Batch {batch_number}")
            print(f"Remaining jobs: {len(remaining_jobs)}")
            print(f"Current batch size: {self.batch_size}")
            
            # Check resources
            if not self.monitor.check_resources(self.cpu_threshold, 
                                              self.memory_threshold):
                print("System overloaded, reducing batch size and cooling down...")
                self.batch_size = max(self.min_batch_size, self.batch_size // 2)
                time.sleep(self.cooldown_time)
                continue
            
            # Process batch
            current_batch = remaining_jobs[:self.batch_size]
            remaining_jobs = remaining_jobs[self.batch_size:]
            
            try:
                # Set up queue processor for this batch
                total_cores = multiprocessing.cpu_count()
                cores_per_job = 2
                n_workers = min(len(current_batch), 
                              max(1, total_cores // cores_per_job))
                
                queue_processor = DockingQueue(n_workers=n_workers, 
                                            cores_per_job=cores_per_job)
                queue_processor.submit_jobs(current_batch)
                
                # Record successful batch
                completed_batches.append({
                    'batch_number': batch_number,
                    'size': len(current_batch),
                    'workers': n_workers,
                    'failed_jobs': queue_processor.failed_jobs
                })
                
                batch_number += 1
                
                # Optional: Save progress after each batch
                self._save_progress(completed_batches)
                
            except Exception as e:
                print(f"Batch {batch_number} failed: {str(e)}")
                # Return failed jobs to queue
                remaining_jobs = current_batch + remaining_jobs
                self.batch_size = max(self.min_batch_size, self.batch_size // 2)
                time.sleep(self.cooldown_time)
    
    def _save_progress(self, completed_batches):
        """Save progress to JSON file"""
        with open('docking_progress.json', 'w') as f:
            json.dump({
                'completed_batches': completed_batches,
                'timestamp': time.strftime("%Y%m%d-%H%M%S")
            }, f, indent=4)


def find_docking_inputs(head_dir: str) -> List[Dict]:
    """
    Find all docking systems from directory structure:
    head_dir/
        ligand_1/
            enzyme_1/
                ligand_aligned.pdbqt
                protein_altloc_removed_com.pdbqt
                active_site.txt
            enzyme_2/
                ligand_aligned.pdbqt
                protein_altloc_removed_com.pdbqt
                active_site.txt
        ligand_2/
            enzyme_1/
                ...
    """
    head_dir = Path(head_dir)
    jobs = []
    
    # Iterate through ligand directories
    for ligand_dir in head_dir.iterdir():
        if not ligand_dir.is_dir():
            continue
            
        # Iterate through enzyme directories under each ligand
        for enzyme_dir in ligand_dir.iterdir():
            if not enzyme_dir.is_dir():
                continue
                
            # Check for required files
            ligand_file = list(enzyme_dir.glob("*_aligned.pdbqt"))
            protein_file = list(enzyme_dir.glob("*_altloc_removed_com.pdbqt"))
            active_site = enzyme_dir / "active_site.txt"
            
            if ligand_file and protein_file and active_site.exists():
                job = {
                    'input_dir': str(enzyme_dir),
                    'box_size': [40, 40, 40],
                    'exhaustiveness': 32,
                    'n_poses': 10,
                    'n_poses_to_save': 10
                }
                jobs.append(job)
                print(f"Found complete system: {ligand_dir.name}/{enzyme_dir.name}")
    
    print(f"\nCreated {len(jobs)} docking jobs")
    return jobs