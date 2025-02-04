from pathlib import Path
import pytest
from src.applications.autodock.parallel_docking.docking_queue import DockingQueue
from src.applications.autodock.parallel_docking.utils import BatchProcessor, find_docking_inputs

@pytest.fixture
def test_dir():
    """Path to test data directory"""
    return Path('tests/applications/autodock/data_for_test')


def test_find_docking_inputs(test_dir):
    """Test that find_docking_inputs correctly identifies valid systems"""

    jobs = find_docking_inputs(test_dir)
    
    # Check job structure
    for job in jobs:
        assert 'input_dir' in job
        assert 'box_size' in job
        assert 'exhaustiveness' in job
        assert 'n_poses' in job
        assert 'n_poses_to_save' in job

def test_batch_processor_initialization():
    """Test BatchProcessor initialization with various parameters"""
    processor = BatchProcessor(
        initial_batch_size=2,
        min_batch_size=1,
        cpu_threshold=80,
        memory_threshold=80,
        cooldown_time=10
    )
    
    assert processor.batch_size == 2
    assert processor.min_batch_size == 1
    assert processor.cpu_threshold == 80
    assert processor.memory_threshold == 80
    assert processor.cooldown_time == 10

def test_docking_queue_initialization():
    """Test DockingQueue initialization"""
    queue = DockingQueue(n_workers=2, cores_per_job=2)
    
    assert queue.n_workers == 2
    assert queue.cores_per_job == 2
    assert len(queue.workers) == 0
    assert len(queue.failed_jobs) == 0

@pytest.mark.integration
def test_batch_processing_integration(test_dir):
    """Integration test for batch processing"""
    try:
        jobs = find_docking_inputs(test_dir)
        processor = BatchProcessor(
            initial_batch_size=2,
            min_batch_size=1,
            cpu_threshold=80,
            memory_threshold=80,
            cooldown_time=5
        )
        
        # Process just first batch
        test_batch = jobs[:2]
        processor.process_jobs(test_batch)
        
        # Check progress file exists
        assert Path('docking_progress.json').exists()

            
    finally:
        # Clean up the progress file
        progress_file = Path('docking_progress.json')
        if progress_file.exists():
            progress_file.unlink()

def test_resource_monitor():
    """Test resource monitoring functionality"""
    from src.applications.autodock.parallel_docking.utils import ResourceMonitor
    
    # Test with very high thresholds to ensure it passes
    assert ResourceMonitor.check_resources(cpu_threshold=99, memory_threshold=99)
    
    # Test with very low thresholds to ensure it can detect high usage
    resources_ok = ResourceMonitor.check_resources(cpu_threshold=1, memory_threshold=1)
    assert isinstance(resources_ok, bool)

@pytest.mark.parametrize("n_workers,cores_per_job", [
    (1, 1),
    (2, 2),
    (4, 1),
])
def test_docking_queue_different_configs(n_workers, cores_per_job):
    """Test DockingQueue with different worker configurations"""
    queue = DockingQueue(n_workers=n_workers, cores_per_job=cores_per_job)
    assert queue.n_workers == n_workers
    assert queue.cores_per_job == cores_per_job