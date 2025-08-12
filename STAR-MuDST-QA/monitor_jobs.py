#!/usr/bin/env python3
"""
Unified MuDST Analysis Job Monitoring Script
Author: Sneha Bhosale (AGH Krakow)
Usage: ./monitor_jobs.py <tag> [--resubmit]
"""

import os
import sys
import subprocess
from pathlib import Path
from glob import glob

def get_job_status(base_dir, tag):
    """Get comprehensive job status"""
    job_dir = Path(base_dir) / tag
    
    # Find submitted jobs from scheduler files
    sched_pattern = str(job_dir / "sched" / "*_*.csh")
    submitted_files = glob(sched_pattern)
    submitted_jobs = []
    for job_file in submitted_files:
        job_id = Path(job_file).stem.replace("sched", "")
        submitted_jobs.append(job_id)
    
    # Get running jobs from condor
    running_jobs = []
    try:
        cmd = ["condor_q", os.getlogin()]
        result = subprocess.run(cmd, capture_output=True, text=True)
        for line in result.stdout.split('\n'):
            if "sched/sched" in line and ".csh" in line:
                start = line.find("sched/sched") + len("sched/sched")
                end = line.find(".csh")
                if end > start:
                    job_id = line[start:end]
                    if job_id in submitted_jobs:
                        running_jobs.append(job_id)
    except Exception as e:
        print(f"Warning: Could not query running jobs: {e}")
    
    # Get completed jobs from output files
    completed_jobs = []
    total_size = 0
    output_pattern = str(job_dir / "rootfiles" / "*.root")
    output_files = glob(output_pattern)
    
    for output_file in output_files:
        try:
            size = os.path.getsize(output_file)
            if size > 0:  # Only count non-empty files
                job_id = Path(output_file).stem
                completed_jobs.append(job_id)
                total_size += size
        except OSError:
            continue
    
    # Find missing/failed jobs
    missing_jobs = []
    for job_id in submitted_jobs:
        if job_id not in running_jobs and job_id not in completed_jobs:
            missing_jobs.append(job_id)
    
    return {
        'submitted': submitted_jobs,
        'running': running_jobs,
        'completed': completed_jobs,
        'missing': missing_jobs,
        'total_size': total_size,
        'job_dir': job_dir
    }

def resubmit_jobs(missing_jobs, job_dir):
    """Resubmit failed jobs"""
    if not missing_jobs:
        print("No jobs to resubmit.")
        return
    
    print(f"Resubmitting {len(missing_jobs)} failed jobs...")
    
    sched_dir = job_dir / "sched"
    logs_dir = job_dir / "logs"
    
    original_dir = os.getcwd()
    resubmitted = 0
    
    try:
        os.chdir(sched_dir)
        
        for job_id in missing_jobs:
            # Clean up old log files
            for log_file in [f"{job_id}.out", f"{job_id}.err"]:
                log_path = logs_dir / log_file
                if log_path.exists():
                    log_path.unlink()
            
            # Parse job ID for resubmission
            parts = job_id.split("_")
            if len(parts) >= 2:
                session_file = f"sched{parts[0]}.session.xml"
                if Path(session_file).exists():
                    cmd = ["star-submit", "-r", parts[1], session_file]
                    try:
                        result = subprocess.run(cmd, capture_output=True, text=True)
                        if result.returncode == 0:
                            resubmitted += 1
                            print(f"✓ Resubmitted job {job_id}")
                        else:
                            print(f"✗ Failed to resubmit job {job_id}: {result.stderr}")
                    except Exception as e:
                        print(f"✗ Error resubmitting job {job_id}: {e}")
                else:
                    print(f"✗ Session file {session_file} not found for job {job_id}")
    
    finally:
        os.chdir(original_dir)
    
    print(f"Resubmitted {resubmitted}/{len(missing_jobs)} jobs successfully.")

def format_size(size_bytes):
    """Format file size in human readable format"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"

def main():
    if len(sys.argv) < 2:
        print("Usage: ./monitor_jobs.py <tag> [--resubmit]")
        print("Examples:")
        print("  ./monitor_jobs.py test1")
        print("  ./monitor_jobs.py test1 --resubmit")
        sys.exit(1)
    
    tag = sys.argv[1]
    should_resubmit = "--resubmit" in sys.argv
    
    # Configuration
    base_dir = "/path/to/your/output/directory"
    
    print("=" * 70)
    print(f"    Job Status Monitor - Tag: {tag}")
    print("=" * 70)
    
    # Get job status
    status = get_job_status(base_dir, tag)
    
    # Print summary
    print(f"Job Directory: {status['job_dir']}")
    print("-" * 70)
    print(f"Submitted Jobs:  {len(status['submitted']):>6}")
    print(f"Running Jobs:    {len(status['running']):>6}")
    print(f"Completed Jobs:  {len(status['completed']):>6}")
    print(f"Failed/Missing:  {len(status['missing']):>6}")
    print("-" * 70)
    
    if status['completed']:
        print(f"Total Output Size: {format_size(status['total_size'])}")
        completion_rate = len(status['completed']) / len(status['submitted']) * 100
        print(f"Completion Rate:   {completion_rate:.1f}%")
    
    # Show progress bar
    if status['submitted']:
        total = len(status['submitted'])
        completed = len(status['completed'])
        running = len(status['running'])
        failed = len(status['missing'])
        
        bar_width = 50
        completed_width = int(bar_width * completed / total)
        running_width = int(bar_width * running / total)
        failed_width = int(bar_width * failed / total)
        remaining_width = bar_width - completed_width - running_width - failed_width
        
        progress_bar = ("█" * completed_width + 
                       "▓" * running_width + 
                       "░" * failed_width + 
                       " " * remaining_width)
        
        print(f"\nProgress: [{progress_bar}]")
        print(f"Legend: ██ Completed  ▓▓ Running  ░░ Failed")
    
    # Show detailed status for failed jobs
    if status['missing']:
        print(f"\nFailed/Missing Jobs ({len(status['missing'])}):")
        for i, job_id in enumerate(status['missing'][:10]):  # Show first 10
            print(f"  {job_id}")
        if len(status['missing']) > 10:
            print(f"  ... and {len(status['missing']) - 10} more")
    
    print("=" * 70)
    
    # Handle resubmission
    if should_resubmit and status['missing']:
        print("\nResubmission requested...")
        resubmit_jobs(status['missing'], status['job_dir'])
    elif status['missing'] and not should_resubmit:
        print(f"\nTo resubmit failed jobs, run:")
        print(f"  ./monitor_jobs.py {tag} --resubmit")
    
    # Recommendations
    print("\nNext Steps:")
    if len(status['completed']) == len(status['submitted']):
        print("  ✓ All jobs completed! Ready to merge results:")
        print(f"    ./merge_results.py {tag}")
    elif status['running']:
        print("  ⏳ Jobs still running. Check again later:")
        print(f"    ./monitor_jobs.py {tag}")
    elif status['missing']:
        print("  ⚠️  Some jobs failed. Consider resubmitting:")
        print(f"    ./monitor_jobs.py {tag} --resubmit")

if __name__ == "__main__":
    main()
