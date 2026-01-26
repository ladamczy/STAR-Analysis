#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path
from glob import glob
import re

def get_condor_jobs():
    """Get running HTCondor jobs"""
    running_jobs = {}
    try:
        cmd = ["condor_q", os.getlogin(), "-nobatch"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        for line in result.stdout.split('\n'):
            # Look for job IDs and file arguments
            if 'run_preselection.sh' in line:
                parts = line.split()
                if len(parts) > 0:
                    cluster_id = parts[0].split('.')[0]
                    # Extract file ID from arguments if visible
                    if '/18' in line:
                        match = re.search(r'/(\d{8})/\1\.root', line)
                        if match:
                            file_id = match.group(1)
                            running_jobs[file_id] = cluster_id
    except Exception as e:
        print(f"Warning: Could not query condor jobs: {e}")
    
    return running_jobs

def get_job_status(work_dir, input_list, output_dir):
    """Get comprehensive job status"""
    work_path = Path(work_dir)
    output_path = Path(output_dir)
    
    # Read input file list
    submitted_jobs = {}
    with open(input_list, 'r') as f:
        for line in f:
            line = line.strip()
            if line and line.endswith('.root'):
                # Extract file ID (e.g., 18080039 from path)
                match = re.search(r'/(\d{8})/\1\.root', line)
                if match:
                    file_id = match.group(1)
                    submitted_jobs[file_id] = line
    
    # Get running jobs
    running_jobs = get_condor_jobs()
    
    # Check completed jobs
    completed_jobs = {}
    failed_jobs = {}
    total_presel_size = 0
    total_hist_size = 0
    
    for file_id in submitted_jobs.keys():
        presel_file = output_path / f"presel_{file_id}.root"
        hist_file = output_path / f"hist_{file_id}.root"
        
        presel_exists = presel_file.exists()
        hist_exists = hist_file.exists()
        
        if presel_exists and hist_exists:
            presel_size = presel_file.stat().st_size
            hist_size = hist_file.stat().st_size
            
            # Check if files are not just empty ROOT files (>1KB)
            if presel_size > 1000 and hist_size > 1000:
                completed_jobs[file_id] = {
                    'presel_size': presel_size,
                    'hist_size': hist_size
                }
                total_presel_size += presel_size
                total_hist_size += hist_size
            else:
                # Empty or corrupted files
                failed_jobs[file_id] = 'empty_output'
        elif file_id not in running_jobs:
            # Not running and no output = failed
            failed_jobs[file_id] = 'no_output'
    
    return {
        'submitted': submitted_jobs,
        'running': running_jobs,
        'completed': completed_jobs,
        'failed': failed_jobs,
        'total_presel_size': total_presel_size,
        'total_hist_size': total_hist_size,
        'work_dir': work_path,
        'output_dir': output_path
    }

def check_logs(work_dir, file_id):
    """Check log files for errors"""
    logs_dir = Path(work_dir) / "logs"
    errors = []
    
    # Find log files matching this file ID
    for err_file in logs_dir.glob(f"Preselection.*.*.err"):
        with open(err_file, 'r') as f:
            content = f.read()
            if file_id in content:
                if 'error' in content.lower() or 'failed' in content.lower():
                    # Extract first error line
                    for line in content.split('\n'):
                        if 'error' in line.lower() or 'failed' in line.lower():
                            errors.append(line.strip()[:80])
                            break
    
    return errors

def resubmit_jobs(failed_jobs, input_list, work_dir):
    """Create resubmission script for failed jobs"""
    if not failed_jobs:
        print("No jobs to resubmit.")
        return
    
    # Read original input list
    all_files = {}
    with open(input_list, 'r') as f:
        for line in f:
            line = line.strip()
            if line and line.endswith('.root'):
                match = re.search(r'/(\d{8})/\1\.root', line)
                if match:
                    file_id = match.group(1)
                    all_files[file_id] = line
    
    # Create resubmission list
    resubmit_file = Path(work_dir) / "resubmit_failed.txt"
    resubmit_sub = Path(work_dir) / "preselection_resubmit.sub"
    
    with open(resubmit_file, 'w') as f:
        for file_id in failed_jobs.keys():
            if file_id in all_files:
                f.write(all_files[file_id] + '\n')
    
    # Create condor submit file for resubmission
    with open(resubmit_sub, 'w') as f:
        f.write("""executable     = run_preselection.sh
arguments      = $(infile)
output         = logs/Preselection.$(ClusterId).$(Process).out
error          = logs/Preselection.$(ClusterId).$(Process).err
log            = logs/Preselection.$(ClusterId).log

should_transfer_files = NO
request_cpus   = 1
request_memory = 4 GB

queue infile from resubmit_failed.txt
""")
    
    print(f"\nCreated resubmission files:")
    print(f"  File list: {resubmit_file}")
    print(f"  Submit file: {resubmit_sub}")
    print(f"\nTo resubmit {len(failed_jobs)} failed jobs, run:")
    print(f"  cd {work_dir}")
    print(f"  condor_submit preselection_resubmit.sub")

def format_size(size_bytes):
    """Format file size in human readable format"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"

def main():
    # Configuration - adjust these paths
    work_dir = "/home/sbhosale/Work/STAR-Analysis/Exclusive_Analysis"
    input_list = "New_star_data_full_part1.txt"
    output_dir = "DataAfterPreselection"
    
    should_resubmit = "--resubmit" in sys.argv
    show_details = "--details" in sys.argv
    
    print("=" * 70)
    print("    UPC DST Analysis Job Monitor")
    print("=" * 70)
    
    # Get job status
    status = get_job_status(work_dir, input_list, output_dir)
    
    # Print summary
    total = len(status['submitted'])
    running = len(status['running'])
    completed = len(status['completed'])
    failed = len(status['failed'])
    
    print(f"Input List:      {input_list}")
    print(f"Output Dir:      {status['output_dir']}")
    print("-" * 70)
    print(f"Total Jobs:      {total:>6}")
    print(f"Running:         {running:>6}")
    print(f"Completed:       {completed:>6}")
    print(f"Failed/Missing:  {failed:>6}")
    print("-" * 70)
    
    if completed > 0:
        print(f"Presel Size:     {format_size(status['total_presel_size'])}")
        print(f"Hist Size:       {format_size(status['total_hist_size'])}")
        completion_rate = completed / total * 100
        print(f"Completion:      {completion_rate:.1f}%")
    
    # Progress bar
    if total > 0:
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
    
    # Show running jobs
    if running > 0:
        print(f"\nCurrently Running ({running}):")
        for file_id, cluster_id in list(status['running'].items())[:5]:
            print(f"  {file_id} (Cluster {cluster_id})")
        if running > 5:
            print(f"  ... and {running - 5} more")
    
    # Show failed jobs
    if failed > 0:
        print(f"\nFailed/Missing Jobs ({failed}):")
        for file_id, reason in list(status['failed'].items())[:10]:
            print(f"  {file_id} - {reason}")
            if show_details:
                errors = check_logs(work_dir, file_id)
                for error in errors:
                    print(f"    Error: {error}")
        if failed > 10:
            print(f"  ... and {failed - 10} more")
    
    print("=" * 70)
    
    # Handle resubmission
    if should_resubmit and failed > 0:
        print("\nPreparing resubmission...")
        resubmit_jobs(status['failed'], input_list, work_dir)
    elif failed > 0 and not should_resubmit:
        print(f"\nTo prepare resubmission files, run:")
        print(f"  ./monitor_upcdst_jobs.py --resubmit")
    
    # Recommendations
    print("\nNext Steps:")
    if completed == total:
        print("  ✓ All jobs completed! Ready for next analysis step")
    elif running > 0:
        print("  ⏳ Jobs still running. Check again later:")
        print(f"    ./monitor_upcdst_jobs.py")
    elif failed > 0:
        print("  ⚠️  Some jobs failed. Check logs for details:")
        print(f"    ./monitor_upcdst_jobs.py --details --resubmit")

if __name__ == "__main__":
    main()
