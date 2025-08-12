#!/usr/bin/env python3
"""
Unified MuDST Analysis Job Submission Helper
Author: Sneha Bhosale (AGH Krakow)
Usage: ./submit_analysis.py <tag> [input_list] 
"""

import os
import sys
import re
import subprocess
from pathlib import Path

def create_directories(base_dir):
    """Create necessary directories for job submission"""
    dirs = ['logs', 'sched', 'merged', 'outfiles', 'errfiles', 'rootfiles']
    for d in dirs:
        Path(base_dir / d).mkdir(parents=True, exist_ok=True)

def create_job_xml(template_file, output_file, replacements):
    """Create job XML from template with replacements"""
    with open(template_file, 'r') as f:
        content = f.read()
    
    for old, new in replacements.items():
        content = content.replace(old, new)
    
    with open(output_file, 'w') as f:
        f.write(content)

def submit_job(xml_file, sched_dir):
    """Submit job using star-submit"""
    original_dir = os.getcwd()
    try:
        os.chdir(sched_dir)
        cmd = ['star-submit', f'../{xml_file.name}']
        result = subprocess.run(cmd, capture_output=True, text=True)
        print("STDOUT:", result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        return result.returncode == 0
    finally:
        os.chdir(original_dir)

def main():
    if len(sys.argv) < 2:
        print("Usage: ./submit_analysis.py <tag> [input_list]")
        print("Examples:")
        print("  ./submit_analysis.py test1")
        print("  ./submit_analysis.py run2024 /path/to/my_files.list")
        sys.exit(1)
    
    # Parse arguments
    tag = sys.argv[1]
    input_list = sys.argv[2] if len(sys.argv) > 2 else "lists_pion/pions_all.list"
    
    # Convert input_list to absolute path
    input_list = os.path.abspath(input_list)
    
    # Configuration
    config = {
        'base_output_dir': '/path/to/your/output/directory',
        'source_dir': os.getcwd(),
        'template_file': 'unified_submit.xml'
    }
    
    # Setup paths
    job_output_dir = Path(config['base_output_dir']) / tag
    source_dir = Path(config['source_dir'])
    template_file = source_dir / config['template_file']
    
    print("=" * 60)
    print("    Unified MuDST Analysis Job Submission")
    print("=" * 60)
    print(f"Tag: {tag}")
    print(f"Input list: {input_list}")
    print(f"Output directory: {job_output_dir}")
    print(f"Source directory: {source_dir}")
    print("=" * 60)
    
    # Check if template exists
    if not template_file.exists():
        print(f"Error: Template file {template_file} not found!")
        print("Please ensure unified_submit.xml is in the current directory.")
        sys.exit(1)
    
    # Check if input list exists
    if not Path(input_list).exists():
        print(f"Error: Input file list {input_list} not found!")
        sys.exit(1)
    
    # Check if UnifiedMuDstQA.C exists
    analysis_code = source_dir / 'UnifiedMuDstQA.C'
    if not analysis_code.exists():
        print(f"Error: Analysis code {analysis_code} not found!")
        print("Please ensure UnifiedMuDstQA.C is in the current directory.")
        sys.exit(1)
    
    # Read and validate file list
    print("Checking file list...")
    with open(input_list, 'r') as f:
        files = [line.strip() for line in f if line.strip() and not line.startswith('#')]
    
    print(f"Found {len(files)} files to process")
    
    # Check if some files exist (sample check)
    missing_files = []
    for file_path in files[:3]:  # Check first 3 files
        if not Path(file_path).exists():
            missing_files.append(file_path)
    
    if missing_files:
        print("WARNING: Some files may not be accessible:")
        for f in missing_files:
            print(f"  {Path(f).name}")
        
        response = input("Continue anyway? (y/n): ")
        if response.lower() != 'y':
            sys.exit(1)
    
    # Create directories
    print("Creating job directories...")
    create_directories(job_output_dir)
    
    # Create job XML using traditional filelist method
    print("Creating job submission file...")
    input_source = f'<input URL="filelist:{input_list}"/>'
    replacements = {
        '__INPUT__': input_source,
        '__BASE_OUTPUT_DIR__': str(job_output_dir),
        '__SOURCE_DIR__': str(source_dir)
    }
    
    job_xml = job_output_dir / f'submit_{tag}.xml'
    create_job_xml(template_file, job_xml, replacements)
    
    # Debug: Show what was created
    print(f"Created job XML: {job_xml}")
    print("XML content preview:")
    with open(job_xml, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines[:15]):  # Show first 15 lines
            print(f"  {i+1:2d}: {line.rstrip()}")
    
    # Show files being processed
    print(f"\nFiles to process:")
    for i, f in enumerate(files[:3]):
        print(f"  {i+1}: {Path(f).name}")
    if len(files) > 3:
        print(f"  ... and {len(files)-3} more files")
    
    # Ask for confirmation  
    response = input(f"\nSubmit {len(files)} job(s)? (y/n): ")
    if response.lower() != 'y':
        print("Submission cancelled")
        sys.exit(0)
    
    # Submit job
    print("Submitting job...")
    sched_dir = job_output_dir / 'sched'
    
    success = submit_job(job_xml, sched_dir)
    
    if success:
        print("✓ Job submitted successfully!")
        print(f"Job files created in: {job_output_dir}")
        print(f"\nMonitor with:")
        print(f"  condor_q")
        print(f"  watch 'condor_q | tail -20'")
        print(f"\nCheck logs:")
        print(f"  tail -f {job_output_dir}/outfiles/*.out")
        print(f"  tail -f {job_output_dir}/errfiles/*.err")
        print(f"\nMerge results with: ./merge_results.py {tag}")
    else:
        print("✗ Job submission failed!")
        print("Check the error messages above.")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
