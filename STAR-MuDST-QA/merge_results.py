#!/usr/bin/env python3
"""
Unified MuDST Analysis Results Merging Script
Author: Sneha Bhosale (AGH Krakow)
Usage: ./merge_results.py <tag> [--chunk-size-mb] 
"""

import os
import sys
import subprocess
from pathlib import Path
from glob import glob

def get_file_info(file_path):
    """Get file size in KB"""
    try:
        # Use ls -s to get size in KB (like original script)
        cmd = ["ls", "-s", str(file_path)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            parts = result.stdout.strip().split()
            if len(parts) >= 2:
                return int(parts[0]), str(file_path)  # size in KB, path
    except Exception:
        pass
    return 0, str(file_path)

def create_merge_macro():
    """Create the ROOT merging macro"""
    macro_content = '''
void MergeFiles(const string tmpnam = "files.tmp", const string outfile = "merged_output.root") {
    // Load to prevent error messages about missing dictionary class
    // gSystem->Load("libstar-upc.so");
    
    // Load files to the merger
    TFileMerger *merg = new TFileMerger();
    cout << "Loading files for merging..." << endl;
    
    ifstream in1(tmpnam.c_str());
    string line;
    Int_t nFiles = 0;
    
    while(getline(in1, line)) {
        if(line.length() > 0) {
            merg->AddFile(line.c_str(), kFALSE);
            nFiles++;
        }
    }
    
    in1.close();
    
    if(nFiles == 0) {
        cout << "No files to merge, exiting." << endl; 
        return;
    }
    
    cout << nFiles << " files loaded for merging." << endl;
    
    // Perform the merging
    merg->OutputFile(outfile.c_str());
    Bool_t stat = merg->Merge();
    
    if(stat) {
        cout << "Successfully merged " << nFiles << " files into " << outfile << endl;
    } else {
        cout << "Error in merging files!" << endl;
    }
    
    delete merg;
}
'''
    
    with open("MergeFiles.C", "w") as f:
        f.write(macro_content)

def merge_chunk(file_list, output_path, chunk_index):
    """Merge a chunk of files using ROOT"""
    temp_file = "files_temp.tmp"
    
    # Write file list to temporary file
    with open(temp_file, "w") as f:
        for file_path in file_list:
            f.write(f"{file_path}\n")
    
    # Format output filename
    output_file = output_path / f"StRP_production_{chunk_index:04d}.root"
    
    # Create merge command
    merge_cmd = [
        "root", "-l", "-b", "-q", 
        f'MergeFiles.C("{temp_file}","{output_file}")'
    ]
    
    print(f"Merging chunk {chunk_index} ({len(file_list)} files)...")
    
    try:
        result = subprocess.run(merge_cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            # Check if output file was created and has reasonable size
            if output_file.exists() and output_file.stat().st_size > 1000:
                file_size = output_file.stat().st_size
                print(f"✓ Created {output_file.name} ({file_size/1024/1024:.1f} MB)")
                return str(output_file)
            else:
                print(f"✗ Merge failed: output file too small or missing")
                return None
        else:
            print(f"✗ ROOT merge command failed:")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            return None
            
    except Exception as e:
        print(f"✗ Error during merge: {e}")
        return None
    
    finally:
        # Clean up temporary file
        if Path(temp_file).exists():
            Path(temp_file).unlink()

def merge_results(job_dir, chunk_size_kb):
    """Main merging function"""
    print(f"Scanning for ROOT files in {job_dir}...")
    
    # Find all ROOT files in rootfiles directory (not base directory)
    root_files = glob(str(job_dir / "rootfiles" / "*.root"))  # Changed this line
    
    if not root_files:
        print("No ROOT files found to merge!")
        return []
    
    print(f"Found {len(root_files)} ROOT files")
    
    # Get file sizes and filter out empty files
    file_info = []
    total_input_size = 0
    
    for root_file in root_files:
        size_kb, path = get_file_info(root_file)
        if size_kb > 0:  # Only include non-empty files
            file_info.append((size_kb, path))
            total_input_size += size_kb
    
    if not file_info:
        print("No non-empty ROOT files found!")
        return []
    
    print(f"Total input size: {total_input_size/1024:.1f} MB")
    print(f"Target chunk size: {chunk_size_kb/1024:.1f} MB")
    
    # Create merged directory
    merged_dir = job_dir / "merged"
    merged_dir.mkdir(exist_ok=True)
    
    # Create merging macro
    create_merge_macro()
    
    # Group files into chunks
    chunks = []
    current_chunk = []
    current_size = 0
    
    for size_kb, file_path in file_info:
        # If adding this file would exceed chunk size and we have files in current chunk
        if current_size + size_kb > chunk_size_kb and current_chunk:
            chunks.append((current_chunk.copy(), current_size))
            current_chunk = []
            current_size = 0
        
        current_chunk.append(file_path)
        current_size += size_kb
    
    # Add the last chunk if it has files
    if current_chunk:
        chunks.append((current_chunk, current_size))
    
    print(f"Creating {len(chunks)} merged chunks...")
    
    # Merge each chunk
    merged_files = []
    chunk_list_file = merged_dir / "StRP_production.list"
    
    with open(chunk_list_file, "w") as list_file:
        for i, (chunk_files, chunk_size) in enumerate(chunks):
            print(f"\nChunk {i}: {len(chunk_files)} files, {chunk_size/1024:.1f} MB")
            
            merged_file = merge_chunk(chunk_files, merged_dir, i)
            if merged_file:
                merged_files.append(merged_file)
                list_file.write(f"{merged_file}\n")
                
                # Show merged file info
                merged_path = Path(merged_file)
                if merged_path.exists():
                    size_mb = merged_path.stat().st_size / 1024 / 1024
                    print(f"  Output: {merged_path.name} ({size_mb:.1f} MB)")
    
    # Clean up merge macro
    if Path("MergeFiles.C").exists():
        Path("MergeFiles.C").unlink()
    
    return merged_files

def main():
    if len(sys.argv) < 2:
        print("Usage: ./merge_results.py <tag> [--chunk-size-mb SIZE]")
        print("Examples:")
        print("  ./merge_results.py test1")
        print("  ./merge_results.py test1 --chunk-size-mb 5.0")
        sys.exit(1)
    
    tag = sys.argv[1]
    
    # Parse chunk size
    chunk_size_mb = 2.5  # Default from original script
    if "--chunk-size-mb" in sys.argv:
        try:
            idx = sys.argv.index("--chunk-size-mb")
            if idx + 1 < len(sys.argv):
                chunk_size_mb = float(sys.argv[idx + 1])
        except (ValueError, IndexError):
            print("Invalid chunk size specified, using default 2.5 MB")
    
    chunk_size_kb = int(chunk_size_mb * 1024)  # Convert to KB
    
    # Configuration
    base_dir = Path("/path/to/your/output/directory")
    job_dir = base_dir / tag
    
    print("=" * 70)
    print(f"    Results Merging - Tag: {tag}")
    print("=" * 70)
    print(f"Job Directory: {job_dir}")
    print(f"Chunk Size: {chunk_size_mb} MB")
    print("=" * 70)
    
    if not job_dir.exists():
        print(f"Error: Job directory {job_dir} does not exist!")
        sys.exit(1)
    
    # Check if we have ROOT in the environment
    try:
        result = subprocess.run(["which", "root"], capture_output=True, text=True)
        if result.returncode != 0:
            print("Error: ROOT not found in PATH!")
            print("Please load the STAR environment: starver SL20c")
            sys.exit(1)
    except Exception:
        print("Error: Could not check for ROOT installation")
        sys.exit(1)
    
    # Perform merging
    merged_files = merge_results(job_dir, chunk_size_kb)
    
    print("=" * 70)
    
    if merged_files:
        print(f"✓ Successfully created {len(merged_files)} merged files:")
        merged_dir = job_dir / "merged"
        
        total_merged_size = 0
        for merged_file in merged_files:
            file_path = Path(merged_file)
            if file_path.exists():
                size_mb = file_path.stat().st_size / 1024 / 1024
                total_merged_size += size_mb
                print(f"  {file_path.name} ({size_mb:.1f} MB)")
        
        print(f"\nTotal merged size: {total_merged_size:.1f} MB")
        print(f"File list: {merged_dir}/StRP_production.list")
        
        print(f"\nMerged files are ready for analysis in:")
        print(f"  {merged_dir}")
        
    else:
        print("✗ No files were successfully merged!")
        print("Check the error messages above for details.")

if __name__ == "__main__":
    main()
