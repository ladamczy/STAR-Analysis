# STAR-MuDST-QA: Unified Quality Analysis for Particle Embedding - User Documentation

Author : Sneha Bhosale (AGH Krakow)

## Overview

The Unified MuDST Analysis project is a comprehensive pipeline for processing and analyzing STAR experiment data, specifically designed for embedded particle physics simulations. This system provides tools for job submission, monitoring, and result merging on computing clusters.

## What This Project Does

This analysis framework processes MuDST (Micro Data Summary Tape) files from the STAR experiment at RHIC (Relativistic Heavy Ion Collider). It performs quality assurance analysis on embedded particle data including:

- **Pion analysis** (π+ and π- particles)
- **Kaon analysis** (K+ and K- particles)  
- **Proton analysis** (protons and antiprotons)

The analysis generates ROOT files containing histograms and data trees for physics studies of particle production and properties in high-energy collisions.

## System Requirements

### Computing Environment
- **STAR Software Framework**: SL20c environment
- **ROOT**: Part of STAR framework (root4star)
- **HTCondor**: Job submission system
- **Python 3**: For management scripts

### Required Files
- `UnifiedMuDstQA.C` - Main analysis code (ROOT macro)
- `unified_submit.xml` - Job submission template
- Input file lists (e.g., `kaons_all.list`, `protons_all.list`)

## Project Structure

```
project/
├── submit_analysis.py      # Job submission script
├── monitor_jobs.py         # Job monitoring script  
├── merge_results.py        # Results merging script
├── unified_submit.xml      # HTCondor job template
├── UnifiedMuDstQA.C        # Analysis ROOT macro           
├── pions_all.list          # Input file 
├── kaons_all.list          # Input file 
└── protons_all.list        # Input file 
```

## Getting Started

### 1. Environment Setup

Load the STAR software environment:
```bash
starver SL20c
```

Verify ROOT is available:
```bash
which root4star
```

### 2. Prepare Input Files

Create or verify your input file lists. Each list should contain absolute paths to MuDST files:
```
/star/data105/embedding/pp500_production_2017/KplusKminus_509_20222902/P20ic.SL20c/2017/130/18130092/st_zerobias_adc_18130092_raw_1000017.MuDst.root
/star/data105/embedding/pp500_production_2017/KplusKminus_509_20222902/P20ic.SL20c/2017/130/18130028/st_zerobias_adc_18130028_raw_1000009.MuDst.root
...
```

### 3. Required Analysis Code

Ensure `UnifiedMuDstQA.C` is in your working directory. This ROOT macro should accept the following parameters:
- Input file path
- Analysis options
- Output ROOT file name

## How to Run Analysis

### Step 0: Test Single File (Optional)

Before submitting large batch jobs, you can test the analysis on a single MuDST file:

```bash
# Load STAR environment
starver SL20c

# Run analysis on one file
root4star -l -b -q 'UnifiedMuDstQA.C("path/to/single/file.MuDst.root",1,0,"test_output.root",true)'
```

**Parameters for UnifiedMuDstQA.C:**
- `"path/to/file.MuDst.root"` - Input MuDST file path
- `1` - Number of events to process (1 = all events)
- `0` - Analysis mode/options
- `"test_output.root"` - Output ROOT file name
- `true` - Enable debug/verbose output

**Example with a real file:**
```bash
root4star -l -b -q 'UnifiedMuDstQA.C("/star/data105/embedding/pp500_production_2017/KplusKminus_509_20222902/P20ic.SL20c/2017/130/18130092/st_zerobias_adc_18130092_raw_1000017.MuDst.root",1,0,"kaon_test.root",true)'
```

This will:
- Process a single MuDST file locally
- Generate output ROOT file for inspection
- Help verify analysis code works correctly
- Useful for debugging before batch submission

### Step 1: Submit Jobs

Submit analysis jobs using a unique tag to identify your run:

```bash
# Using default pion file list
./submit_analysis.py my_run_tag

# Using custom file list
./submit_analysis.py kaon_analysis kaons_all.list
```

**Parameters:**
- `tag`: Unique identifier for your analysis run
- `input_list`: Path to file list (optional, defaults to `pions_all.list`)

**Example:**
```bash
./submit_analysis.py kaon_run_2024 kaons_all.list
```

The script will:
- Validate input files and analysis code
- Create job directories under `/path/to/your/output/directory/[tag]/`
- Generate HTCondor submission files
- Submit jobs to the queue

### Step 2: Monitor Job Progress

Check the status of your submitted jobs:

```bash
./monitor_jobs.py my_run_tag
```

**Output includes:**
- Number of submitted, running, completed, and failed jobs
- Progress bar visualization  
- Total output file size
- List of failed jobs (if any)

**Automatic resubmission of failed jobs:**
```bash
./monitor_jobs.py my_run_tag --resubmit
```

### Step 3: Merge Results

Once all jobs are complete, merge the output ROOT files:

```bash
# Using default chunk size (2.5 MB)
./merge_results.py my_run_tag

# Using custom chunk size
./merge_results.py my_run_tag --chunk-size-mb 5.0
```

**Parameters:**
- `tag`: Same tag used for job submission
- `--chunk-size-mb`: Target size for merged files (default: 2.5 MB)

The merging process:
- Scans output ROOT files in `/path/to/your/output/directory/[tag]/rootfiles/`
- Groups files into appropriately-sized chunks
- Creates merged files: `StRP_production_0000.root`, `StRP_production_0001.root`, etc.
- Generates file list: `StRP_production.list`

## Directory Structure After Execution

After running analysis, your output directory will contain:

```
/path/to/your/output/directory/[tag]/
├── sched/              # Job scheduler files
├── outfiles/           # Job stdout logs  
├── errfiles/           # Job stderr logs
├── rootfiles/          # Individual analysis output files
├── merged/             # Merged ROOT files
│   ├── StRP_production_0000.root
│   ├── StRP_production_0001.root
│   └── StRP_production.list
└── submit_[tag].xml    # Job submission file
```

## Troubleshooting

### Common Issues

**1. "Template file unified_submit.xml not found"**
- Ensure `unified_submit.xml` is in your working directory

**2. "Analysis code UnifiedMuDstQA.C not found"**  
- Verify the ROOT macro is present and properly named

**3. "ROOT not found in PATH"**
- Load STAR environment: `starver SL20c`

**4. Jobs failing during execution**
- Check error logs in `errfiles/[jobid].err`
- Verify input MuDST files are accessible
- Check disk space in output directory

**5. No files to merge**
- Ensure jobs completed successfully
- Check `rootfiles/` directory for output files
- Verify ROOT files are not empty


## Advanced Usage

### Custom Analysis Parameters

Modify `UnifiedMuDstQA.C` to accept additional parameters or change analysis settings.

### Batch Processing Multiple Particle Types

```bash
# Submit multiple analyses
./submit_analysis.py pion_analysis pions_all.list
./submit_analysis.py kaon_analysis kaons_all.list  
./submit_analysis.py proton_analysis protons_all.list

# Monitor all runs
./monitor_jobs.py pion_analysis
./monitor_jobs.py kaon_analysis
./monitor_jobs.py proton_analysis
```

### Customizing Output Location

Edit the `base_output_dir` configuration in `submit_analysis.py`:
```python
config = {
    'base_output_dir': '/path/to/your/output/directory',
    'source_dir': os.getcwd(),
    'template_file': 'unified_submit.xml'
}
```

## Data Analysis Workflow

1. **Job Submission**: Process individual MuDST files in parallel
2. **Quality Monitoring**: Track job completion and identify issues
3. **Result Merging**: Combine output files for efficient analysis
4. **Physics Analysis**: Use merged ROOT files for final physics studies

## Performance Tips

- **Chunk Size**: Adjust merge chunk size based on available memory and analysis needs
- **File Lists**: Split large file lists for better job distribution
- **Resource Limits**: Monitor job wall time and adjust if needed
- **Storage**: Ensure sufficient disk space for output files

## Support and Troubleshooting

For technical issues:
1. Check job logs in `outfiles/` and `errfiles/` directories
2. Verify STAR environment is properly loaded
3. Ensure adequate disk space and file permissions
4. Contact system administrators for cluster-specific issues

## Example Complete Workflow

```bash
# 1. Setup environment
starver SL20c

# 2. Submit kaon analysis
./submit_analysis.py kaon_run_dec2024 kaons_all.list

# 3. Monitor progress (repeat as needed)
./monitor_jobs.py kaon_run_dec2024

# 4. Merge results when complete
./merge_results.py kaon_run_dec2024

# 5. Analysis complete - files ready in:
# /path/to/your/output/directory/kaon_run_dec2024/merged/
```

Your merged ROOT files are now ready for physics analysis!
