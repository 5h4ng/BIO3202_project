# PeCorA (Peptide Correlation Analysis) - MATLAB Implementation

## BIO3202 Course Project

This project is a MATLAB implementation of the PeCorA algorithm as part of the BIO3202 course. The goal is to reproduce the peptide correlation analysis method originally developed in R for detecting discordant peptide quantities in shotgun proteomics data.

## Background

PeCorA is a statistical method designed to identify peptides that show abundance changes inconsistent with their parent protein's overall expression pattern. This is particularly useful for detecting:

- Post-translational modifications (PTMs)
- Protein isoform-specific changes
- Differential peptide regulation
- Proteoform-level differences

The original work was published by Dermit et al. (2020) and implemented as an R package. This project reproduces the key functionality in MATLAB.

## Usage

### Method 1: Automated Execution (Recommended)

```bash
# Run complete workflow with default settings
./run_pecora_workflow.sh

# Specify custom MATLAB path
./run_pecora_workflow.sh /Applications/MATLAB_R2023b.app/bin/matlab

# Use custom input data
./run_pecora_workflow.sh matlab data/your_data.csv
```

### Method 2: Manual Step-by-Step

```matlab
% 1. Main analysis (preprocessing + statistics)
run_PeCorA

% 2. Individual peptide visualizations
run('plot_PeCorA_results')

% 3. Workflow demonstration figures
run('generate_figure2_workflow')

% 4. Comprehensive statistical summary
run('generate_summary_statistics')
```

## References

- **Original Publication:** Dermit, M., et al. (2020). "Peptide correlation analysis (PeCorA) reveals differential proteoform regulation." *Molecular & Cellular Proteomics*, 19(1), 135-146.
- **Original R Package:** https://github.com/jessegmeyerlab/PeCorA
- **Course:** BIO3202

## Attribution

This MATLAB implementation was developed as part of BIO3202 coursework. The original algorithm and statistical methods are credited to Dermit et al. (2020).

---

**To run the complete analysis:** Execute `./run_pecora_workflow.sh` in the project directory. 