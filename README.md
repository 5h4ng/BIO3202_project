# PeCorA (Peptide Correlation Analysis) - MATLAB Implementation

This is a MATLAB implementation of the PeCorA tool for detecting discordant peptide quantities in shotgun proteomics data.

## Directory Structure
```
PeCorA/
├── src/                    # Source code directory
│   ├── PeCorA_preprocessing.m
│   ├── PeCorA.m
│   └── PeCorA_plotting.m
├── data/                   # Data directory
│   └── PeCorA_noZ.csv     # Your input data file
├── results/               # Results output directory (created automatically)
└── run_PeCorA.m          # Main script
```

## Requirements
- MATLAB R2019b or newer
- Statistics and Machine Learning Toolbox

## Data Format
Your input CSV file should contain the following columns:
- `Peptide_Modified_Sequence`: Peptide sequence (including modifications)
- `Condition`: Experimental condition/group
- `BioReplicate`: Biological replicate number
- `Protein`: Protein identifier
- At least one column containing peak area data

## Usage

1. **Prepare Your Data**
   - Place your CSV data file in the `data` directory
   - Ensure all required columns are present
   - Data should not contain missing values

2. **Run the Analysis**
   ```matlab
   % Method 1: Run the main script directly
   run('run_PeCorA.m')
   
   % Method 2: In MATLAB command window
   >> cd /path/to/PeCorA
   >> run_PeCorA
   ```

3. **View Results**
   Results will be saved in the `results` directory:
   - `PeCorA_results.csv`: Contains all analysis results
   - `PeCorA_plot_1.png` to `PeCorA_plot_5.png`: Visualizations of top 5 significant peptides

## Customizing Parameters

You can modify parameters in `run_PeCorA.m`:

```matlab
% Modify preprocessing parameters
scaled_peptides = PeCorA_preprocessing(t, ...
    8, ...    % Change peak area column number
    100, ...  % Change filtering threshold
    'cntrl'); % Change control group name

% Modify number of visualizations
num_plots = min(5, height(disagree_peptides)); % Change number of plots to generate
```

## Output Description

1. **Preprocessed Data** (`scaled_peptides`):
   - Original data columns
   - `ms1log2`: Log2-transformed peak areas
   - `ms1scaled`: Normalized values
   - `ms1adj`: Values adjusted relative to control

2. **Analysis Results** (`disagree_peptides`):
   - Significantly different peptides (adj_pval ≤ 0.01)
   - Sorted by adjusted p-value

3. **Visualization**:
   - Boxplots showing peptide expression distribution
   - Gray: Other peptides
   - Green: Differential peptide
   - Individual data points included

## Troubleshooting

1. **Data File Issues**
   - Ensure the data file exists in the correct location
   - Check that column names match exactly
   - Verify data format is correct

2. **Memory Issues**
   - For large datasets, ensure sufficient system memory
   - Consider reducing data size or increasing system virtual memory

3. **MATLAB Errors**
   - Verify MATLAB version and required toolboxes
   - Check for any missing dependencies

## References
- Original R package: [PeCorA GitHub Repository](https://github.com/jessegmeyerlab/PeCorA)
- Documentation: [PeCorA Vignette](docs/PeCorA_vignette.pdf)
- Research Paper: [Dermit et al. 2020](docs/dermit-et-al-2020-peptide-correlation-analysis-(pecora)-reveals-differential-proteoform-regulation.pdf) 