#!/bin/bash

# PeCorA Complete Workflow Script
# Usage: ./run_pecora_workflow.sh [matlab_path] [input_data]
# 
# This script runs the complete PeCorA (Peptide Correlation Analysis) workflow
# from raw data to final results and visualizations.

# Set default parameters
MATLAB_PATH=${1:-"matlab"}
INPUT_DATA=${2:-"data/PeCorA_noZ.csv"}

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if file exists
check_file() {
    if [ ! -f "$1" ]; then
        print_error "File not found: $1"
        exit 1
    fi
}

# Function to run MATLAB command
run_matlab() {
    local cmd="$1"
    local desc="$2"
    
    print_status "Running: $desc"
    echo "MATLAB Command: $cmd"
    
    $MATLAB_PATH -nodisplay -nosplash -nodesktop -r "try; $cmd; catch ME; fprintf('ERROR: %s\n', ME.message); exit(1); end; exit(0);" 2>&1
    
    if [ $? -eq 0 ]; then
        print_success "$desc completed"
    else
        print_error "$desc failed"
        exit 1
    fi
}

# Print header
echo "========================================================================"
echo "  PeCorA (Peptide Correlation Analysis) Complete Workflow"
echo "========================================================================"
echo ""
echo "Parameters:"
echo "  MATLAB Path: $MATLAB_PATH"
echo "  Input Data:  $INPUT_DATA"
echo ""

# Step 0: Check prerequisites
print_status "Checking prerequisites..."

# Check if MATLAB is available
if ! command -v $MATLAB_PATH &> /dev/null; then
    print_error "MATLAB not found at: $MATLAB_PATH"
    print_error "Please install MATLAB or provide correct path as first argument"
    exit 1
fi

# Check if input data file exists
check_file "$INPUT_DATA"
print_success "Prerequisites check passed"

# Create results directory if it doesn't exist
if [ ! -d "results" ]; then
    mkdir results
    print_status "Created results directory"
fi

echo ""
echo "========================================================================"
echo "  STEP 1: Data Preprocessing and Main Analysis"
echo "========================================================================"

# Step 1: Run main PeCorA analysis (includes preprocessing)
run_matlab "addpath('src'); \
            opts = detectImportOptions('$INPUT_DATA'); \
            opts.VariableNamingRule = 'preserve'; \
            t = readtable('$INPUT_DATA', opts); \
            scaled_peptides = PeCorA_preprocessing(t, 'Normalized Area', 100, 'cntrl'); \
            disagree_peptides = PeCorA(scaled_peptides); \
            if ~exist('results', 'dir'), mkdir('results'); end; \
            writetable(disagree_peptides, 'results/PeCorA_results.csv'); \
            writetable(scaled_peptides, 'results/PeCorA_scaled_peptides.csv'); \
            fprintf('\\nPeCorA Analysis Complete: %d peptides analyzed\\n', height(disagree_peptides));" \
           "Main PeCorA Analysis"

echo ""
echo "========================================================================"
echo "  STEP 2: Individual Peptide Visualizations"
echo "========================================================================"

# Step 2: Generate individual peptide plots
run_matlab "run('plot_PeCorA_results')" \
           "Individual Peptide Plots Generation"

echo ""
echo "========================================================================"
echo "  STEP 3: Workflow Demonstration Figures"
echo "========================================================================"

# Step 3: Generate workflow figures
run_matlab "run('generate_figure2_workflow')" \
           "Workflow Demonstration Figures"

echo ""
echo "========================================================================"
echo "  STEP 4: Comprehensive Statistical Summary"
echo "========================================================================"

# Step 4: Generate comprehensive statistics
run_matlab "run('generate_summary_statistics')" \
           "Comprehensive Statistical Summary"

echo ""
echo "========================================================================"
echo "  WORKFLOW COMPLETION SUMMARY"
echo "========================================================================"

# Check and summarize results
print_status "Verifying generated files..."

# Essential output files
EXPECTED_FILES=(
    "results/PeCorA_results.csv"
    "results/PeCorA_scaled_peptides.csv"
    "results/Summary_Statistics.csv"
    "results/Summary_Statistics.txt"
    "results/Figure_2A_Raw_Data.png"
    "results/Figure_2B_Global_Scaling.png"
    "results/Figure_2C_Scaling_Comparison.png"
    "results/README.md"
)

missing_files=0
for file in "${EXPECTED_FILES[@]}"; do
    if [ -f "$file" ]; then
        print_success "✓ $file"
    else
        print_warning "✗ $file (missing)"
        ((missing_files++))
    fi
done

# Count individual peptide plots
plot_count=$(ls results/PeCorA_plot_*.png 2>/dev/null | wc -l)
if [ $plot_count -gt 0 ]; then
    print_success "✓ $plot_count individual peptide plots generated"
else
    print_warning "✗ No individual peptide plots found"
fi

# Final summary
echo ""
if [ $missing_files -eq 0 ]; then
    print_success "WORKFLOW COMPLETED SUCCESSFULLY!"
    echo ""
    echo "Generated Results:"
    echo "  • Main analysis results: results/PeCorA_results.csv"
    echo "  • Scaled peptide data: results/PeCorA_scaled_peptides.csv"
    echo "  • Statistical summaries: results/Summary_Statistics.csv/.txt"
    echo "  • Workflow figures: results/Figure_2*.png"
    echo "  • Individual plots: results/PeCorA_plot_*.png ($plot_count files)"
    echo "  • Documentation: results/README.md"
    echo ""
    echo "Next Steps:"
    echo "  1. Review results/README.md for detailed analysis interpretation"
    echo "  2. Examine results/Summary_Statistics.txt for key findings"
    echo "  3. View workflow figures to understand the analysis process"
    echo "  4. Check individual peptide plots for specific discoveries"
else
    print_warning "WORKFLOW COMPLETED WITH $missing_files MISSING FILES"
    print_warning "Some steps may have failed. Check MATLAB output above for errors."
fi

echo ""
echo "Analysis timestamp: $(date)"
echo "========================================================================" 