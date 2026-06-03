#!/usr/bin/env bash

# ==============================================================================
# GOParsimony Package Wrapper
# ==============================================================================
# 
# Script to run the GOParsimony Package
# 
# Automates the generation of visually distinct, precedence-filtered GO 
# Directed Acyclic Graphs (DAGs) and hierarchical heatmaps.
# 
# 
# Copyright 2026. Paul Martin Harrison. 
# License: 3-clause BSD license 
# 
#
# Core Inputs:
#   -b : Background annotation file (Simple 2-column or GAF)
#   -s : Sample annotation file (Simple 1-column list of IDs or accessions)
#   -o : GO ontology file (Defaults to go-basic.obo)
#
# Output:
#   A suite of publication-ready PDFs including Precedence DAGs, concise 
#   skip-edge graphs, Jaccard-distance heatmaps, and comparative statistics.
# 
#  To run and get help: 
#    sh run-GOParsimony.sh -h 
#     OR
#    chmod u+x run-GOParsimony.sh; ./run-GOParsimony.sh -h 
# 
# ==============================================================================

# Stop execution immediately if any command fails
set -e

# ==============================================================================
# DEFAULT VARIABLES
# ==============================================================================
OBO_FILE=""
SLIMS_FILE=""
SAMPLE_FILE=""
BACKGROUND_FILE=""
HIGHLIGHT_FILE=""
PREFIX=""
MAX_PREC=5
OUTPUT_EXPANDED=false

# Name of your compiled C executable and Python script
C_EXEC="./GOParsimony" 
PY_SCRIPT="precedence-graphs.py" 

# ==============================================================================
# FIND PYTHON 3 EXECUTABLE
# ==============================================================================
# 1. First, check if python3 is in the user's PATH (respects conda/brew/venv)
if command -v python3 &>/dev/null; then
    PYTHON_EXEC="python3"
# 2. Next, check if 'python' defaults to Python 3
elif command -v python &>/dev/null && python -c 'import sys; exit(0 if sys.version_info >= (3,0) else 1)' &>/dev/null; then
    PYTHON_EXEC="python"
# 3. Fallback: manually search /usr/bin/ for a python3 executable
else
    FOUND_PY=$(ls /usr/bin/python3* 2>/dev/null | grep -E '^/usr/bin/python3(\.[0-9]+)?$' | head -n 1 || true)
    if [ -n "$FOUND_PY" ]; then
        PYTHON_EXEC="$FOUND_PY"
    else
        echo "Error: Python 3 executable not found in PATH or /usr/bin/. Please install Python 3."
        exit 1
    fi
fi

# ==============================================================================
# USAGE FUNCTION
# ==============================================================================
usage() {
    echo "======================================================================="
    echo " script to run GOParsimony complete pipeline " 
    echo ""     
    echo " Usage: $0 -s <sample_file> -b <background_file> [OPTIONS]"
    echo "======================================================================="
    echo " Required:"
    echo "   -b FILE    Path to the background annotations file"
    echo "   -s FILE    Path to the list of sample IDs/accessions (one-column format, must be consistent with the background annotations file)"
    echo ""
    echo " Optional:"
    echo "   -o FILE    Path to the go.obo file (Default: looks for go-basic.obo)"
    echo "   -m FILE    Path to the GOSlims file"
    echo "   -l FILE    Path to a list of GO terms to highlight"
    echo "   -p STR     Output file prefix (Default: None)"
    echo "   -x INT     Maximum precedence threshold for clustering (Default: 5)"
    echo "   -e         Output the expanded annotation lists (Default: deleted after run)"
    echo "   -h         Display this help message"
    echo "======================================================================="
    exit 1
}

# ==============================================================================
# PARSE COMMAND LINE ARGUMENTS
# ==============================================================================
while getopts "o:s:b:m:l:p:x:he" opt; do
    case "$opt" in
        o) OBO_FILE="$OPTARG" ;;
        s) SAMPLE_FILE="$OPTARG" ;;
        b) BACKGROUND_FILE="$OPTARG" ;;
        m) SLIMS_FILE="$OPTARG" ;;
        l) HIGHLIGHT_FILE="$OPTARG" ;;
        p) PREFIX="$OPTARG" ;;
        x) MAX_PREC="$OPTARG" ;;
        e) OUTPUT_EXPANDED=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

# ==============================================================================
# VALIDATION
# ==============================================================================
# Check if required arguments are provided
if [ -z "$SAMPLE_FILE" ] || [ -z "$BACKGROUND_FILE" ]; then
    echo "Error: Missing required sample or background files."
    usage
fi

# Check if the C executable exists
if [ ! -f "$C_EXEC" ]; then
    echo "Error: C executable '$C_EXEC' not found! Please compile your C program first."
    exit 1
fi

# Check if the Python script exists
if [ ! -f "$PY_SCRIPT" ]; then
    echo "Error: Python script '$PY_SCRIPT' not found in the current directory."
    exit 1
fi

# ==============================================================================
# 1. RUN THE C PROGRAM
# ==============================================================================
echo "-------------------------------------------------------------------------"
echo " Step 1: Running C Parser..."
echo "-------------------------------------------------------------------------"

# Construct the C command dynamically using explicit flags
C_CMD="$C_EXEC -b \"$BACKGROUND_FILE\" -s \"$SAMPLE_FILE\""

if [ -n "$OBO_FILE" ]; then
    C_CMD="$C_CMD -o \"$OBO_FILE\""
fi

if [ -n "$SLIMS_FILE" ]; then
    C_CMD="$C_CMD -m \"$SLIMS_FILE\""
fi

if [ -n "$PREFIX" ]; then
    C_CMD="$C_CMD -p \"$PREFIX\""
fi

if [ "$OUTPUT_EXPANDED" = true ]; then
    C_CMD="$C_CMD -e"
fi

echo "Executing: $C_CMD"
eval $C_CMD

# ==============================================================================
# 2. RUN THE PYTHON VISUALIZATION SCRIPT
# ==============================================================================
echo ""
echo "-------------------------------------------------------------------------"
echo " Step 2: Running Python Visualization..."
echo " Using Python Executable: $PYTHON_EXEC"
echo "-------------------------------------------------------------------------"

# Format the expected filenames based on whether a prefix was provided
if [ -n "$PREFIX" ]; then
    case "$PREFIX" in
        *_) PREFIX_SEP="$PREFIX" ;;
        *)  PREFIX_SEP="${PREFIX}_" ;;
    esac
else
    PREFIX_SEP=""
fi

RES_FILE="${PREFIX_SEP}Enrichment_Results.tsv"
EDGES_FILE="${PREFIX_SEP}Precedence_Edges.tsv"
SIG_PROTS_FILE="${PREFIX_SEP}Significantly_Annotated_List.tsv"

# Verify that the C program successfully created the expected output files
if [ ! -f "$RES_FILE" ] || [ ! -f "$EDGES_FILE" ] || [ ! -f "$SIG_PROTS_FILE" ]; then
    echo "Error: The C program did not generate the expected .tsv files:"
    echo "  Expected: $RES_FILE, $EDGES_FILE, $SIG_PROTS_FILE"
    exit 1
fi

# Construct the base Python command feeding the prefixed filenames into it
PY_CMD="$PYTHON_EXEC $PY_SCRIPT --results \"$RES_FILE\" --edges \"$EDGES_FILE\" --sig_prots \"$SIG_PROTS_FILE\" --out \"$PREFIX\" --max_prec $MAX_PREC"

# Append highlight file if provided
if [ -n "$HIGHLIGHT_FILE" ]; then
    PY_CMD="$PY_CMD --highlight \"$HIGHLIGHT_FILE\""
fi

echo "Executing: $PY_CMD"
eval $PY_CMD

echo ""
echo "======================================================================="
echo " Pipeline Complete! All outputs have been saved successfully."
echo "======================================================================="

