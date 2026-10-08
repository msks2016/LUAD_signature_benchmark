#!/bin/bash
set -e

echo "=========================================="
echo "LUAD lncRNA Signature Pipeline"
echo "Started: $(date)"
echo "=========================================="

echo ""
echo "Step 1/4: Environment setup..."
Rscript scripts/R/00_setup.R

echo ""
echo "Step 2/4: Parsing signature files..."
Rscript scripts/R/01_parse_signatures.R

echo ""
echo "Step 3/4: Generating PRISMA flow diagram..."
Rscript scripts/R/03_prisma_diagram.R

echo ""
echo "Step 4/4: Generating Supplementary Table 1..."
Rscript scripts/R/04_supplementary_table1.R

echo ""
echo "=========================================="
echo "PIPELINE COMPLETE"
echo "Finished: $(date)"
echo "=========================================="
echo ""
echo "Outputs:"
echo "  - results/figures/PRISMA_flow_diagram.pdf"
echo "  - results/tables/Supplementary_Table_1.xlsx"
echo "  - data/processed/all_signatures_combined.csv"
