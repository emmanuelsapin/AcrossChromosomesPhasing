# User Manual - Haplotype Phasing Engine

Complete user guide for the Haplotype Phasing Engine software.

## Table of Contents

1. [Introduction](#introduction)
2. [Installation Guide](#installation-guide)
3. [Getting Started](#getting-started)
4. [Command Line Reference](#command-line-reference)
5. [Input Data Format](#input-data-format)
6. [Output Data Format](#output-data-format)
7. [Configuration Options](#configuration-options)
8. [Advanced Usage](#advanced-usage)
9. [Troubleshooting](#troubleshooting)
10. [FAQ](#frequently-asked-questions)

## Introduction

### What is Haplotype Phasing?

Haplotype phasing is the process of determining which alleles on homologous chromosomes are inherited together. This software uses relative-based inference algorithms to phase haplotypes by analyzing relationships between individuals in a population.

### Software Overview

The Haplotype Phasing Engine is a professional-grade tool for:
- Processing large-scale genomic datasets
- Identifying genetic relationships (PIHAT computation)
- Phasing haplotypes using relative information
- Generating phased output in standard formats

### Key Features

- **Scalability**: Handles datasets with up to 100,000+ individuals
- **Accuracy**: Advanced algorithms for high-quality phasing
- **Performance**: Multi-threaded processing with OpenMP
- **Flexibility**: Multiple input/output formats
- **Reliability**: Comprehensive error handling and validation

## Installation Guide

### Step 1: Verify Prerequisites

Check your system has the required tools:

```bash
# Check C++ compiler
g++ --version  # Should be 4.9 or higher

# Check OpenMP support
g++ -fopenmp --version

# Check Make
make --version
```

### Step 2: Download/Clone the Project

```bash
# If using git
git clone [repository-url]
cd Fileforgithub

# Or extract from archive
unzip haplotype-phasing-engine.zip
cd Fileforgithub
```

### Step 3: Compile the Software

**Option A: Modular Build (Recommended)**
```bash
make -f Makefile_modular
```

**Option B: Object-Oriented Build**
```bash
make -f Makefile_oo
```

**Expected Output:**
```
g++ -std=c++11 -O0 -g -fopenmp -Wall -Wextra -I./include -c main_oo.cpp -o main_oo.o
...
g++ -std=c++11 -O0 -g -fopenmp -Wall -Wextra -I./include -o ProgramPhasing_Modular ...
```

### Step 4: Verify Installation

```bash
./ProgramPhasing_Modular --help
# Or simply run without arguments to see error message (confirms it's working)
```

## Getting Started

### Your First Run

1. **Prepare Input Data**
   - Organize your HAP files: `1.hap`, `2.hap`, ..., `22.hap`
   - Place them in a directory (e.g., `./input/`)

2. **Run the Program**
   ```bash
   ./ProgramPhasing_Modular \
     -NbIndiv 100 \
     -PathInput ./input/ \
     -PathOutput ./output/
   ```

3. **Check Results**
   - Output files will be in `./output/` directory
   - Files: `1.ped`, `2.ped`, ..., `22.ped`

### Quick Example

```bash
# Create directories
mkdir -p input output

# Copy your HAP files to input/
cp your_data/*.hap input/

# Run phasing
./ProgramPhasing_Modular \
  -NbIndiv 500 \
  -PathInput ./input/ \
  -PathOutput ./output/ \
  -Verbose 1

# Check results
ls -lh output/
```

## Command Line Reference

### Syntax

```bash
./ProgramPhasing_Modular [REQUIRED_ARGS] [OPTIONAL_ARGS]
```

### Required Arguments

#### `-NbIndiv <number>`
- **Description**: Number of individuals in the dataset
- **Type**: Integer
- **Range**: 1 to 100,000
- **Example**: `-NbIndiv 5000`

#### `-PathInput <path>`
- **Description**: Path prefix for input HAP files
- **Type**: String (directory path)
- **Note**: Should end with `/` or `\` depending on OS
- **Example**: `-PathInput /data/genomes/`
- **File Naming**: Program expects `<PathInput><chromosome>.hap`
  - For `-PathInput ./data/`, files should be `./data/1.hap`, `./data/2.hap`, etc.

#### `-PathOutput <path>`
- **Description**: Path prefix for output PED files
- **Type**: String (directory path)
- **Note**: Directory will be created if it doesn't exist
- **Example**: `-PathOutput ./results/`
- **File Naming**: Output files will be `<PathOutput><chromosome>.ped`

### Optional Arguments

#### `-PathMAF <path>`
- **Description**: Path to Minor Allele Frequency (MAF) file
- **Type**: String (file path)
- **Default**: `/pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt`
- **Example**: `-PathMAF ./MAF.txt`
- **Format**: Text file with chromosome, SNP index, and MAF value per line

#### `-Verbose <0|1>`
- **Description**: Enable verbose output
- **Type**: Integer (0 or 1)
- **Default**: 1 (enabled)
- **Example**: `-Verbose 0` for quiet mode
- **Verbose Output Includes**:
  - Progress messages
  - Individual processing status
  - Relative identification results
  - Execution statistics

#### `-ListIndiv <path>`
- **Description**: Path to a file containing a list of individual IDs to process
- **Type**: String (file path)
- **Default**: Not specified (all individuals are processed)
- **Example**: `-ListIndiv ./liste_individus.txt`
- **File Format**: Text file with one individual ID per line (0-based indexing)
- **Features**:
  - Empty lines are ignored
  - Comments (lines starting with `#`) are ignored
  - Individual IDs must be in range [0, NbIndiv)
  - Invalid IDs are skipped with a warning
- **Use Cases**:
  - Process only a subset of individuals for testing
  - Reduce computation time by focusing on specific individuals
  - Reproduce analyses on specific individuals
- **Example File Content**:
  ```
  # List of individuals to process
  0
  10
  20
  30
  50
  ```
- **Note**: When specified, only the listed individuals are processed and written to output files. All genomic data is still loaded in memory (required for relationship calculations).

### Complete Example

```bash
./ProgramPhasing_Modular \
  -NbIndiv 10000 \
  -PathInput /large/dataset/genomes/ \
  -PathOutput /results/phased_genomes/ \
  -PathMAF /data/MAF.txt \
  -Verbose 1
```

**Example with Individual List:**
```bash
./ProgramPhasing_Modular \
  -NbIndiv 10000 \
  -PathInput /large/dataset/genomes/ \
  -PathOutput /results/phased_genomes/ \
  -ListIndiv ./liste_individus.txt \
  -Verbose 1
```

## Input Data Format

### HAP File Format

The program expects input files in HAP format (similar to IMPUTE2 format):

**File Structure:**
```
rsID position allele0 allele1 individual1_allele0 individual1_allele1 ...
```

**Example:**
```
rs123 1000 A T 0 1
rs456 2000 G C 1 1
rs789 3000 T A 0 0
```

**File Naming Convention:**
- Chromosome 1: `<PathInput>1.hap`
- Chromosome 2: `<PathInput>2.hap`
- ...
- Chromosome 22: `<PathInput>22.hap`

### MAF File Format (Optional)

If using `-PathMAF`, the file should contain:
```
chromosome snp_index maf_value
```

**Example:**
```
1 0 0.25
1 1 0.15
2 0 0.30
...
```

### Data Requirements

- **Chromosomes**: 1-22 (autosomal)
- **SNPs**: Variable per chromosome
- **Individuals**: Must match `-NbIndiv` parameter
- **Encoding**: Alleles encoded as 0/1 (reference/alternate)

## Output Data Format

### PED File Format

Output files are in PLINK PED format:

**Format:**
```
FamilyID IndividualID PaternalID MaternalID Sex Phenotype Genotype1_Allele1 Genotype1_Allele2 ...
```

**Example:**
```
0 0 0 0 -9 0 0 1 1 1 0 0 ...
0 1 0 0 -9 0 1 0 0 1 1 1 ...
```

**File Naming:**
- Chromosome 1: `<PathOutput>1.ped`
- Chromosome 2: `<PathOutput>2.ped`
- ...
- Chromosome 22: `<PathOutput>22.ped`

### Output Interpretation

- **FamilyID**: Set to individual index (0-based)
- **IndividualID**: Set to individual index (0-based)
- **PaternalID/MaternalID**: Set to 0 (not used)
- **Sex**: Set to -9 (unknown)
- **Phenotype**: Set to 0 (not used)
- **Genotypes**: Phased alleles (0/1 encoding)

## Configuration Options

### Environment Variables

**OpenMP Thread Configuration:**
```bash
export OMP_NUM_THREADS=8  # Use 8 threads
```

**Memory Configuration:**
```bash
# For large datasets, ensure sufficient memory
ulimit -v unlimited
```

### Algorithm Parameters

Current implementation uses:
- **PIHAT Threshold**: 0.33 (default)
- **Maximum Generations**: 1 (configurable in code)
- **Breakpoint Count**: 1000 (configurable)

*Note: Advanced parameters require code modification*

## Advanced Usage

### Batch Processing

Process multiple datasets:

```bash
#!/bin/bash
for dataset in dataset1 dataset2 dataset3; do
  ./ProgramPhasing_Modular \
    -NbIndiv 5000 \
    -PathInput ./${dataset}/input/ \
    -PathOutput ./${dataset}/output/ \
    -Verbose 0
done
```

### Parallel Execution

Run multiple chromosomes in parallel:

```bash
# Process chromosomes 1-11
for chr in {1..11}; do
  ./ProgramPhasing_Modular \
    -NbIndiv 10000 \
    -PathInput ./input/ \
    -PathOutput ./output/chr${chr}/ \
    -Verbose 0 &
done
wait
```

### Processing a Subset of Individuals

Use `-ListIndiv` to process only specific individuals. This is useful for:
- Testing and debugging on a small subset
- Focusing on specific individuals of interest
- Reducing computation time
- Reproducing analyses on selected individuals

**Step 1: Create a list file** (`liste_individus.txt`):
```
# List of individuals to process
# Format: One ID per line (0-based indexing)
# Comments and empty lines are ignored

0
10
20
30
50
100
```

**Step 2: Run with the list:**
```bash
./ProgramPhasing_Modular \
  -NbIndiv 1000 \
  -PathInput ./input/ \
  -PathOutput ./output/ \
  -ListIndiv ./liste_individus.txt \
  -Verbose 1
```

**Output:**
- Only the specified individuals are processed
- Output files contain only these individuals
- Processing time is reduced proportionally
- Progress shows: `start individual: X (Y/Z)` where Y is current position and Z is total

**Important Notes:**
- All genomic data is still loaded in memory (required for relationship calculations)
- Individual IDs must be in range [0, NbIndiv)
- Invalid IDs are skipped with a warning
- If no valid IDs are found, the program exits with an error

### Logging Output

Save execution log:

```bash
./ProgramPhasing_Modular \
  -NbIndiv 5000 \
  -PathInput ./input/ \
  -PathOutput ./output/ \
  -Verbose 1 2>&1 | tee execution.log
```

## Troubleshooting

### Problem: "File not found" Error

**Symptoms:**
```
Error: Cannot access input path: ./input/
```

**Solutions:**
1. Check path exists: `ls -la ./input/`
2. Verify file naming: Files should be `1.hap`, `2.hap`, etc.
3. Check permissions: `chmod +r ./input/*.hap`
4. Use absolute paths if relative paths fail

### Problem: "Memory allocation failed"

**Symptoms:**
```
Error: Memory allocation failed for chromosome X
```

**Solutions:**
1. Reduce number of individuals processed at once
2. Increase system swap space
3. Process chromosomes separately
4. Use a machine with more RAM

### Problem: "Invalid number of individuals"

**Symptoms:**
```
ERROR: Number of individuals is zero or undefined
```

**Solutions:**
1. Check `-NbIndiv` argument is provided
2. Verify value is between 1 and 100,000
3. Ensure value matches actual data

### Problem: "Failed to read individual list file"

**Symptoms:**
```
ERROR: Could not open file ./liste_individus.txt for reading individual list
ERROR: Failed to read individual list file. Exiting.
```

**Solutions:**
1. Verify the file path is correct
2. Check file permissions: `chmod +r ./liste_individus.txt`
3. Ensure the file exists: `ls -la ./liste_individus.txt`
4. Check file format: one ID per line, starting from 0

### Problem: "No valid individual IDs found"

**Symptoms:**
```
ERROR: No valid individual IDs found in file ./liste_individus.txt
```

**Solutions:**
1. Verify IDs are in range [0, NbIndiv)
2. Check file format (one ID per line)
3. Ensure IDs are numeric (not text)
4. Remove any invalid entries from the list file

### Problem: Slow Performance

**Symptoms:**
- Processing takes much longer than expected

**Solutions:**
1. Enable OpenMP: `export OMP_NUM_THREADS=8`
2. Use `-Verbose 0` to reduce I/O overhead
3. Ensure input files are on fast storage (SSD)
4. Check system resources aren't constrained

### Problem: Incorrect Output

**Symptoms:**
- Output files are empty or malformed

**Solutions:**
1. Verify input file format matches HAP format
2. Check chromosome numbering (1-22)
3. Ensure SNP counts are consistent
4. Validate with smaller test dataset first

## Frequently Asked Questions

### Q: What is the maximum number of individuals supported?

**A:** The software is configured for up to 100,000 individuals (NBINDIVMAX). For larger datasets, modify `Constants.h` and recompile.

### Q: Can I process only specific chromosomes?

**A:** The current version processes all chromosomes 1-22. For selective processing, modify the source code or use separate input directories.

### Q: Can I process only specific individuals?

**A:** Yes! Use the `-ListIndiv` parameter with a file containing individual IDs (one per line). This allows you to:
- Test on a small subset before processing the full dataset
- Focus on specific individuals of interest
- Reduce computation time
- Reproduce analyses on selected individuals

Example:
```bash
# Create liste_individus.txt with IDs: 0, 10, 20, 30
./ProgramPhasing_Modular \
  -NbIndiv 1000 \
  -PathInput ./input/ \
  -PathOutput ./output/ \
  -ListIndiv ./liste_individus.txt
```

### Q: How long does phasing take?

**A:** Processing time depends on:
- Number of individuals (linear scaling)
- Number of SNPs per chromosome
- Hardware (CPU cores, memory, storage speed)
- Typical: 1,000 individuals ≈ 5-10 minutes, 10,000 ≈ 30-60 minutes

### Q: What if I don't have a MAF file?

**A:** The program will attempt to use the default MAF file path. If not found, MAF values will be computed from the input data (slower but functional).

### Q: Can I use VCF or PED input formats?

**A:** The current version supports HAP format. VCF/PED support is planned for future releases. Use conversion tools to convert your data to HAP format.

### Q: How do I interpret PIHAT values?

**A:** PIHAT values indicate genetic relatedness:
- 0.0: Unrelated
- 0.25: First-degree relatives (parent-child, siblings)
- 0.5: Identical twins or duplicate samples
- Higher values: Closer relationships

### Q: Is the output deterministic?

**A:** Yes, given the same input and parameters, the output should be identical. Random number generation is only used in specific test scenarios.

### Q: Can I resume interrupted processing?

**A:** Currently, the program processes all individuals in one run. For resumability, implement checkpointing or process individuals in batches.

### Q: How much disk space do I need?

**A:** Approximately:
- Input: ~1-2 GB per 10,000 individuals (depends on SNP density)
- Output: Similar size to input
- Temporary: Minimal (in-memory processing)

### Q: What operating systems are supported?

**A:** 
- Linux: Fully supported
- macOS: Supported (may need OpenMP installation)
- Windows: Supported via WSL or MinGW

## Additional Resources

### Documentation Files

- `README.md`: Project overview and quick start
- `README_ARCHITECTURE.md`: Detailed architecture documentation
- `README_MODULAR.md`: Modular structure guide
- `STRUCTURE_MODULAIRE.md`: Module dependencies

### Getting Help

1. Check this manual first
2. Review error messages (they're descriptive)
3. Check documentation files
4. Verify your data format matches specifications
5. Test with a small dataset first

### Performance Tips

1. **Use SSD storage** for input/output files
2. **Configure OpenMP threads** appropriately for your CPU
3. **Process in batches** for very large datasets
4. **Use quiet mode** (`-Verbose 0`) for production runs
5. **Monitor system resources** during processing

---

**Manual Version**: 1.0  
**Last Updated**: 2024  
**Software Version**: 2.0


