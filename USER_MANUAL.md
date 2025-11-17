# User Manual - Across Chromosomes Phasing Program

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Getting Started](#getting-started)
4. [Command Line Interface](#command-line-interface)
5. [Input Files](#input-files)
6. [Output Files](#output-files)
7. [Advanced Usage](#advanced-usage)
8. [File Format Specifications](#file-format-specifications)
9. [Algorithm Overview](#algorithm-overview)
10. [Troubleshooting](#troubleshooting)
11. [Performance Tips](#performance-tips)
12. [Examples](#examples)

---

## Introduction

The Across Chromosomes Phasing program is designed to perform haplotype phasing across multiple chromosomes. It can process genetic data in HAP or PED format and output phased haplotypes.

### What is Haplotype Phasing?

Haplotype phasing is the process of determining which alleles belong to the same chromosome (haplotype). This program uses cross-chromosome information and, when available, parent-child relationships to improve phasing accuracy.

### Key Capabilities

- **Cross-chromosome phasing**: Uses relationships between chromosomes to improve accuracy
- **Parent-aware processing**: Leverages parent information when available
- **Flexible input/output**: Supports both HAP and PED file formats
- **Selective processing**: Process only specific individuals from your dataset
- **Parallel processing**: Utilizes multiple CPU cores for faster computation

---

## Installation

### System Requirements

- **Operating System**: Linux, macOS, or Windows (with MinGW/MSYS2)
- **CPU**: Multi-core processor recommended
- **Memory**: At least 4GB RAM (more for large datasets)
- **Disk Space**: Sufficient space for input and output files

### Compilation Steps

1. **Ensure prerequisites are installed**:
   ```bash
   # Check for g++
   g++ --version
   
   # Check for make
   make --version
   ```

2. **Compile the program**:
   ```bash
   make all
   ```

3. **Verify installation**:
   ```bash
   ./ProgramPhasing
   # Should show usage information
   
   ./GenerateRandomData
   # Should show usage information
   ```

### Windows Installation

On Windows, you have several options:

1. **Use WSL (Windows Subsystem for Linux)**: Recommended
2. **Install MinGW-w64**: Provides `g++` and `make`
3. **Use MSYS2**: Complete Unix-like environment

---

## Getting Started

### Step 1: Prepare Your Data

Your input data should be organized as:
- **HAP format**: `<base_path>1.hap`, `<base_path>2.hap`, ..., `<base_path>22.hap`
- **PED format**: `<base_path>1.ped`, `<base_path>2.ped`, ..., `<base_path>22.ped`

### Step 2: Run the Program

Basic command:
```bash
./ProgramPhasing -NbIndiv <number> -PathInput <input_path> -PathOutput <output_path>
```

### Step 3: Check Results

Output files will be created at:
- `<output_path>1.ped`, `<output_path>2.ped`, ..., `<output_path>22.ped`

---

## Command Line Interface

### Required Parameters

#### `-NbIndiv <number>`
- **Description**: Number of individuals in your dataset
- **Type**: Integer
- **Range**: 1 to 100,000 (defined by `NBINDIVMAX`)
- **Example**: `-NbIndiv 1000`
- **Notes**: Must match the actual number of individuals in your input files

#### `-PathInput <path>`
- **Description**: Base path for input files
- **Type**: String
- **Format**: Path without chromosome number or extension
- **Example**: `-PathInput data/chr`
- **Notes**: 
  - For HAP: Program looks for `data/chr1.hap`, `data/chr2.hap`, etc.
  - For PED: Program looks for `data/chr1.ped`, `data/chr2.ped`, etc.

#### `-PathOutput <path>`
- **Description**: Base path for output files
- **Type**: String
- **Format**: Path without chromosome number or extension
- **Example**: `-PathOutput results/phased`
- **Notes**: Output files will be `results/phased1.ped`, `results/phased2.ped`, etc.

### Optional Parameters

#### `-InputFormat <HAP|PED>`
- **Description**: Format of input files
- **Type**: String
- **Options**: `HAP` or `PED` (case-insensitive)
- **Default**: `HAP`
- **Example**: `-InputFormat PED`
- **Notes**: 
  - `HAP`: SNP-centric format (one SNP per line)
  - `PED`: Individual-centric format (one individual per line)

#### `-OutputFormat <HAP|PED>`
- **Description**: Format of output files
- **Type**: String
- **Options**: `HAP` or `PED` (case-insensitive)
- **Default**: `PED`
- **Example**: `-OutputFormat PED`
- **Notes**: Currently, HAP output is implemented as PED format

#### `-PathParentInfo <filename>`
- **Description**: Path to file containing parent information
- **Type**: String
- **Format**: Plain text file
- **Example**: `-PathParentInfo parents.txt`
- **File Format**:
  ```
  Individual_ID Parent1_ID Parent2_ID
  0 5 10
  1 -1 -1
  2 15 20
  ```
- **Notes**:
  - Use `-1` to indicate unknown parent
  - Individual IDs must be valid (0 to NbIndiv-1)
  - When parents are known, the program uses `predict_phasing_with_parents_providing_GT`
  - When parents are unknown, the program uses `predict_phasing_without_parents`

#### `-ListIndiv <filename>`
- **Description**: Path to file containing list of individuals to process
- **Type**: String
- **Format**: Plain text file, one ID per line
- **Example**: `-ListIndiv individuals.txt`
- **File Format**:
  ```
  0
  5
  10
  25
  ```
- **Notes**:
  - Individual IDs must be in range [0, NbIndiv-1)
  - Lines starting with `#` are treated as comments
  - Empty lines are ignored
  - If not specified, all individuals are processed

### Complete Command Example

```bash
./ProgramPhasing \
    -NbIndiv 1000 \
    -PathInput data/chr \
    -PathOutput results/phased \
    -InputFormat PED \
    -OutputFormat PED \
    -PathParentInfo parents.txt \
    -ListIndiv individuals.txt
```

---

## Input Files

### HAP Format Files

**File naming**: `<PathInput><chromosome>.hap`
- Example: If `PathInput = "data/chr"`, files are `data/chr1.hap`, `data/chr2.hap`, etc.

**File structure**:
- One SNP per line
- Format: `rsID position allele1_indiv1 allele2_indiv1 allele1_indiv2 allele2_indiv2 ...`
- Alleles are separated by spaces with specific spacing (see format specification)

**Example line**:
```
rs1_0 1234567 0 0 0 1 1 1 0 0 1 1 ...
```

**Reading process**:
1. Program reads rsID and position
2. For each individual, reads two alleles
3. Encodes genotype as 2-bit value (0=00, 1=01, 2=10, 3=11)
4. Stores in compressed format (2 bits per genotype)

### PED Format Files

**File naming**: `<PathInput><chromosome>.ped`
- Example: If `PathInput = "data/chr"`, files are `data/chr1.ped`, `data/chr2.ped`, etc.

**File structure**:
- One individual per line
- Format: `Family_ID Individual_ID Father_ID Mother_ID Sex Phenotype Allele1_SNP1 Allele2_SNP1 Allele1_SNP2 Allele2_SNP2 ...`
- First 6 columns are metadata
- Subsequent columns are allele pairs (one pair per SNP)

**Example line**:
```
0 0 0 0 -9 0 0 0 0 1 1 1 0 1 1 0
```

**Column meanings**:
1. **Family_ID**: Family identifier (can be 0 or any integer)
2. **Individual_ID**: Individual identifier (must be unique)
3. **Father_ID**: Father's Individual_ID (0 if unknown)
4. **Mother_ID**: Mother's Individual_ID (0 if unknown)
5. **Sex**: 1=male, 2=female, -9=unknown
6. **Phenotype**: Phenotype value (0=unknown, -9=missing)
7. **Allele pairs**: Two alleles per SNP (0=reference, 1=alternate)

**Reading process**:
1. Program skips first 6 columns
2. Counts SNPs from first line (number of allele columns / 2)
3. Reads allele pairs for each individual
4. Encodes and stores in same compressed format as HAP

### Parent Information File

**File format**: Plain text, space-separated
```
Individual_ID Parent1_ID Parent2_ID
```

**Rules**:
- One entry per line
- Use `-1` for unknown parent
- Individual IDs must be valid (0 to NbIndiv-1)
- Parent IDs must be valid or -1

**Example**:
```
0 5 10
1 -1 -1
2 15 20
3 25 30
4 -1 -1
```

**Behavior**:
- If both parents are known and valid: Uses `predict_phasing_with_parents_providing_GT`
- If any parent is unknown (-1) or invalid: Uses `predict_phasing_without_parents`

### Individual List File

**File format**: Plain text, one ID per line
```
Individual_ID
```

**Rules**:
- One individual ID per line
- IDs must be in range [0, NbIndiv-1)
- Lines starting with `#` are comments
- Empty lines are ignored

**Example**:
```
# Process these individuals
0
5
10
25
100
```

**Behavior**:
- Only listed individuals are processed
- If file not provided, all individuals are processed

---

## Output Files

### PED Format Output (Default)

**File naming**: `<PathOutput><chromosome>.ped`
- Example: `results/out1.ped`, `results/out2.ped`, etc.

**File structure**:
- One individual per line
- Format: `Family_ID Individual_ID Father_ID Mother_ID Sex Phenotype Allele1_SNP1 Allele2_SNP1 ...`
- Family_ID and Individual_ID are set to the individual's index
- Father_ID, Mother_ID, Sex, Phenotype are set to default values (0, 0, -9, 0)
- Alleles represent phased haplotypes

**Example output line**:
```
0 0 0 0 -9 0 0 0 0 1 1 1 0 1 1 0
```

### Output File Contents

Each output file contains:
- **Header information**: Implicit in PED format (no explicit header)
- **Phased genotypes**: For all processed individuals
- **All SNPs**: For the specified chromosome

### Output Validation

After running the program, verify:
1. **File existence**: All 22 chromosome files should be created
2. **File size**: Should be reasonable (not empty, not corrupted)
3. **Format**: Should match PED format specification
4. **Individual count**: Should match number of processed individuals

---

## Advanced Usage

### Processing with Parent Information

When parent information is available, the program can use it to improve phasing accuracy:

```bash
# Create parent info file
cat > parents.txt << EOF
0 5 10
1 -1 -1
2 15 20
EOF

# Run with parent info
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out \
    -PathParentInfo parents.txt
```

**How it works**:
- For each individual, program checks if parents are known
- If both parents are known: Uses parent-aware phasing algorithm
- If parents are unknown: Falls back to parent-agnostic algorithm

### Selective Individual Processing

Process only specific individuals:

```bash
# Create individual list
cat > individuals.txt << EOF
0
5
10
25
100
EOF

# Process only listed individuals
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out \
    -ListIndiv individuals.txt
```

**Use cases**:
- Testing on a subset of data
- Processing specific individuals of interest
- Debugging specific cases

### Format Conversion

Convert between HAP and PED formats:

```bash
# Read HAP, write PED
./ProgramPhasing -NbIndiv 1000 -PathInput data/hap_chr -PathOutput data/ped_chr \
    -InputFormat HAP -OutputFormat PED

# Read PED, write PED (re-phasing)
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/phased \
    -InputFormat PED -OutputFormat PED
```

### Combining Options

Use all features together:

```bash
./ProgramPhasing \
    -NbIndiv 1000 \
    -PathInput data/chr \
    -PathOutput results/phased \
    -InputFormat PED \
    -OutputFormat PED \
    -PathParentInfo parents.txt \
    -ListIndiv individuals.txt
```

---

## File Format Specifications

### HAP Format Detailed Specification

**File extension**: `.hap`

**Line format**:
```
<rsID> <position> <spaces_and_alleles_for_all_individuals>
```

**Detailed structure**:
- **rsID**: SNP identifier (e.g., `rs1_0`, `rs2_15`)
- **position**: Genomic position in base pairs
- **Individual data**: For each individual:
  - Space character
  - Allele 1 (character '0' or '1')
  - Three space characters
  - Allele 2 (character '0' or '1')
  - Two space characters

**Example** (for 3 individuals):
```
rs1_0 1234567 0   0   0   1   1   1  
```

**Encoding**:
- Genotype 0/0 (homozygous reference) → alleles: 0, 0
- Genotype 0/1 (heterozygous) → alleles: 0, 1
- Genotype 1/1 (homozygous alternate) → alleles: 1, 1
- Missing data → randomly assigned or treated as 0

### PED Format Detailed Specification

**File extension**: `.ped`

**Line format**:
```
<Family_ID> <Individual_ID> <Father_ID> <Mother_ID> <Sex> <Phenotype> <Allele1_SNP1> <Allele2_SNP1> <Allele1_SNP2> <Allele2_SNP2> ...
```

**Column specifications**:

| Column | Type | Description | Example Values |
|--------|------|-------------|----------------|
| Family_ID | Integer | Family identifier | 0, 1, 2, ... |
| Individual_ID | Integer | Individual identifier | 0, 1, 2, ... |
| Father_ID | Integer | Father's ID (0=unknown) | 0, -1, or valid ID |
| Mother_ID | Integer | Mother's ID (0=unknown) | 0, -1, or valid ID |
| Sex | Integer | Sex code | 1 (male), 2 (female), -9 (unknown) |
| Phenotype | Integer | Phenotype value | 0 (unknown), -9 (missing), or value |
| Allele columns | Integer | Allele values | 0 (reference), 1 (alternate) |

**Example line** (4 SNPs):
```
0 0 0 0 -9 0 0 0 0 1 1 1 0 1
```

**Interpretation**:
- Individual 0, no family/parent info
- SNP 1: 0/0 (homozygous reference)
- SNP 2: 0/1 (heterozygous)
- SNP 3: 1/1 (homozygous alternate)
- SNP 4: 0/1 (heterozygous)

### MAF Format Specification

**File extension**: `.maf`

**Format**: Tab-separated values with header

**Structure**:
```
# Chromosome <N> - Minor Allele Frequency file
# SNP_ID	Position	MAF
rs1_0	1234567	0.234567
rs1_1	2345678	0.123456
...
```

**Columns**:
1. **SNP_ID**: SNP identifier (rsID format)
2. **Position**: Genomic position in base pairs
3. **MAF**: Minor Allele Frequency (0.0 to 0.5)

---

## Algorithm Overview

### Phasing Methods

The program implements two main phasing methods:

#### 1. Phasing Without Parents (`predict_phasing_without_parents`)

**When used**:
- No parent information file provided
- Parent information file provided but parents are unknown for the individual
- Invalid parent IDs

**Algorithm**:
- Uses cross-chromosome correlations
- Analyzes relationships between chromosome segments
- Uses population-level information
- Implements iterative refinement

#### 2. Phasing With Parents (`predict_phasing_with_parents_providing_GT`)

**When used**:
- Parent information file provided
- Both parents are known and valid for the individual
- Parent IDs are different from individual ID

**Algorithm**:
- Uses Mendelian inheritance rules
- Leverages parent-child relationships
- Cross-validates with cross-chromosome information
- More accurate when parent data is available

### Processing Pipeline

1. **Initialization**:
   - Parse command line arguments
   - Initialize data structures
   - Set up chromosome dividers

2. **Data Loading**:
   - Read input files (HAP or PED format)
   - Load parent information (if provided)
   - Load individual list (if provided)
   - Allocate memory for genomes

3. **For Each Individual**:
   - Load segments and calculate relationships
   - Determine phasing method (with/without parents)
   - Execute phasing algorithm
   - Store phased haplotypes

4. **Output Generation**:
   - Write phased data to output files
   - Format according to specified output format

### Key Algorithms

- **Relationship detection**: Calculates PI-HAT (relatedness coefficient) between individuals
- **Segment identification**: Identifies shared segments between relatives
- **Cross-chromosome correlation**: Analyzes correlations between chromosome segments
- **Iterative refinement**: Improves phasing through multiple iterations

---

## Troubleshooting

### Common Errors and Solutions

#### Error: "Number of individuals is zero or undefined"

**Cause**: `-NbIndiv` parameter not provided or set to 0

**Solution**:
```bash
# Ensure -NbIndiv is specified
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out
```

#### Error: "Could not open file" or "file end is not found"

**Cause**: Input files don't exist or path is incorrect

**Solution**:
1. Check that input files exist:
   ```bash
   ls data/chr*.hap  # For HAP format
   ls data/chr*.ped  # For PED format
   ```

2. Verify the path is correct:
   ```bash
   # If files are in current directory
   ./ProgramPhasing -NbIndiv 1000 -PathInput ./chr -PathOutput results/out
   
   # If files are in subdirectory
   ./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out
   ```

3. Check file permissions:
   ```bash
   chmod 644 data/chr*.hap  # Make files readable
   ```

#### Error: "Failed to read parent info file"

**Cause**: Parent info file format is incorrect or file doesn't exist

**Solution**:
1. Verify file exists and is readable
2. Check file format (see Parent Information File section)
3. Ensure individual IDs are valid (0 to NbIndiv-1)
4. Check for syntax errors (extra spaces, invalid characters)

#### Error: "Failed to read individual list file"

**Cause**: Individual list file format is incorrect

**Solution**:
1. Verify file format (one ID per line)
2. Ensure IDs are in valid range [0, NbIndiv-1)
3. Check for empty file or invalid entries

#### Error: "Number of individuals is higher than 100000"

**Cause**: `-NbIndiv` exceeds maximum allowed value

**Solution**:
- Reduce number of individuals
- Or modify `NBINDIVMAX` in `types.h` and recompile (not recommended)

#### Warning: "Unknown input format" or "Unknown output format"

**Cause**: Format specification is incorrect

**Solution**:
- Use exactly `HAP` or `PED` (case-insensitive)
- Check spelling: `-InputFormat PED` (not `PEDFILE` or `pedfile`)

### Performance Issues

#### Program runs very slowly

**Possible causes and solutions**:
1. **Large dataset**: Normal for large numbers of individuals/SNPs
   - Consider processing subsets using `-ListIndiv`
   - Use more CPU cores (OpenMP should detect automatically)

2. **Insufficient memory**: Program may swap to disk
   - Close other applications
   - Reduce number of individuals or SNPs
   - Add more RAM

3. **Single-threaded execution**: OpenMP not working
   - Check OpenMP installation: `echo $OMP_NUM_THREADS`
   - Set threads: `export OMP_NUM_THREADS=4`

#### Out of memory errors

**Solution**:
- Reduce `-NbIndiv` to a smaller number
- Process chromosomes separately (modify code)
- Use a machine with more RAM

### File Format Issues

#### PED file reading errors

**Symptoms**: "WARNING: Line X has insufficient columns"

**Solution**:
1. Verify all lines have the same number of columns
2. Check for missing alleles (should be pairs)
3. Ensure proper delimiter (spaces or tabs)
4. Check for special characters or encoding issues

#### HAP file reading errors

**Symptoms**: "file end is not found" or incorrect SNP counts

**Solution**:
1. Verify file format matches specification exactly
2. Check spacing between alleles (critical for HAP format)
3. Ensure all individuals have data for all SNPs
4. Check for line ending issues (Windows vs Unix)

---

## Performance Tips

### Optimization Strategies

1. **Use appropriate file format**:
   - HAP format: Faster to read for SNP-centric operations
   - PED format: Better for individual-centric operations

2. **Process subsets**:
   - Use `-ListIndiv` to process only needed individuals
   - Reduces memory usage and processing time

3. **Parallel processing**:
   - Ensure OpenMP is enabled (compiled with `-fopenmp`)
   - Set `OMP_NUM_THREADS` environment variable:
     ```bash
     export OMP_NUM_THREADS=8
     ./ProgramPhasing ...
     ```

4. **Memory management**:
   - Close other applications
   - Process in batches if dataset is very large

5. **File I/O optimization**:
   - Use fast storage (SSD preferred)
   - Ensure sufficient disk space
   - Use local storage (not network drives)

### Expected Performance

**Typical processing times** (approximate, depends on hardware):
- 100 individuals, 1000 SNPs/chr: ~1-5 minutes
- 1000 individuals, 1000 SNPs/chr: ~10-30 minutes
- 10000 individuals, 1000 SNPs/chr: ~1-3 hours

**Memory usage**:
- Approximately: `NbIndiv × SNPs_per_chr × 2 bits / 8 bytes`
- Plus overhead for algorithm data structures

---

## Examples

### Example 1: Basic Phasing

**Scenario**: You have HAP format files and want to phase all individuals.

**Files**:
- `data/chr1.hap`, `data/chr2.hap`, ..., `data/chr22.hap`
- 1000 individuals

**Command**:
```bash
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/phased
```

**Output**:
- `results/phased1.ped`, `results/phased2.ped`, ..., `results/phased22.ped`

### Example 2: Phasing with Parent Information

**Scenario**: You have parent information for some individuals.

**Files**:
- Input: `data/chr1.hap`, etc.
- Parent info: `parents.txt`
  ```
  0 5 10
  1 -1 -1
  2 15 20
  ```

**Command**:
```bash
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/phased \
    -PathParentInfo parents.txt
```

**Result**:
- Individual 0: Uses parent-aware phasing (parents 5 and 10)
- Individual 1: Uses standard phasing (no parents)
- Individual 2: Uses parent-aware phasing (parents 15 and 20)

### Example 3: Selective Processing

**Scenario**: You want to phase only specific individuals.

**Files**:
- Input: `data/chr1.hap`, etc.
- Individual list: `individuals.txt`
  ```
  0
  5
  10
  25
  ```

**Command**:
```bash
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/phased \
    -ListIndiv individuals.txt
```

**Result**:
- Only individuals 0, 5, 10, and 25 are processed
- Output files contain only these individuals

### Example 4: PED Format Workflow

**Scenario**: Complete workflow with PED format.

**Step 1 - Generate PED data**:
```bash
./GenerateRandomData test_data/chr PED
```

**Step 2 - Phase the data**:
```bash
./ProgramPhasing -NbIndiv 1000 -PathInput test_data/chr -PathOutput results/phased \
    -InputFormat PED -OutputFormat PED
```

**Step 3 - Verify results**:
```bash
head -5 results/phased1.ped
```

### Example 5: Complex Workflow

**Scenario**: Full workflow with all features.

**Step 1 - Generate test data**:
```bash
./GenerateRandomData test_data/chr PED
```

**Step 2 - Create parent info file**:
```bash
cat > parents.txt << EOF
0 5 10
1 -1 -1
2 15 20
3 25 30
EOF
```

**Step 3 - Create individual list**:
```bash
cat > individuals.txt << EOF
0
1
2
3
EOF
```

**Step 4 - Run phasing**:
```bash
./ProgramPhasing \
    -NbIndiv 1000 \
    -PathInput test_data/chr \
    -PathOutput results/phased \
    -InputFormat PED \
    -OutputFormat PED \
    -PathParentInfo parents.txt \
    -ListIndiv individuals.txt
```

### Example 6: Format Conversion

**Scenario**: Convert HAP files to PED format.

**Command**:
```bash
./ProgramPhasing -NbIndiv 1000 -PathInput data/hap_chr -PathOutput data/ped_chr \
    -InputFormat HAP -OutputFormat PED
```

**Note**: This reads HAP files and writes PED files, effectively converting the format.

---

## Additional Resources

### Related Documentation

- **README.md**: Project overview and quick start
- **README_GENERATOR.md**: Data generator documentation
- **README_PED_FORMAT.md**: Detailed PED format specification

### Getting Help

If you encounter issues:
1. Check this manual's Troubleshooting section
2. Verify file formats match specifications
3. Check command line syntax
4. Review example commands

### Reporting Issues

When reporting issues, include:
- Command used (exact command line)
- Error messages (complete output)
- Input file format and size
- System information (OS, compiler version)
- Sample input files (if possible)

---

## Appendix

### A. File Size Estimates

**HAP files**:
- Size ≈ `(SNPs × (rsID_length + position_length + individuals × 8_bytes))`
- Example: 1000 SNPs, 1000 individuals ≈ 8-10 MB per chromosome

**PED files**:
- Size ≈ `(individuals × (6_metadata_columns + SNPs × 2_alleles) × bytes_per_value)`
- Example: 1000 individuals, 1000 SNPs ≈ 8-10 MB per chromosome

### B. Memory Requirements

**Approximate memory usage**:
```
Base memory + (NbIndiv × SNPs_per_chr × 2 bits / 8) + algorithm overhead
```

**Example**: 1000 individuals, 1000 SNPs/chr
- Genome data: ~250 KB per chromosome
- 22 chromosomes: ~5.5 MB
- Algorithm overhead: ~50-100 MB
- **Total**: ~100-200 MB

### C. Supported Chromosomes

The program processes chromosomes 1-22 (autosomes). Chromosome 0 is used for metadata.

### D. Constants and Limits

Defined in `types.h`:
- `NSNPPERCHR`: Maximum SNPs per chromosome (60000)
- `MAXPOP`: Maximum population size (435188)
- `NBINDIVMAX`: Maximum individuals (100000)
- `MAXCLOSERELAT`: Maximum close relatives (6000)

---

**End of User Manual**
