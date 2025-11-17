# Random Data Generator - User Manual

## Description

The `GenerateRandomData` program generates random HAP, PED, and MAF files for testing the haplotype phasing algorithm. It creates synthetic genetic data with realistic properties following population genetics principles.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Output Files](#output-files)
- [File Format Details](#file-format-details)
- [Parameters](#parameters)
- [Examples](#examples)
- [Integration with Phasing Program](#integration-with-phasing-program)
- [Technical Details](#technical-details)
- [Troubleshooting](#troubleshooting)

## Overview

The data generator creates synthetic genetic data that mimics real population genetics data. It generates:

- **Genotype data**: For a specified number of individuals and SNPs
- **Multiple formats**: HAP or PED format files
- **MAF information**: Minor Allele Frequency files for each chromosome
- **Realistic distributions**: Uses population genetics principles

## Features

- **Configurable individuals**: Generates data for 1000 individuals (configurable in code)
- **Configurable SNPs**: Creates 1000 SNPs per chromosome (configurable in code)
- **22 chromosomes**: Generates data for chromosomes 1-22 (autosomes)
- **Realistic allele frequency distribution**: Minor Allele Frequencies (MAF) are directly generated in the range [0.05, 0.5], ensuring all SNPs have a MAF of at least 5%
- **Hardy-Weinberg equilibrium**: Genotypes are generated according to Hardy-Weinberg equilibrium principles
- **Multiple output formats**: Supports both HAP and PED file formats
- **Consistent data**: Same MAF values used for both genotype and MAF files

## Compilation

To compile the generator:

```bash
make generator
```

Or compile everything including the main program:

```bash
make all
```

## Usage

```bash
./GenerateRandomData <output_path> [format]
```

### Parameters

- `output_path`: Base path for output files (e.g., `data/chr` or `test_data/chr`)
- `format`: Optional format specification - `HAP` (default) or `PED`

### Examples

#### Generate HAP files (default)
```bash
./GenerateRandomData data/chr
```
or
```bash
./GenerateRandomData data/chr HAP
```

This will generate:
- `data/chr1.hap`, `data/chr2.hap`, ..., `data/chr22.hap` (HAP files)
- `data/chr1.maf`, `data/chr2.maf`, ..., `data/chr22.maf` (MAF files)

#### Generate PED files
```bash
./GenerateRandomData data/chr PED
```

This will generate:
- `data/chr1.ped`, `data/chr2.ped`, ..., `data/chr22.ped` (PED files)
- `data/chr1.maf`, `data/chr2.maf`, ..., `data/chr22.maf` (MAF files)

## Output Files

### HAP Files

HAP files contain haplotype data in the following format:
- Each line represents one SNP
- Format: `rsID position allele1_individual1 allele2_individual1 allele1_individual2 allele2_individual2 ...`
- Alleles are encoded as '0' (reference) or '1' (alternate)

Example line:
```
rs1_0 1234567 0 0 0 1 1 1 ...
```

### PED Files

PED files contain genotype data in the following format:
- Each line represents one individual
- Format: `Family_ID Individual_ID Father_ID Mother_ID Sex Phenotype Allele1_SNP1 Allele2_SNP1 Allele1_SNP2 Allele2_SNP2 ...`
- First 6 columns are metadata (family, individual, parents, sex, phenotype)
- Subsequent columns are allele pairs (one pair per SNP)
- Alleles are encoded as `0` (reference) or `1` (alternate)

Example line:
```
0 0 0 0 -9 0 0 0 0 1 1 1 0 1 1 0
```
This represents individual 0 with:
- SNP 1: 0/0 (homozygous reference)
- SNP 2: 0/1 (heterozygous)
- SNP 3: 1/1 (homozygous alternate)
- SNP 4: 0/1 (heterozygous)

### MAF Files

MAF files contain Minor Allele Frequency information:
- Header line with column names
- Format: `SNP_ID\tPosition\tMAF`
- Tab-separated values

Example:
```
# Chromosome 1 - Minor Allele Frequency file
# SNP_ID	Position	MAF
rs1_0	1234567	0.234567
rs1_1	2345678	0.123456
...
```

## File Format Details

### HAP File Format

The HAP file format matches what the phasing program expects:
- SNP identifier (rsID)
- Position in base pairs
- For each individual: two alleles (0 or 1) separated by spaces

The reading code expects:
- First character: space or SNP identifier
- For each individual: space, allele1, space, space, space, allele2, space, space

### Genotype Generation

Genotypes are generated using Hardy-Weinberg equilibrium:
- **Homozygous reference (0/0)**: Probability = p² where p = 1 - MAF
- **Heterozygous (0/1)**: Probability = 2pq where q = MAF
- **Homozygous alternate (1/1)**: Probability = q²

## Integration with Phasing Program

After generating the data, you can use it with the phasing program:

### Using HAP format (default)
```bash
# Generate HAP data
./GenerateRandomData data/chr HAP

# Run phasing program with HAP input (default)
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out
```

### Using PED format
```bash
# Generate PED data
./GenerateRandomData data/chr PED

# Run phasing program with PED input
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out -InputFormat PED
```

### Mixed formats
```bash
# Generate HAP data
./GenerateRandomData data/chr HAP

# Run phasing program with HAP input, PED output (default)
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out

# Or explicitly specify formats
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out -InputFormat HAP -OutputFormat PED
```

## Parameters

### Command Line Arguments

1. **`output_path`** (required):
   - Base path for output files
   - Example: `data/chr` or `test_data/chr`
   - Files will be named: `<output_path>1.hap`, `<output_path>2.hap`, etc.

2. **`format`** (optional):
   - Output format: `HAP` (default) or `PED`
   - Case-insensitive
   - Example: `PED` or `ped`

### Internal Constants

Defined in `GenerateRandomData.cpp`:
- `NUM_INDIVIDUALS`: 1000 (number of individuals)
- `NUM_SNPS_PER_CHR`: 1000 (SNPs per chromosome)
- `NUM_CHROMOSOMES`: 22 (chromosomes 1-22)
- `MIN_ALLELE_FREQ`: 0.05 (minimum allele frequency for alternate allele)
- `MAX_ALLELE_FREQ`: 0.95 (maximum allele frequency for alternate allele)

**Note on Allele Frequencies**: The generator directly creates Minor Allele Frequency (MAF) values in the range [0.05, 0.5]. This ensures that **all SNPs have a MAF of at least 5%**, which is important for statistical power in phasing algorithms. The MAF is generated uniformly in this range, avoiding very rare variants (MAF < 0.05) that can be problematic for phasing analysis.

To change these values, edit the constants in the source file and recompile.

## Technical Details

### Random Number Generation

- **Seed**: Current time (ensures different data each run)
- **Generator**: Standard C `rand()` function
- **Distribution**: Uniform for MAF, Hardy-Weinberg for genotypes

### Allele Frequency Distribution

- **MAF range**: 0.05 to 0.5 (uniform random, directly generated)
- **Distribution**: Uniform random for MAF values
- **Guarantee**: All SNPs have a MAF of at least 5% (0.05), ensuring no SNPs have a MAF lower than 5%
- **Rationale**: This range avoids very rare variants (MAF < 0.05) which can be problematic for phasing algorithms, while still covering a wide spectrum of allele frequencies including both common and less common variants. By generating MAF directly, we ensure that every SNP meets the minimum 5% threshold for reliable phasing analysis.
- **Usage**: Same MAF values used for:
  - Genotype generation (HAP/PED files)
  - MAF file output

### Genotype Generation

Genotypes follow **Hardy-Weinberg equilibrium**:

For a SNP with MAF = q and major allele frequency p = 1 - q:

- **P(homozygous reference 0/0)**: p²
- **P(heterozygous 0/1)**: 2pq
- **P(homozygous alternate 1/1)**: q²
- **P(missing)**: Small probability (encoded as random assignment)

### SNP Position Generation

- **Distribution**: Evenly distributed across chromosome length
- **Variation**: Small random variation (±10,000 bp)
- **Chromosome lengths**: Based on approximate human chromosome sizes

### Data Consistency

- **Same MAF**: All files (HAP/PED and MAF) use the same MAF values
- **Consistent genotypes**: Genotypes match the specified MAF
- **Reproducibility**: Same seed produces same data (if you set seed manually)

## Notes

- The random number generator is seeded with the current time, so each run produces different data
- To get reproducible results, modify the code to use a fixed seed: `srand(12345);`
- MAF values (Minor Allele Frequency) are directly generated uniformly in the range [0.05, 0.5]
- **All SNPs are guaranteed to have a MAF of at least 5%** (0.05), ensuring no SNPs have a MAF lower than 5%
- This avoids very rare variants (MAF < 0.05) which are often filtered out in real analyses and ensures sufficient statistical power for phasing algorithms
- SNP positions are distributed evenly across each chromosome with small random variations
- Missing data is rare but can occur (encoded as genotype 3, then randomly assigned)
- All individuals are independent (no family structure in generated data)
- No linkage disequilibrium is modeled (SNPs are independent)

## Troubleshooting

### Common Issues

#### Error: "Could not open file for writing"

**Causes and solutions**:
1. **Directory doesn't exist**:
   ```bash
   # Create directory first
   mkdir -p data
   ./GenerateRandomData data/chr
   ```

2. **Permission denied**:
   ```bash
   # Check permissions
   ls -ld data/
   # Fix if needed
   chmod 755 data/
   ```

3. **Invalid path**:
   - Use relative paths: `data/chr`
   - Or absolute paths: `/home/user/data/chr`
   - Avoid special characters in paths

#### Files are empty or incomplete

**Causes and solutions**:
1. **Insufficient disk space**:
   ```bash
   # Check available space
   df -h .
   # Free up space if needed
   ```

2. **Interrupted execution**:
   - Re-run the generator
   - Check for error messages during generation

3. **Path issues**:
   - Verify output path is correct
   - Check that path doesn't contain invalid characters

#### Compilation errors

**Common issues**:
1. **Missing compiler**:
   ```bash
   # Install g++
   sudo apt-get install g++  # Ubuntu/Debian
   ```

2. **OpenMP not found**:
   - Usually not required for generator
   - Can compile without OpenMP if needed

#### Generated data seems unrealistic

**Notes**:
- Data is synthetic and simplified
- No linkage disequilibrium is modeled
- No population structure is included
- For more realistic data, use real datasets or specialized simulators

### Performance

**Generation time** (approximate):
- 1000 individuals, 1000 SNPs/chr: ~1-5 minutes
- Depends on disk I/O speed

**File sizes**:
- HAP files: ~8-10 MB per chromosome
- PED files: ~8-10 MB per chromosome
- MAF files: ~50-100 KB per chromosome
- **Total**: ~200-250 MB for all chromosomes

### Validation

After generation, verify files:

```bash
# Check file existence
ls -lh data/chr*.hap  # or .ped

# Check file sizes (should be similar)
ls -lh data/chr*.hap | head -5

# Check file format (first few lines)
head -3 data/chr1.hap
head -3 data/chr1.ped
head -5 data/chr1.maf
```

