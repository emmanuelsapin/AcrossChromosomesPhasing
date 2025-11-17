# Across Chromosomes Phasing

A comprehensive haplotype phasing program that performs phasing across multiple chromosomes using advanced algorithms. The program supports both HAP and PED file formats for input and output.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [File Formats](#file-formats)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Overview

This program implements a sophisticated haplotype phasing algorithm that can:
- Phase haplotypes across multiple chromosomes (1-22)
- Process individuals with or without known parent information
- Support selective processing of specific individuals
- Handle both HAP and PED file formats
- Generate random test data for validation

The program uses object-oriented design with modular libraries for file I/O and phasing algorithms.

## Features

### Core Functionality
- **Cross-chromosome phasing**: Analyzes relationships between chromosomes to improve phasing accuracy
- **Parent-aware phasing**: Uses parent information when available for more accurate results
- **Selective processing**: Process only specified individuals from a list
- **Multiple file formats**: Support for HAP and PED input/output formats
- **Parallel processing**: Uses OpenMP for multi-threaded computation

### Data Generation
- **Random data generator**: Creates synthetic genetic data for testing
- **Configurable parameters**: 1000 individuals, 1000 SNPs per chromosome (default)
- **Realistic distributions**: Uses Hardy-Weinberg equilibrium for genotype generation
- **Multiple formats**: Generates HAP or PED format files

## Project Structure

```
.
├── main.cpp                 # Main entry point
├── phasing_program.h        # Main program class declaration
├── phasing_program.cpp      # Main program implementation
├── types.h                  # Global type definitions and structures
├── GenerateRandomData.cpp   # Random data generator
├── Makefile                 # Build configuration
│
├── libfileio/               # File I/O library
│   ├── fileio.h            # File I/O function declarations
│   ├── fileio.cpp          # File I/O implementations (HAP, PED, etc.)
│   └── Makefile            # Library build configuration
│
├── libphasing/              # Phasing algorithms library
│   ├── phasing.h           # Phasing function declarations
│   ├── phasing.cpp         # Phasing algorithm implementations
│   └── Makefile            # Library build configuration
│
├── libtypes/                # Types library (optional)
│   └── types.h             # Type definitions (currently in root)
│
├── README.md                # This file
├── USER_MANUAL.md           # Comprehensive user manual
├── README_GENERATOR.md      # Data generator documentation
└── README_PED_FORMAT.md     # PED format documentation
```

## Installation

### Prerequisites

- **C++ Compiler**: GCC 4.9+ or compatible (supports C++11)
- **OpenMP**: For parallel processing
- **Make**: For building the project
- **Math library**: Standard math library (linked automatically)

### Compilation

#### Compile Everything
```bash
make all
```

This will:
1. Build the `libfileio` library
2. Build the `libphasing` library
3. Compile the main program (`ProgramPhasing`)
4. Compile the data generator (`GenerateRandomData`)

#### Compile Individual Components

**Main program only:**
```bash
make ProgramPhasing
```

**Data generator only:**
```bash
make generator
```

**Libraries only:**
```bash
make libs
```

#### Clean Build Artifacts
```bash
make clean
```

### Windows Notes

On Windows, you may need:
- MinGW-w64 or MSYS2 for `make` and `g++`
- Or use the provided PowerShell scripts (if available)

## Quick Start

### 1. Generate Test Data

```bash
# Generate HAP format data (1000 individuals, 1000 SNPs/chr)
./GenerateRandomData data/chr HAP

# Or generate PED format data
./GenerateRandomData data/chr PED
```

### 2. Run Phasing

```bash
# Basic usage with HAP input (default)
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out

# With PED input
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out -InputFormat PED
```

### 3. View Results

Results will be written to:
- `results/out1.ped`, `results/out2.ped`, ..., `results/out22.ped`

## Documentation

- **[USER_MANUAL.md](USER_MANUAL.md)**: Comprehensive user manual with detailed explanations
- **[README_GENERATOR.md](README_GENERATOR.md)**: Data generator documentation
- **[README_PED_FORMAT.md](README_PED_FORMAT.md)**: PED file format specifications

## File Formats

### HAP Format
- **Extension**: `.hap`
- **Structure**: One SNP per line
- **Format**: `rsID position allele1_indiv1 allele2_indiv1 ...`
- **Use case**: Efficient storage, SNP-centric view

### PED Format
- **Extension**: `.ped`
- **Structure**: One individual per line
- **Format**: `Family_ID Individual_ID Father_ID Mother_ID Sex Phenotype Allele1_SNP1 Allele2_SNP1 ...`
- **Use case**: Standard format, individual-centric view, compatible with PLINK

### MAF Format
- **Extension**: `.maf`
- **Structure**: Tab-separated values
- **Format**: `SNP_ID\tPosition\tMAF`
- **Use case**: Minor Allele Frequency information

## Examples

### Example 1: Basic Phasing with HAP Files

```bash
# Generate data
./GenerateRandomData test_data/chr HAP

# Run phasing
./ProgramPhasing -NbIndiv 1000 -PathInput test_data/chr -PathOutput phased/out
```

### Example 2: Phasing with Parent Information

```bash
# Create parent info file (example_parents.txt)
# Format: Individual_ID Parent1_ID Parent2_ID
# 0 5 10
# 1 -1 -1
# 2 15 20

# Run with parent info
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out \
    -PathParentInfo example_parents.txt
```

### Example 3: Selective Individual Processing

```bash
# Create individual list (list_indiv.txt)
# Format: One individual ID per line
# 0
# 5
# 10
# 25

# Process only listed individuals
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out \
    -ListIndiv list_indiv.txt
```

### Example 4: PED Input and Output

```bash
# Generate PED data
./GenerateRandomData data/chr PED

# Process PED files
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out \
    -InputFormat PED -OutputFormat PED
```

### Example 5: Mixed Formats

```bash
# Read HAP, write PED (default behavior)
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out \
    -InputFormat HAP -OutputFormat PED
```

## Command Line Options

### Main Program Options

| Option | Description | Required | Default |
|--------|-------------|----------|---------|
| `-NbIndiv` | Number of individuals | Yes | - |
| `-PathInput` | Base path for input files | Yes | - |
| `-PathOutput` | Base path for output files | Yes | - |
| `-PathParentInfo` | Path to parent information file | No | - |
| `-ListIndiv` | Path to list of individuals to process | No | All individuals |
| `-InputFormat` | Input format: HAP or PED | No | HAP |
| `-OutputFormat` | Output format: HAP or PED | No | PED |

### Data Generator Options

```bash
./GenerateRandomData <output_path> [format]
```

- `output_path`: Base path for output files
- `format`: Optional - `HAP` (default) or `PED`

## Output Files

The program generates output files for each chromosome (1-22):
- **PED format**: `<PathOutput>1.ped`, `<PathOutput>2.ped`, ..., `<PathOutput>22.ped`
- **HAP format**: Currently writes as PED format

Each output file contains phased haplotypes for all processed individuals.

## Performance Considerations

- **Memory**: The program allocates memory based on `NbIndiv` and number of SNPs
- **Parallelization**: Uses OpenMP for parallel processing where applicable
- **File I/O**: Large files may take time to read/write
- **Processing time**: Depends on number of individuals, SNPs, and available CPU cores

## Troubleshooting

### Common Issues

**Error: "Number of individuals is zero or undefined"**
- Solution: Ensure `-NbIndiv` is specified and greater than 0

**Error: "Could not open file"**
- Solution: Check that input files exist and paths are correct
- For HAP: Files should be named `<PathInput>1.hap`, `<PathInput>2.hap`, etc.
- For PED: Files should be named `<PathInput>1.ped`, `<PathInput>2.ped`, etc.

**Error: "Failed to read parent info file"**
- Solution: Check file format (see USER_MANUAL.md for format specification)
- Ensure file exists and is readable

**Compilation errors on Windows**
- Solution: Install MinGW-w64 or use WSL (Windows Subsystem for Linux)
- Ensure OpenMP support is enabled in your compiler

**Memory errors with large datasets**
- Solution: Reduce number of individuals or SNPs per chromosome
- Check available system memory

## Contributing

This is a research/academic project. For contributions:
1. Follow the existing code style
2. Maintain backward compatibility
3. Update documentation for new features
4. Test with both HAP and PED formats

## License

[Specify your license here]

## Citation

If you use this software in your research, please cite:

[Add citation information]

## Contact

[Add contact information]

## Version History

- **Current Version**: [Version number]
- **Features**: 
  - Support for HAP and PED file formats
  - Parent-aware phasing
  - Selective individual processing
  - Modular library architecture
  - Random data generator
