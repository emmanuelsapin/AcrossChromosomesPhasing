# Haplotype Phasing Engine - Professional Architecture

## Overview

This project implements a genetic haplotype phasing engine using a professional object-oriented architecture with advanced design patterns.

## Architecture

### Namespaces

- **`PhasingEngine`**: Main namespace containing all classes
- **`PhasingEngine::Constants`**: System constants
- **`PhasingEngine::ErrorCodes`**: Standardized error codes

### Abstract Interfaces

1. **`IGenomeDataAccessor`**: Interface for genome data access
2. **`IRelativeFinder`**: Interface for relative identification
3. **`IPhasingStrategy`**: Strategy pattern for different phasing algorithms

### Main Classes

#### `GenomeDataManager`
- Complete genome data management
- Index validation (chromosome, individual, SNP)
- Automatic MAF (Minor Allele Frequencies) calculation
- Robust memory management with destructors

#### `RelativeIdentificationEngine`
- PIHAT matrix computation
- Best relative identification
- Sorting and ranking of relationships

#### `PhasingAlgorithmEngine`
- Main phasing execution engine
- Multiple strategy support (Strategy Pattern)
- Chromosomal window merging
- Phasing correction application

#### `ConfigurationManager`
- Centralized configuration management
- Command-line argument parsing
- Parameter validation

#### `HaplotypePhasingProgram`
- Main class orchestrating the pipeline
- Complete lifecycle management
- Execution statistics

### Error Handling

- **`PhasingException`**: Custom exception with error codes
- Systematic input validation
- Descriptive error messages

### Implemented Design Patterns

1. **Factory Pattern**: `GenomeFileLoader` for different formats
2. **Strategy Pattern**: `IPhasingStrategy` for multiple algorithms
3. **RAII**: Automatic resource management
4. **Interface Segregation**: Specialized interfaces

## Compilation

```bash
make -f Makefile_oo
```

## Usage

```bash
./ProgramPhasing_OO -NbIndiv <N> -PathInput <prefix> -PathOutput <prefix> [-PathMAF <path>] [-Verbose <0|1>]
```

## Documentation

Doxygen documentation is integrated in header files. Generate with:

```bash
doxygen Doxyfile
```
