# Haplotype Phasing Engine

A professional haplotype phasing engine for genetic data analysis using relative-based inference algorithms. A video-based visualization illustrating a single execution of the algorithm is available, along with the source code used to generate the video. This project provides a comprehensive object-oriented framework for processing genomic data and identifying haplotype phases. 

## Features

- **Modular Architecture**: Clean separation of concerns with independent modules
- **Professional Design**: Object-oriented design with design patterns (Factory, Strategy, RAII)
- **Error Handling**: Comprehensive exception handling with custom error codes
- **Parallel Processing**: OpenMP support for multi-threaded computations
- **Multiple File Formats**: Support for HAP, PED, and VCF formats
- **PIHAT Computation**: Advanced relative identification using PIHAT matrices
- **Configuration Management**: Flexible command-line configuration
- **Data Generator**: Random data generator for testing (1000 individuals, 1000 SNPs per chromosome)
- **Video Presentation**: Possible executions of the algorithm are depicted in brief video demonstrations, illustrating the interactions among the various constituent elements.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Compilation](#compilation)
- [Usage](#usage)
- [Architecture](#architecture)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Prerequisites

- **C++ Compiler**: GCC 4.9+ or Clang 3.5+ with C++11 support
- **OpenMP**: For parallel processing support
- **Make**: Build automation tool
- **Linux/Unix**: Windows support via WSL or MinGW

### System Requirements

- Minimum 4GB RAM (8GB+ recommended for large datasets)
- Disk space: 500MB+ for installation and temporary files
- Multi-core processor recommended for optimal performance

## Quick Start

1. **Clone or download the repository**
   ```bash
   cd Fileforgithub
   ```

2. **Compile the project**
   ```bash
   make -f Makefile_modular
   ```

3. **Generate test data (optional)**
   ```bash
   # Using PowerShell script (Windows)
   powershell -ExecutionPolicy Bypass -File .\GenererDonnees.ps1
   
   # Or compile and run the C++ generator
   g++ -std=c++11 -O2 -o GenerateurDonneesAleatoires GenerateurDonneesAleatoires.cpp
   ./GenerateurDonneesAleatoires ./donnees_aleatoires/
   ```

4. **Run with sample data**
   ```bash
   ./ProgramPhasing_Modular -NbIndiv 1000 -PathInput ./donnees_aleatoires/ -PathOutput ./results/ -Verbose 1
   ```

## Compilation

### Modular Build (Recommended)

```bash
make -f Makefile_modular
```

This creates `ProgramPhasing_Modular` executable.

### Object-Oriented Build

```bash
make -f Makefile_oo
```

This creates `ProgramPhasing_OO` executable.

### Clean Build

```bash
make -f Makefile_modular clean
```

### Compilation Options

Edit the Makefile to customize:
- `CXXFLAGS`: Compiler flags (optimization, warnings)
- `LDFLAGS`: Linker flags
- `CXX`: Compiler selection

## Usage

### Basic Command Line

```bash
./ProgramPhasing_Modular [OPTIONS]
```

### Required Arguments

- `-NbIndiv <N>`: Number of individuals to process
- `-PathInput <path>`: Input file path prefix (without chromosome number)
- `-PathOutput <path>`: Output file path prefix

### Optional Arguments

- `-PathMAF <path>`: Path to MAF (Minor Allele Frequency) file
- `-Verbose <0|1>`: Enable verbose output (default: 1)

### Example Commands

**Basic usage:**
```bash
./ProgramPhasing_Modular -NbIndiv 5000 -PathInput /data/genomes/ -PathOutput /results/phased/
```

**With MAF file:**
```bash
./ProgramPhasing_Modular -NbIndiv 10000 -PathInput ./input/ -PathOutput ./output/ -PathMAF ./MAF.txt -Verbose 1
```

**Quiet mode:**
```bash
./ProgramPhasing_Modular -NbIndiv 2000 -PathInput ./data/ -PathOutput ./results/ -Verbose 0
```

### Input File Format

The program expects input files in HAP format:
- Files named: `<PathInput><chromosome>.hap`
- Example: If `PathInput` is `./data/`, files should be `./data/1.hap`, `./data/2.hap`, etc.

### Output Format

Output files are generated in PED format:
- Files named: `<PathOutput><chromosome>.ped`
- Each file contains phased genotypes for all individuals

## Architecture

### Module Structure

```
include/
├── Constants.h              # System constants
├── ErrorCodes.h             # Error code definitions
├── Exceptions.h             # Custom exception classes
├── Interfaces.h             # Abstract interfaces
├── ChromosomeDivider.h      # Chromosomal segment structure
├── GenomeDataManager.h      # Genome data management
├── RelativeIdentificationEngine.h  # Relative identification
├── GenomeFileLoader.h       # File loading utilities
├── PhasingAlgorithmEngine.h # Core phasing algorithms
├── OutputFileWriter.h       # Output generation
├── ConfigurationManager.h   # Configuration management
└── HaplotypePhasingProgram.h # Main application class

src/
├── GenomeDataManager.cpp
├── RelativeIdentificationEngine.cpp
├── GenomeFileLoader.cpp
├── ConfigurationManager.cpp
├── OutputFileWriter.cpp
├── PhasingAlgorithmEngine.cpp
└── HaplotypePhasingProgram.cpp
```

### Key Components

- **GenomeDataManager**: Manages genomic data with validation and caching
- **RelativeIdentificationEngine**: Computes PIHAT matrices and identifies relatives
- **PhasingAlgorithmEngine**: Executes phasing algorithms with multiple strategies
- **ConfigurationManager**: Handles command-line arguments and configuration
- **OutputFileWriter**: Generates output files in various formats

## Documentation

### API Documentation

Detailed API documentation is available in header files with Doxygen comments. Generate documentation with:

```bash
doxygen Doxyfile
```

### Additional Documentation

- `README_ARCHITECTURE.md`: Detailed architecture overview
- `README_MODULAR.md`: Modular structure documentation
- `STRUCTURE_MODULAIRE.md`: Module dependency structure
- `STRUCTURE_DEPENDANCES.md`: Dependency organization
- `USER_MANUAL.md`: Complete user manual
- `README_Generateur.md`: Data generator documentation
- `GUIDE_GITHUB.md`: Guide for transferring project to GitHub

## Troubleshooting

### Common Issues

**Compilation Errors:**
- Ensure C++11 support is enabled
- Check OpenMP installation: `gcc -fopenmp --version`
- Verify all dependencies are present

**Runtime Errors:**
- Check file paths are correct and accessible
- Verify input file format matches expected HAP format
- Ensure sufficient memory for dataset size

**Performance Issues:**
- Use `-Verbose 0` to reduce output overhead
- Ensure OpenMP is properly configured
- Check system resources (CPU, memory)

### Error Codes

- `FILE_NOT_FOUND`: Input file cannot be accessed
- `INVALID_INPUT`: Invalid command-line arguments
- `MEMORY_ALLOCATION_FAILED`: Insufficient memory
- `INVALID_CHROMOSOME`: Invalid chromosome index
- `INVALID_INDIVIDUAL`: Invalid individual index
- `INVALID_SNP_INDEX`: Invalid SNP index
- `DATA_CORRUPTION`: Data integrity issue

## Performance

### Optimization Tips

1. **Parallel Processing**: Ensure OpenMP threads are configured
   ```bash
   export OMP_NUM_THREADS=8
   ```

2. **Memory Management**: For large datasets, process in batches

3. **File I/O**: Use fast storage (SSD) for input/output files

### Benchmarks

Typical performance on modern hardware:
- 1,000 individuals: ~5-10 minutes
- 10,000 individuals: ~30-60 minutes
- 100,000 individuals: ~4-8 hours

*Performance varies based on dataset characteristics and hardware*

## Contributing

Contributions are welcome! Please follow these guidelines:

1. Maintain code style and documentation standards
2. Add tests for new features
3. Update documentation as needed
4. Follow the modular architecture principles

## License

[Specify your license here]

## Contact

For questions, issues, or contributions, please [specify contact method].

## Acknowledgments

This project implements advanced haplotype phasing algorithms for genetic data analysis.

---

**Version**: 2.0  
**Last Updated**: 2024


