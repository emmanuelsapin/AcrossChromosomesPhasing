# PED File Format Support

## Overview

The phasing program now supports both HAP and PED file formats for input and output.

## File Formats

### HAP Format (Input)
- **Extension**: `.hap`
- **Structure**: One SNP per line
- **Format**: `rsID position allele1_indiv1 allele2_indiv1 allele1_indiv2 allele2_indiv2 ...`
- **Example**: `rs1_0 1234567 0 0 0 1 1 1 ...`

### PED Format (Input/Output)
- **Extension**: `.ped`
- **Structure**: One individual per line
- **Format**: `Family_ID Individual_ID Father_ID Mother_ID Sex Phenotype Allele1_SNP1 Allele2_SNP1 Allele1_SNP2 Allele2_SNP2 ...`
- **Example**: `0 0 0 0 -9 0 0 0 0 1 1 1 ...`

## Command Line Options

### Input Format
```
-InputFormat HAP|PED
```
- **Default**: HAP
- Specifies the format of input files
- Options: `HAP` or `PED` (case-insensitive)

### Output Format
```
-OutputFormat HAP|PED
```
- **Default**: PED
- Specifies the format of output files
- Options: `HAP` or `PED` (case-insensitive)

## Usage Examples

### Example 1: HAP Input, PED Output (Default)
```bash
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out
```
This reads HAP files (`data/chr1.hap`, `data/chr2.hap`, etc.) and writes PED files (`results/out1.ped`, `results/out2.ped`, etc.).

### Example 2: PED Input, PED Output
```bash
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out -InputFormat PED
```
This reads PED files (`data/chr1.ped`, `data/chr2.ped`, etc.) and writes PED files.

### Example 3: HAP Input, HAP Output
```bash
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out -OutputFormat HAP
```
This reads HAP files and writes HAP files (note: HAP output is currently implemented as PED output).

### Example 4: PED Input, HAP Output
```bash
./ProgramPhasing -NbIndiv 1000 -PathInput data/chr -PathOutput results/out -InputFormat PED -OutputFormat HAP
```
This reads PED files and writes HAP files (note: HAP output is currently implemented as PED output).

## PED File Format Details

### Header Columns (First 6 columns)
1. **Family_ID**: Family identifier (can be 0 or any integer)
2. **Individual_ID**: Individual identifier (must be unique)
3. **Father_ID**: Father's Individual_ID (0 if unknown)
4. **Mother_ID**: Mother's Individual_ID (0 if unknown)
5. **Sex**: Sex code (1=male, 2=female, -9=unknown)
6. **Phenotype**: Phenotype value (0=unknown, -9=missing)

### Genotype Data
- Each subsequent pair of columns represents one SNP
- First allele of the pair is the first allele
- Second allele of the pair is the second allele
- Alleles are encoded as: `0` (reference), `1` (alternate), or missing markers (`-`, `N`)

### Example PED File Line
```
0 0 0 0 -9 0 0 0 0 1 1 1 0 1 1 0
```
This represents:
- Individual 0, no family/parent info
- SNP 1: 0/0 (homozygous reference)
- SNP 2: 0/1 (heterozygous)
- SNP 3: 1/1 (homozygous alternate)
- SNP 4: 0/1 (heterozygous)

## Implementation Notes

- PED files are read chromosome by chromosome
- The program automatically detects the number of SNPs from the first line
- Missing data in PED files is treated as 0 (reference allele)
- The program supports large PED files with up to 1,000,000 characters per line
- Progress indicators are shown every 100 individuals

## File Naming Convention

### Input Files
- **HAP format**: `<PathInput><chromosome>.hap` (e.g., `data/chr1.hap`)
- **PED format**: `<PathInput><chromosome>.ped` (e.g., `data/chr1.ped`)

### Output Files
- **PED format**: `<PathOutput><chromosome>.ped` (e.g., `results/out1.ped`)
- **HAP format**: Currently writes as PED format

## Compatibility

- The program maintains backward compatibility with HAP input format
- Default behavior (HAP input, PED output) matches previous versions
- All existing command-line options continue to work

