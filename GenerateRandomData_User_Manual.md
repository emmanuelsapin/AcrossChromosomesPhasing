# GenerateRandomData — User Manual

## Overview

**GenerateRandomData** is a synthetic data generator that creates HAP, PED, and MAF files for testing the across-chromosome phasing program (acrossCHRphasing). It produces genotype data with a predefined trio structure: one offspring and two parents.

---

## 1. Building

### Prerequisites

- **GCC** (g++) with C++ support

### Build command

```bash
g++ -o GenerateRandomData GenerateRandomData.cpp -O2
```

### Windows

```bat
g++ -o GenerateRandomData.exe GenerateRandomData.cpp -O2
```

---

## 2. Usage

```bash
./GenerateRandomData <output_path> [format]
```

| Argument | Description |
|----------|-------------|
| `output_path` | Base path for output files (e.g. `data/chr`). **Required.** |
| `format` | Optional: `HAP`, `PED`, or `BOTH` (default). HAP and PED contain the same data when `BOTH` is used. |

---

## 3. What It Generates

| Output | Description |
|--------|--------------|
| **HAP files** | `{output_path}1.hap`, `{output_path}2.hap`, … `{output_path}22.hap` |
| **MAF files** | `{output_path}1.maf`, … `{output_path}22.maf` |
| **PED files** | `{output_path}1.ped`, … `{output_path}22.ped` (when format is PED or BOTH) |
| **Window file** | `recombinaisonWindows.txt` in the same directory as the HAP files |

The window file defines recombination segment boundaries per chromosome and is required by acrossCHRphasing.

---

## 4. Data Layout

| Parameter | Value |
|-----------|-------|
| Individuals | 1000 |
| SNPs per chromosome | 1000 |
| Chromosomes | 1–22 |
| Individual 0 | Offspring |
| Individual 1 | Parent 1 |
| Individual 2 | Parent 2 |
| Recombination | One crossover event per chromosome at a random SNP |
| Minor allele frequency (MAF) | Between 0.05 and 0.95 |

---

## 5. Output File Formats

### HAP file format

One file per chromosome: `{output_path}{chr}.haps` (e.g. `data/chr1.haps`).

**Format**: One line per SNP, space/tab separated:

```
chr rsid position ref alt allele1_0 allele2_0 allele1_1 allele2_1 ...
```

| Field | Description |
|-------|-------------|
| `chr` | Chromosome number (1–22) |
| `rsid` | SNP identifier (e.g. rs1_0) |
| `position` | Position in base pairs |
| `ref` | Reference allele |
| `alt` | Alternate allele |
| `allele1_i allele2_i` | Allele pair for individual *i* (0 or 1) |

### MAF file format

One file per chromosome: `{output_path}{chr}.maf`.

Contains SNP_ID, position, and minor allele frequency per SNP.

### PED file format

One file per chromosome: `{output_path}{chr}.ped`.

Pedigree format with family/individual IDs and genotype columns. Individual 0 (offspring) has Parent 1 and Parent 2 as father/mother.

### Window file format (recombinaisonWindows.txt)

One line per chromosome (0–22):

```
chr bound0 bound1 bound2 bound3 bound4
```

- `bound` values are 0-based SNP indices.
- `-1` means end of chromosome.
- Default layout: windows 0–500 and 500–1000 SNPs per chromosome.

---

## 6. Examples

### Generate HAP and PED files (default)

```bash
./GenerateRandomData data/chr BOTH
```

### Generate HAP files only

```bash
./GenerateRandomData data/chr HAP
```

### Generate PED files only

```bash
./GenerateRandomData data/chr PED
```

---

## 7. Troubleshooting

| Issue | Action |
|-------|--------|
| `ERROR: Cannot create ...` | Check write permissions and disk space |
| `ERROR: Could not open file for writing` | Ensure the output directory exists or can be created |
| Unknown format warning | Use `HAP`, `PED`, or `BOTH` (case-insensitive) |

---

*Last updated: February 2025.*
