# acrossCHRphasing — User Manual

## Overview

**acrossCHRphasing** performs **across-chromosome haplotype phasing** — determining which allele at each SNP came from which parent (parent-of-origin phasing). It is based on the method described in Sapin, Kelly, Keller (2020).

The program reads genotype data from HAP files and outputs phased haplotypes (alleles assigned to Parent 1 vs Parent 2) for a focal individual.

---

## 1. Building

### Prerequisites

- **GCC** (g++) with C++ support
- **OpenMP** support (included with GCC on Linux; MinGW on Windows)

### Build (Windows)

```bat
build_acrossCHRphasing.bat
```

Or manually:

```bat
g++ -o acrossCHRphasing.exe acrossCHRphasing.cpp -O2 -fopenmp -lm -Wl,--stack,67108864
```

### Build (Linux)

```bash
g++ -o acrossCHRphasing acrossCHRphasing.cpp -O2 -fopenmp -lm "-Wl,--stack,67108864"
```

The `-Wl,--stack,67108864` flag sets a 64 MB stack size to avoid stack overflow.

---

## 2. Usage

```bash
./acrossCHRphasing pathHap=<path> pathwindow=<path> posOffspring=<N> [posParent1=<N>] [posParent2=<N>] pathresult=<path>
```

Arguments are passed as `key=value` and may be in any order.

### Required arguments

| Argument | Description |
|----------|-------------|
| `pathHap=` | Base path for HAP files. Program expects `{pathHap}1.hap`, `{pathHap}2.hap`, … `{pathHap}22.hap` |
| `pathwindow=` | Full path to the window coordinates file (e.g. `data/recombinaisonWindows.txt`) |
| `posOffspring=` | 0-based index of the individual to phase in the HAP file |
| `pathresult=` | Full path to the output file (will be overwritten) |

### Optional arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `posParent1=` | -1 | 0-based index of Parent 1 in the HAP file. Use -1 when no Parent 1 |
| `posParent2=` | -1 | 0-based index of Parent 2 in the HAP file. Use -1 when no Parent 2 |

When a parent is missing (-1), the program uses the offspring ID instead (effectively no trio information for that parent).

### Help

```bash
./acrossCHRphasing --help
# or
./acrossCHRphasing -h
# or
./acrossCHRphasing help
```

---

## 3. Input File Formats

### HAP file format

One file per chromosome: `{pathHap}{chr}.hap` (e.g. `data/chr1.hap`).

**Format**: One line per SNP, space/tab separated:

```
chr rsid position ref alt allele1_0 allele2_0 allele1_1 allele2_1 ...
```

| Field | Description |
|-------|-------------|
| `chr` | Chromosome number (1–22) |
| `rsid` | SNP identifier (e.g. rs123) |
| `position` | Position in base pairs |
| `ref` | Reference allele (0) |
| `alt` | Alternate allele (1) |
| `allele1_i allele2_i` | Allele pair for individual *i* (0 or 1) |

**Example line** (3 individuals):

```
1 rs1_0 2489564 A G 0 1 1 0 0 0
```

Genotype encoding used internally: 0 = hom ref, 1/2 = het, 3 = hom alt.

### Window file format (recombinaisonWindows.txt)

One line per chromosome (0–22). Each line defines segment boundaries:

```
chr bound0 bound1 bound2 bound3 bound4
```

- `bound` values are 0-based SNP indices.
- `-1` means end of chromosome.
- The program constructs non-overlapping windows (e.g. 0–500, 500–1000 per chromosome).

**Example** (for 1000 SNPs per chr):

```
# chr bound1 bound2 bound3 bound4 bound5
0 0 0 0 0
1 500 -1 0 0 0
2 500 -1 0 0 0
...
22 500 -1 0 0 0
```

### Hotspot file format (pathHotspots)

When using hotspot-based segments, provide a file with columns: Chr, cM, cMperMb, bp. **pathMap** is also required (for SNP indexing). The program does not derive hotspots from the map file; you must supply the hotspot file.

---

## 4. Output Format

The program writes one output file at `pathresult`.

**Header**:

```
#chr snp allelevalue
```

**Body**: One line per SNP:

```
chr snp allelevalue
```

| Field | Description |
|-------|-------------|
| `chr` | Chromosome (1–22) |
| `snp` | 0-based SNP index |
| `allelevalue` | Phased genotype: 0 = hom ref, 1 = het (alt on haplotype 1), 2 = het (alt on haplotype 2), 3 = hom alt |

Haplotype 1 is inferred as originating from Parent 1; haplotype 2 from Parent 2. At heterozygous sites, `1` and `2` indicate which haplotype carries the alternate allele.

---

## 5. Algorithm Overview

1. **Load genomes**: Reads HAP files for chromosomes 1–22, packs genotypes into a binary format.
2. **Allele frequencies**: Computes minor allele frequency (MAF) per SNP.
3. **Pi-hat calculation**: Computes pi-hat between the focal individual and all others to identify related individuals (pi-hat > 0.33 excluded).
4. **Across-chromosome phasing**: Uses segment correlations and window boundaries to assign alleles to parent-of-origin haplotypes.
5. **Output**: Writes phased genotypes to the output file.

---

## 6. System Requirements and Limits

### Memory

- **~32 GB** for the `pihatagainstallchrMPphaseerror` array
- Additional memory for genomes and other structures
- Runs best with **≥64 GB RAM** for typical datasets

**Note**: Run on a compute node with sufficient memory, not on a shared login node.

### Limits (compile-time constants)

| Constant | Value | Description |
|----------|-------|-------------|
| MAXPOP | 435188 | Max individuals |
| MAXSNP | 330005 | Max SNPs total |
| NSNPPERCHR | 26230 | Max SNPs per chromosome |
| NCHR | 23 | Chromosomes (1–22 used) |

### Parallelization

Uses OpenMP. Set thread count via `OMP_NUM_THREADS`:

```bash
export OMP_NUM_THREADS=8   # Linux
set OMP_NUM_THREADS=8      # Windows CMD
```

---

## 7. Example Workflow

Using data from **GenerateRandomData** (see `GenerateRandomData_User_Manual.md`):

```bash
# Step 1: Generate test data
./GenerateRandomData data/chr BOTH

# Step 2: Run phasing (individual 0 = offspring, 1 = Parent 1, 2 = Parent 2)
./acrossCHRphasing pathHap=data/chr pathwindow=data/recombinaisonWindows.txt posOffspring=0 posParent1=1 posParent2=2 pathresult=./phased_genome.txt

# Step 3: Inspect output
head -20 phased_genome.txt
```

---

## 8. Troubleshooting

| Issue | Possible cause | Action |
|-------|----------------|--------|
| `Cannot open file ... Check pathHap` | Wrong HAP path or missing files | Ensure `{pathHap}1.hap` … `{pathHap}22.hap` exist |
| `Cannot open ... Check pathwindow` | Wrong or missing window file | Ensure `pathwindow` points to a valid file |
| `RETURN(1): pathHap=, pathwindow=, and posOffspring= are required` | Missing required arguments | Set `pathHap`, `pathwindow`, and `posOffspring` |
| Stack overflow / crash before “Start main phasing program” | Stack too small | Rebuild with `-Wl,--stack,67108864` |
| `Failed to allocate ... (~32 GB)` | Not enough RAM | Use a compute node with ≥64 GB RAM |
| `invalid n_individuals` | HAP format error | Check that each line has chr, rsid, pos, ref, alt, plus 2×*N* allele values |

---

## 9. Reference

Sapin, E., Kelly, J.K., Keller, M.P. (2020). Across-chromosome phasing of haplotypes in a pedigreed population.

---

*Last updated: February 2025.*
