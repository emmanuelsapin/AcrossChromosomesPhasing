/**
 * @file GenerateRandomData.cpp
 * @brief Random data generator for haplotype phasing program
 * 
 * This program generates random HAP, PED, and MAF files for testing the phasing algorithm.
 * It generates data for 1000 individuals and 1000 SNPs per chromosome (chromosomes 1-22).
 * Individual 0 is the offspring of individuals 1 and 2 (trio structure).
 * 
 * Usage:
 *   ./GenerateRandomData <output_path> [format]
 * 
 * Parameters:
 *   output_path: Base path for output files (e.g., "data/chr")
 *   format: Optional - "HAP", "PED", or "BOTH" (default) - HAP and PED represent same data
 * 
 * Output files (HAP and PED represent the same data when BOTH):
 *     - <output_path>1.haps, <output_path>2.haps, ..., <output_path>22.haps (HAP files)
 *     - <output_path>1.maf, <output_path>2.maf, ..., <output_path>22.maf (MAF files)
 *   PED format:
 *     - <output_path>1.ped, <output_path>2.ped, ..., <output_path>22.ped (PED files)
 *     - <output_path>1.maf, <output_path>2.maf, ..., <output_path>22.maf (MAF files)
 *   Window file (always): recombinaisonWindows.txt in same directory (windows 0-500, 500-1000 per chr)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#ifdef _WIN32
#include <direct.h>  // For _mkdir on Windows
#define mkdir(path, mode) _mkdir(path)
#else
#include <sys/stat.h>  // For mkdir on Linux/Unix
#include <sys/types.h>
#endif

// Constants for data generation
#define NUM_INDIVIDUALS 1000
#define NUM_SNPS_PER_CHR 1000
#define NUM_CHROMOSOMES 22
#define MIN_ALLELE_FREQ 0.05  // Minimum allele frequency (for alternate allele)
#define MAX_ALLELE_FREQ 0.95  // Maximum allele frequency (for alternate allele)
#define BOUNDS_BREAK_AT_SNP 500  // Segment break for recombinaisonWindows.txt (windows 0-500, 500-1000)

// Trio: individual 0 = offspring of individuals 1 and 2 (parents)
#define IDX_OFFSPRING 0
#define IDX_PARENT1   1
#define IDX_PARENT2   2

/**
 * @brief Generate a random double between 0 and 1
 * @return Random double value
 */
double random_double()
{
	return (double)rand() / (double)RAND_MAX;
}

/**
 * @brief Generate a random integer between min and max (inclusive)
 * @param min Minimum value
 * @param max Maximum value
 * @return Random integer
 */
int random_int(int min, int max)
{
	return min + rand() % (max - min + 1);
}

/**
 * @brief Generate a random SNP ID (rsID format)
 * @param chr Chromosome number
 * @param snp_index SNP index within chromosome
 * @param buffer Buffer to store the rsID (must be at least 20 characters)
 */
void generate_rsid(int chr, int snp_index, char* buffer)
{
	sprintf(buffer, "rs%d_%d", chr, snp_index);
}

/**
 * @brief Generate a random position on a chromosome
 * @param chr Chromosome number (1-22)
 * @param snp_index SNP index (0-based)
 * @return Random position in base pairs
 */
int generate_position(int chr, int snp_index)
{
	// Approximate chromosome lengths in base pairs
	int chr_lengths[23] = {0, 
		248956422, 242193529, 198295559, 190214555, 181538259,
		170805979, 159345973, 145138636, 138394717, 133797422,
		135086622, 133275309, 114364328, 107043718, 101991189,
		90338345, 83257441, 80373285, 58617616, 64444167,
		46709983, 50818468};
	
	// Distribute SNPs evenly across the chromosome
	int base_pos = (snp_index * chr_lengths[chr]) / NUM_SNPS_PER_CHR;
	// Add some random variation
	int variation = random_int(-10000, 10000);
	return base_pos + variation;
}

/**
 * @brief Generate a random allele frequency and compute MAF
 * @return Minor Allele Frequency (MAF) between 0.05 and 0.5 (inclusive)
 * 
 * This function generates a Minor Allele Frequency (MAF) that is guaranteed to be
 * at least 5% (0.05) and at most 50% (0.5). This ensures that all SNPs have
 * sufficient minor allele frequency for reliable phasing analysis.
 * 
 * The MAF is generated directly in the range [0.05, 0.5] to ensure no SNPs have
 * a MAF lower than 5%, which is important for statistical power in phasing algorithms.
 */
double generate_allele_frequency()
{
	// Generate MAF directly in the range [0.05, 0.5]
	// This ensures that no SNP has a MAF lower than 5%
	// MAF range: [MIN_ALLELE_FREQ, 0.5]
	double maf = MIN_ALLELE_FREQ + random_double() * (0.5 - MIN_ALLELE_FREQ);
	
	// Ensure MAF is at least MIN_ALLELE_FREQ (safety check)
	if (maf < MIN_ALLELE_FREQ)
	{
		maf = MIN_ALLELE_FREQ;
	}
	
	// Ensure MAF is at most 0.5 (safety check)
	if (maf > 0.5)
	{
		maf = 0.5;
	}
	
	return maf;
}

/**
 * @brief Generate a random genotype for an individual based on MAF
 * @param maf Minor Allele Frequency (between 0 and 0.5)
 * @return Genotype: 0 (homozygous reference), 1 (heterozygous), 2 (homozygous alternate), 3 (missing)
 */
int generate_genotype(double maf)
{
	double r = random_double();
	
	// Calculate genotype probabilities based on Hardy-Weinberg equilibrium
	// p = 1 - maf (major allele frequency)
	// q = maf (minor allele frequency)
	double p = 1.0 - maf;
	double q = maf;
	
	// P(homozygous reference) = p^2
	// P(heterozygous) = 2pq
	// P(homozygous alternate) = q^2
	double prob_hom_ref = p * p;
	double prob_het = 2.0 * p * q;
	double prob_hom_alt = q * q;
	
	if (r < prob_hom_ref)
	{
		return 0;  // Homozygous reference (0/0)
	}
	else if (r < prob_hom_ref + prob_het)
	{	if (rand()%2)
		return 1;  // Heterozygous (0/1)
		else
		return 2;  // Heterozygous (1/0)
	}
	else if (r < prob_hom_ref + prob_het + prob_hom_alt)
	{
		return 3;  // Homozygous alternate (1/1)
	}
	else
	{
		// Small probability of missing data
		return 4;  // Missing
	}
}

/**
 * @brief Convert genotype to allele pair for HAP file format
 * @param genotype Genotype value (0, 1, 2, or 3)
 * @param allele1 Output: first allele (0 or 1)
 * @param allele2 Output: second allele (0 or 1)
 */
void genotype_to_alleles(int genotype, int* allele1, int* allele2)
{
	switch (genotype)
	{
		case 0:  // Homozygous reference
			*allele1 = 0;
			*allele2 = 0;
			break;
		case 1:  // Heterozygous
			*allele1 = 0;
			*allele2 = 1;
			break;
		case 2:  // Heterozygous
			*allele1 = 1;
			*allele2 = 0;
			break;
		case 3:  // Homozygous alternate
			*allele1 = 1;
			*allele2 = 1;
			break;
		case 4:  // Missing - randomly assign
		default:
			*allele1 = random_int(0, 1);
			*allele2 = random_int(0, 1);
			break;
	}
}

int create_directory_if_needed(const char* path);  /* forward */

/**
 * @brief Generate genotype data for one chromosome (same data used for HAP and PED).
 * Individual 0 = offspring of 1 and 2. One recombination event per chr: at crossover_snp
 * the transmitted haplotype from each parent switches.
 */
void generate_chromosome_data(int chr, double* maf_array, unsigned char (*data)[NUM_SNPS_PER_CHR][2])
{
	int which_hap_from_p1 = rand() % 2;
	int which_hap_from_p2 = rand() % 2;
	/* One recombination: crossover between SNP 1 and NUM_SNPS_PER_CHR-1 */
	int crossover_snp = 1 + rand() % (NUM_SNPS_PER_CHR - 2);
	
	for (int snp = 0; snp < NUM_SNPS_PER_CHR; snp++)
	{
		int geno_parent1 = generate_genotype(maf_array[snp]);
		int geno_parent2 = generate_genotype(maf_array[snp]);
		int p1_a1, p1_a2, p2_a1, p2_a2;
		genotype_to_alleles(geno_parent1, &p1_a1, &p1_a2);
		genotype_to_alleles(geno_parent2, &p2_a1, &p2_a2);
		
		/* After crossover, switch which haplotype each parent transmits */
		int h1 = which_hap_from_p1;
		int h2 = which_hap_from_p2;
		if (snp >= crossover_snp) {
			h1 = 1 - which_hap_from_p1;
			h2 = 1 - which_hap_from_p2;
		}
		
		for (int indiv = 0; indiv < NUM_INDIVIDUALS; indiv++)
		{
			if (indiv == IDX_OFFSPRING)
			{
				data[indiv][snp][0] = (unsigned char)(h1 ? p1_a2 : p1_a1);
				data[indiv][snp][1] = (unsigned char)(h2 ? p2_a2 : p2_a1);
			}
			else if (indiv == IDX_PARENT1)
			{
				data[indiv][snp][0] = (unsigned char)p1_a1;
				data[indiv][snp][1] = (unsigned char)p1_a2;
			}
			else if (indiv == IDX_PARENT2)
			{
				data[indiv][snp][0] = (unsigned char)p2_a1;
				data[indiv][snp][1] = (unsigned char)p2_a2;
			}
			else
			{
				int g = generate_genotype(maf_array[snp]);
				int a1, a2;
				genotype_to_alleles(g, &a1, &a2);
				data[indiv][snp][0] = (unsigned char)a1;
				data[indiv][snp][1] = (unsigned char)a2;
			}
		}
	}
}

/**
 * @brief Write recombinaisonWindows.txt (window coordinates for acrossCHRphasing).
 * Format: one line per chr (0-22), "chr b0 b1 b2 b3 b4". -1 = end of chr.
 * Windows: 0-500 and 500-1000 per chromosome.
 */
static int write_window_file(const char* output_path)
{
	char filename[512];
	const char *last = strrchr(output_path, '/');
#ifdef _WIN32
	if (!last) last = strrchr(output_path, '\\');
#endif
	if (last) {
		size_t len = (size_t)(last - output_path + 1);
		if (len + 30 < sizeof(filename)) {
			strncpy(filename, output_path, len);
			filename[len] = '\0';
			strcat(filename, "recombinaisonWindows.txt");
		} else
			strcpy(filename, "recombinaisonWindows.txt");
	} else
		strcpy(filename, "recombinaisonWindows.txt");

	create_directory_if_needed(filename);
	FILE *fp = fopen(filename, "w");
	if (!fp) {
		printf("ERROR: Cannot create %s\n", filename);
		return 1;
	}
	fprintf(fp, "# Segment boundary positions. Format: chr bound1 bound2 bound3 bound4 bound5\n");
	fprintf(fp, "# -1 = end of chromosome. Windows: 0-%d and %d-%d\n",
		BOUNDS_BREAK_AT_SNP, BOUNDS_BREAK_AT_SNP, NUM_SNPS_PER_CHR);
	fprintf(fp, "0 0 0 0 0\n");
	for (int chr = 1; chr <= NUM_CHROMOSOMES; chr++)
		fprintf(fp, "%d %d -1 0 0 0\n", chr, BOUNDS_BREAK_AT_SNP);
	fclose(fp);
	printf("Wrote window file: %s (windows 0-%d, %d-%d per chr)\n", filename, BOUNDS_BREAK_AT_SNP, BOUNDS_BREAK_AT_SNP, NUM_SNPS_PER_CHR);
	return 0;
}

/**
 * @brief Create directory if it doesn't exist
 * @param path File path (directory will be extracted from this)
 * @return 0 on success, 1 on error
 */
int create_directory_if_needed(const char* path)
{
	// Extract directory path from file path
	char dir_path[300];
	strcpy(dir_path, path);
	
	// Find last '/' or '\' separator
	char* last_sep = NULL;
#ifdef _WIN32
	last_sep = strrchr(dir_path, '\\');
	if (last_sep == NULL)
		last_sep = strrchr(dir_path, '/');
#else
	last_sep = strrchr(dir_path, '/');
#endif
	
	if (last_sep != NULL)
	{
		// Terminate string at separator to get directory path
		*last_sep = '\0';
		
		// Try to create directory (mkdir returns 0 on success, -1 if exists or error)
		if (mkdir(dir_path, 0755) != 0)
		{
			// Check if error is because directory already exists
			// On most systems, EEXIST means directory exists (which is OK)
			// We'll ignore the error if directory already exists
		}
	}
	
	return 0;
}

/**
 * @brief Write a HAP file from pre-generated genotype data
 */
int write_hap_file(const char* output_path, int chr, double* maf_array, unsigned char (*data)[NUM_SNPS_PER_CHR][2])
{
	char filename[300];
	FILE* fp;
	
	strcpy(filename, output_path);
	char chr_str[10];
	sprintf(chr_str, "%d", chr);
	strcat(filename, chr_str);
	strcat(filename, ".hap");
	
	create_directory_if_needed(filename);
	fp = fopen(filename, "w");
	if (fp == NULL)
	{
		printf("ERROR: Could not open file %s for writing\n", filename);
		return 1;
	}
	
	printf("Writing HAP file: %s (chromosome %d, %d SNPs, %d individuals)\n", 
		filename, chr, NUM_SNPS_PER_CHR, NUM_INDIVIDUALS);
	
	for (int snp = 0; snp < NUM_SNPS_PER_CHR; snp++)
	{
		char rsid[50];
		generate_rsid(chr, snp, rsid);
		int position = generate_position(chr, snp);
		char ancestral_allele = 'A';
		char alternative_allele = 'G';
		
		fprintf(fp, "%d %s %d %c %c", chr, rsid, position, ancestral_allele, alternative_allele);
		
		for (int indiv = 0; indiv < NUM_INDIVIDUALS; indiv++)
		{
			fprintf(fp, " %d %d", data[indiv][snp][0], data[indiv][snp][1]);
		}
		fprintf(fp, "\n");
		
		if ((snp + 1) % 100 == 0)
			printf("  HAP progress: %d/%d SNPs\n", snp + 1, NUM_SNPS_PER_CHR);
	}
	
	fclose(fp);
	printf("Successfully generated HAP file: %s\n", filename);
	return 0;
}

/**
 * @brief Generate a MAF file for a specific chromosome
 * @param output_path Base path for output files
 * @param chr Chromosome number (1-22)
 * @param maf_array Optional: Pre-computed MAF values (NULL to generate randomly)
 * @return 0 on success, 1 on error
 */
int generate_maf_file(const char* output_path, int chr, double* maf_array)
{
	char filename[300];
	FILE* fp;
	
	// Construct filename: <output_path><chr>.maf
	strcpy(filename, output_path);
	char chr_str[10];
	sprintf(chr_str, "%d", chr);
	strcat(filename, chr_str);
	strcat(filename, ".maf");
	
	// Create directory if needed
	create_directory_if_needed(filename);
	
	// Open file for writing
	fp = fopen(filename, "w");
	if (fp == NULL)
	{
		printf("ERROR: Could not open file %s for writing\n", filename);
		return 1;
	}
	
	printf("Generating MAF file: %s (chromosome %d, %d SNPs)\n", 
		filename, chr, NUM_SNPS_PER_CHR);
	
	// Write header
	fprintf(fp, "# Chromosome %d - Minor Allele Frequency file\n", chr);
	fprintf(fp, "# SNP_ID\tPosition\tMAF\n");
	
	// Generate MAF values for each SNP
	for (int snp = 0; snp < NUM_SNPS_PER_CHR; snp++)
	{
		char rsid[50];
		generate_rsid(chr, snp, rsid);
		int position = generate_position(chr, snp);
		
		// Use pre-computed MAF if available, otherwise generate randomly
		double maf;
		if (maf_array != NULL)
		{
			maf = maf_array[snp];
		}
		else
		{
			// Generate allele frequency between 0.05 and 0.95, then compute MAF
			maf = generate_allele_frequency();
		}
		
		// Write SNP information
		fprintf(fp, "%s\t%d\t%.6f\n", rsid, position, maf);
	}
	
	fclose(fp);
	printf("Successfully generated MAF file: %s\n", filename);
	return 0;
}

/**
 * @brief Write a PED file from pre-generated genotype data
 */
int write_ped_file(const char* output_path, int chr, unsigned char (*data)[NUM_SNPS_PER_CHR][2])
{
	char filename[300];
	FILE* fp;
	
	strcpy(filename, output_path);
	char chr_str[10];
	sprintf(chr_str, "%d", chr);
	strcat(filename, chr_str);
	strcat(filename, ".ped");
	
	create_directory_if_needed(filename);
	fp = fopen(filename, "w");
	if (fp == NULL)
	{
		printf("ERROR: Could not open file %s for writing\n", filename);
		return 1;
	}
	
	printf("Writing PED file: %s (chromosome %d, %d SNPs, %d individuals)\n", 
		filename, chr, NUM_SNPS_PER_CHR, NUM_INDIVIDUALS);
	
	for (int indiv = 0; indiv < NUM_INDIVIDUALS; indiv++)
	{
		int father_id = (indiv == IDX_OFFSPRING) ? IDX_PARENT1 : 0;
		int mother_id = (indiv == IDX_OFFSPRING) ? IDX_PARENT2 : 0;
		fprintf(fp, "%d %d %d %d -9 0", indiv, indiv, father_id, mother_id);
		for (int snp = 0; snp < NUM_SNPS_PER_CHR; snp++)
			fprintf(fp, " %d %d", data[indiv][snp][0], data[indiv][snp][1]);
		fprintf(fp, "\n");
		
		if ((indiv + 1) % 100 == 0)
			printf("  PED progress: %d/%d individuals\n", indiv + 1, NUM_INDIVIDUALS);
	}
	
	fclose(fp);
	printf("Successfully wrote PED file: %s\n", filename);
	return 0;
}

/**
 * @brief Main function
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return 0 on success, 1 on error
 */
int main(int argc, char* argv[])
{
	// Check command line arguments
	if (argc < 2)
	{
		printf("Usage: %s <output_path> [format]\n", argv[0]);
		printf("  output_path: Base path for output files (e.g., 'data/chr')\n");
		printf("  format: Optional - 'HAP', 'PED', or 'BOTH' (default) - HAP and PED = same data\n");
		printf("\n");
		printf("This program generates random data for:\n");
		printf("  - %d individuals\n", NUM_INDIVIDUALS);
		printf("  - %d SNPs per chromosome\n", NUM_SNPS_PER_CHR);
		printf("  - Chromosomes 1-%d\n", NUM_CHROMOSOMES);
		return 1;
	}
	
	const char* output_path = argv[1];
	int output_format = 2;  // 0 = HAP only, 1 = PED only, 2 = BOTH (default)
	
	// Check for format argument
	if (argc >= 3)
	{
		if (strcmp(argv[2], "PED") == 0 || strcmp(argv[2], "ped") == 0)
			output_format = 1;
		else if (strcmp(argv[2], "HAP") == 0 || strcmp(argv[2], "hap") == 0)
			output_format = 0;
		else if (strcmp(argv[2], "BOTH") == 0 || strcmp(argv[2], "both") == 0)
			output_format = 2;
		else
			printf("WARNING: Unknown format '%s', using default (BOTH)\n", argv[2]);
	}
	
	// Initialize random number generator
	srand((unsigned int)time(NULL));
	
	printf("========================================\n");
	printf("Random Data Generator for Phasing Program\n");
	printf("========================================\n");
	printf("Output path: %s\n", output_path);
	printf("Output format: %s\n", output_format == 0 ? "HAP" : (output_format == 1 ? "PED" : "BOTH (HAP + PED, same data)"));
	printf("Individuals: %d\n", NUM_INDIVIDUALS);
	printf("SNPs per chromosome: %d\n", NUM_SNPS_PER_CHR);
	printf("Chromosomes: 1-%d\n", NUM_CHROMOSOMES);
	printf("Trio: individual 0 = offspring of individuals 1 and 2\n");
	printf("  Detail: Parent1=indiv1, Parent2=indiv2. Per chromosome: one recombination event;\n");
	printf("  at a random crossover SNP, the transmitted haplotype from each parent switches.\n");
	printf("  Before crossover: offspring gets hap0 or hap1 from P1 and P2.\n");
	printf("  After crossover: offspring gets the other haplotype from each parent.\n");
	printf("========================================\n\n");

	// Write window coordinates file (recombinaisonWindows.txt: windows 0-500, 500-1000)
	if (write_window_file(output_path) != 0) {
		printf("WARNING: Failed to write window file\n");
	}
	
	// Generate files for each chromosome
	int error_count = 0;
	unsigned char chrom_data[NUM_INDIVIDUALS][NUM_SNPS_PER_CHR][2];
	
	for (int chr = 1; chr <= NUM_CHROMOSOMES; chr++)
	{
		printf("\n--- Processing Chromosome %d ---\n", chr);
		
		// Pre-compute MAF values for this chromosome
		double maf_array[NUM_SNPS_PER_CHR];
		for (int snp = 0; snp < NUM_SNPS_PER_CHR; snp++)
			maf_array[snp] = generate_allele_frequency();
		
		// Generate genotype data once (same data for HAP and PED)
		printf("  Generating genotype data...\n");
		generate_chromosome_data(chr, maf_array, chrom_data);
		
		// Write HAP file if requested
		if (output_format == 0 || output_format == 2)
		{
			if (write_hap_file(output_path, chr, maf_array, chrom_data) != 0)
			{
				error_count++;
				printf("ERROR: Failed to write HAP file for chromosome %d\n", chr);
			}
		}
		
		// Write PED file if requested
		if (output_format == 1 || output_format == 2)
		{
			if (write_ped_file(output_path, chr, chrom_data) != 0)
			{
				error_count++;
				printf("ERROR: Failed to write PED file for chromosome %d\n", chr);
			}
		}
		
		// Generate MAF file (always)
		if (generate_maf_file(output_path, chr, maf_array) != 0)
		{
			error_count++;
			printf("ERROR: Failed to generate MAF file for chromosome %d\n", chr);
		}
		
		printf("Chromosome %d: Complete\n", chr);
	}
	
	printf("\n========================================\n");
	if (error_count == 0)
	{
		printf("SUCCESS: All files generated successfully!\n");
		if (output_format == 2)
			printf("Generated %d HAP files, %d PED files, %d MAF files (HAP and PED = same data)\n",
				NUM_CHROMOSOMES, NUM_CHROMOSOMES, NUM_CHROMOSOMES);
		else if (output_format == 0)
			printf("Generated %d HAP files and %d MAF files\n", NUM_CHROMOSOMES, NUM_CHROMOSOMES);
		else
			printf("Generated %d PED files and %d MAF files\n", NUM_CHROMOSOMES, NUM_CHROMOSOMES);
	}
	else
	{
		printf("WARNING: %d error(s) occurred during generation\n", error_count);
	}
	printf("========================================\n");
	
	if (error_count == 0) {
		printf("RETURN(0): GenerateRandomData completed successfully.\n");
		return 0;
	} else {
		printf("RETURN(1): GenerateRandomData failed with %d errors.\n", error_count);
		return 1;
	}
}

