/**
 * @file GenerateRandomData.cpp
 * @brief Random data generator for haplotype phasing program
 * 
 * This program generates random HAP, PED, and MAF files for testing the phasing algorithm.
 * It generates data for 1000 individuals and 1000 SNPs per chromosome (chromosomes 1-22).
 * 
 * Usage:
 *   ./GenerateRandomData <output_path> [format]
 * 
 * Parameters:
 *   output_path: Base path for output files (e.g., "data/chr")
 *   format: Optional - "HAP" (default) or "PED" to specify output format
 * 
 * Output files:
 *   HAP format (default):
 *     - <output_path>1.hap, <output_path>2.hap, ..., <output_path>22.hap (HAP files)
 *     - <output_path>1.maf, <output_path>2.maf, ..., <output_path>22.maf (MAF files)
 *   PED format:
 *     - <output_path>1.ped, <output_path>2.ped, ..., <output_path>22.ped (PED files)
 *     - <output_path>1.maf, <output_path>2.maf, ..., <output_path>22.maf (MAF files)
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
 * @brief Generate a HAP file for a specific chromosome
 * @param output_path Base path for output files
 * @param chr Chromosome number (1-22)
 * @param maf_array Optional: Pre-computed MAF values (NULL to generate randomly)
 * @return 0 on success, 1 on error
 */
int generate_hap_file(const char* output_path, int chr, double* maf_array)
{
	char filename[300];
	FILE* fp;
	
	// Construct filename: <output_path><chr>.hap
	strcpy(filename, output_path);
	char chr_str[10];
	sprintf(chr_str, "%d", chr);
	strcat(filename, chr_str);
	strcat(filename, ".hap");
	
	// Create directory if needed
	create_directory_if_needed(filename);
	
	// Open file for writing
	fp = fopen(filename, "w");
	if (fp == NULL)
	{
		printf("ERROR: Could not open file %s for writing\n", filename);
		return 1;
	}
	
	printf("Generating HAP file: %s (chromosome %d, %d SNPs, %d individuals)\n", 
		filename, chr, NUM_SNPS_PER_CHR, NUM_INDIVIDUALS);
	
	// Generate or use pre-computed MAF for this chromosome
	double maf[NUM_SNPS_PER_CHR];
	if (maf_array != NULL)
	{
		// Use pre-computed MAF values
		for (int snp = 0; snp < NUM_SNPS_PER_CHR; snp++)
		{
			maf[snp] = maf_array[snp];
		}
	}
	else
	{
		// Generate random allele frequencies and compute MAF
		// Allele frequencies are between MIN_ALLELE_FREQ (0.05) and MAX_ALLELE_FREQ (0.95)
		// MAF will be between 0.05 and 0.5
		for (int snp = 0; snp < NUM_SNPS_PER_CHR; snp++)
		{
			maf[snp] = generate_allele_frequency();
		}
	}
	
	// Generate each SNP in Relate .haps format
	// Format: Chromosome SNP_ID Position Ancestral_allele Alternative_allele allele1_indiv1 allele2_indiv1 ...
	// Reference: https://myersgroup.github.io/relate/input_data.html
	for (int snp = 0; snp < NUM_SNPS_PER_CHR; snp++)
	{
		// Generate rsID
		char rsid[50];
		generate_rsid(chr, snp, rsid);
		
		// Generate position
		int position = generate_position(chr, snp);
		
		// Generate ancestral and alternative alleles (randomly choose A, C, G, T)
		// For simplicity, we'll use A as ancestral and G as alternative
		// In real data, these would come from reference genome
		char ancestral_allele = 'A';
		char alternative_allele = 'G';
		
		// Write Relate .haps format:
		// Column 1: Chromosome number [integer]
		// Column 2: SNP ID [string]
		// Column 3: SNP position [integer]
		// Column 4: Ancestral allele [char] (represented as 0 in haplotypes)
		// Column 5: Alternative allele [char] (represented as 1 in haplotypes)
		// Columns 6+: Two alleles per individual (0 = ancestral, 1 = alternative)
		fprintf(fp, "%d %s %d %c %c", chr, rsid, position, ancestral_allele, alternative_allele);
		
		// Debug: count genotypes for first SNP
		int debug_geno_count[5] = {0, 0, 0, 0, 0};
		int debug_alt_alleles = 0;
		
		// Generate and write haplotypes for all individuals
		for (int indiv = 0; indiv < NUM_INDIVIDUALS; indiv++)
		{
			// Generate genotype based on MAF
			int genotype = generate_genotype(maf[snp]);
			
			// Debug: count genotypes
			if (snp == 0 && chr == 1 && genotype < 5)
			{
				debug_geno_count[genotype]++;
			}
			
			// Convert to allele pair
			int allele1, allele2;
			genotype_to_alleles(genotype, &allele1, &allele2);
			
			// Debug: count alternate alleles
			if (snp == 0 && chr == 1)
			{
				debug_alt_alleles += allele1 + allele2;
			}
			
			// In Relate format: 0 = ancestral allele, 1 = alternative allele
			// Write both alleles separated by space
			fprintf(fp, " %d %d", allele1, allele2);
		}
		
		// Debug: print statistics for first SNP
		if (snp == 0 && chr == 1)
		{
			printf("  DEBUG SNP 0: MAF=%.6f, genotypes: 0=%d, 1=%d, 2=%d, 3=%d, 4=%d, alt_alleles=%d/2000\n",
				maf[snp], debug_geno_count[0], debug_geno_count[1], debug_geno_count[2], 
				debug_geno_count[3], debug_geno_count[4], debug_alt_alleles);
		}
		
		fprintf(fp, "\n");
		
		// Progress indicator
		if ((snp + 1) % 100 == 0)
		{
			printf("  Progress: %d/%d SNPs\n", snp + 1, NUM_SNPS_PER_CHR);
		}
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
 * @brief Generate a PED file for a specific chromosome
 * @param output_path Base path for output files
 * @param chr Chromosome number (1-22)
 * @param maf_array Pre-computed MAF values for each SNP
 * @return 0 on success, 1 on error
 */
int generate_ped_file(const char* output_path, int chr, double maf_array[])
{
	char filename[300];
	FILE* fp;
	
	// Construct filename: <output_path><chr>.ped
	strcpy(filename, output_path);
	char chr_str[10];
	sprintf(chr_str, "%d", chr);
	strcat(filename, chr_str);
	strcat(filename, ".ped");
	
	// Create directory if needed
	create_directory_if_needed(filename);
	
	// Open file for writing
	fp = fopen(filename, "w");
	if (fp == NULL)
	{
		printf("ERROR: Could not open file %s for writing\n", filename);
		return 1;
	}
	
	printf("Generating PED file: %s (chromosome %d, %d SNPs, %d individuals)\n", 
		filename, chr, NUM_SNPS_PER_CHR, NUM_INDIVIDUALS);
	
	// Generate each individual
	for (int indiv = 0; indiv < NUM_INDIVIDUALS; indiv++)
	{
		// Write PED header: Family_ID Individual_ID Father_ID Mother_ID Sex Phenotype
		fprintf(fp, "%d %d 0 0 -9 0", indiv, indiv);
		
		// Generate and write genotypes for all SNPs
		for (int snp = 0; snp < NUM_SNPS_PER_CHR; snp++)
		{
			// Generate genotype based on MAF
			int genotype = generate_genotype(maf_array[snp]);
			
			// Convert to allele pair
			int allele1, allele2;
			genotype_to_alleles(genotype, &allele1, &allele2);
			
			// Write alleles in PED format: space allele1 space allele2
			fprintf(fp, " %d %d", allele1, allele2);
		}
		
		fprintf(fp, "\n");
		
		// Progress indicator
		if ((indiv + 1) % 100 == 0)
		{
			printf("  Progress: %d/%d individuals\n", indiv + 1, NUM_INDIVIDUALS);
		}
	}
	
	fclose(fp);
	printf("Successfully generated PED file: %s\n", filename);
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
		printf("  format: Optional - 'HAP' (default) or 'PED' to specify output format\n");
		printf("\n");
		printf("This program generates random data for:\n");
		printf("  - %d individuals\n", NUM_INDIVIDUALS);
		printf("  - %d SNPs per chromosome\n", NUM_SNPS_PER_CHR);
		printf("  - Chromosomes 1-%d\n", NUM_CHROMOSOMES);
		return 1;
	}
	
	const char* output_path = argv[1];
	int output_format = 0;  // 0 = HAP, 1 = PED
	
	// Check for format argument
	if (argc >= 3)
	{
		if (strcmp(argv[2], "PED") == 0 || strcmp(argv[2], "ped") == 0)
			output_format = 1;
		else if (strcmp(argv[2], "HAP") == 0 || strcmp(argv[2], "hap") == 0)
			output_format = 0;
		else
		{
			printf("WARNING: Unknown format '%s', using default (HAP)\n", argv[2]);
		}
	}
	
	// Initialize random number generator
	srand((unsigned int)time(NULL));
	
	printf("========================================\n");
	printf("Random Data Generator for Phasing Program\n");
	printf("========================================\n");
	printf("Output path: %s\n", output_path);
	printf("Output format: %s\n", output_format == 0 ? "HAP" : "PED");
	printf("Individuals: %d\n", NUM_INDIVIDUALS);
	printf("SNPs per chromosome: %d\n", NUM_SNPS_PER_CHR);
	printf("Chromosomes: 1-%d\n", NUM_CHROMOSOMES);
	printf("========================================\n\n");
	
	// Generate files for each chromosome
	int error_count = 0;
	for (int chr = 1; chr <= NUM_CHROMOSOMES; chr++)
	{
		printf("\n--- Processing Chromosome %d ---\n", chr);
		
		// Pre-compute MAF values for this chromosome (used for both HAP and PED)
		// Allele frequencies are generated between 0.05 and 0.95, then converted to MAF (0.05 to 0.5)
		double maf_array[NUM_SNPS_PER_CHR];
		for (int snp = 0; snp < NUM_SNPS_PER_CHR; snp++)
		{
			maf_array[snp] = generate_allele_frequency();
		}
		
		// Generate genotype file based on format
		if (output_format == 0)  // HAP format
		{
			// Generate HAP file using pre-computed MAF
			if (generate_hap_file(output_path, chr, maf_array) != 0)
			{
				error_count++;
				printf("ERROR: Failed to generate HAP file for chromosome %d\n", chr);
				continue;
			}
		}
		else  // PED format
		{
			// Generate PED file using pre-computed MAF
			if (generate_ped_file(output_path, chr, maf_array) != 0)
			{
				error_count++;
				printf("ERROR: Failed to generate PED file for chromosome %d\n", chr);
				continue;
			}
		}
		
		// Generate MAF file (always generated, using pre-computed MAF values)
		if (generate_maf_file(output_path, chr, maf_array) != 0)
		{
			error_count++;
			printf("ERROR: Failed to generate MAF file for chromosome %d\n", chr);
			continue;
		}
		
		printf("Chromosome %d: Complete\n", chr);
	}
	
	printf("\n========================================\n");
	if (error_count == 0)
	{
		printf("SUCCESS: All files generated successfully!\n");
		if (output_format == 0)
		{
			printf("Generated %d HAP files and %d MAF files\n", 
				NUM_CHROMOSOMES, NUM_CHROMOSOMES);
		}
		else
		{
			printf("Generated %d PED files and %d MAF files\n", 
				NUM_CHROMOSOMES, NUM_CHROMOSOMES);
		}
	}
	else
	{
		printf("WARNING: %d error(s) occurred during generation\n", error_count);
	}
	printf("========================================\n");
	
	return (error_count == 0) ? 0 : 1;
}

