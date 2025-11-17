/**
 * @file phasing_program.cpp
 * @brief Implementation of the PhasingProgram class
 * 
 * This file contains the implementation of all methods for the PhasingProgram class,
 * which is the main class for performing haplotype phasing across chromosomes.
 * 
 * The class handles:
 * - Reading input files in HAP or PED format
 * - Loading parent information and individual lists
 * - Performing phasing algorithms (with or without parent information)
 * - Writing output files
 * - Orchestrating the entire phasing workflow
 */

#include "phasing_program.h"
#include "libfileio/fileio.h"
#include "libphasing/phasing.h"

/**
 * @brief Constructor for PhasingProgram class
 * 
 * Initializes all member variables to their default values:
 * - Relationship tracking variables set to 0 or -1
 * - File format flags: HAP input (0), PED output (1)
 * - All arrays initialized to zero
 * - Genome pointers set to nullptr (memory allocated later)
 * 
 * The constructor uses an initializer list for efficient initialization
 * of primitive types, then initializes arrays in the constructor body.
 */
PhasingProgram::PhasingProgram() : IDjob(0), placefirttoconsider(0), printdetail(0), NbIndiv(0),
	nbparentinfo(0), nbIndivToProcess(0), useListIndiv(0), nbbreak(0), nbrelatpihat(0),
	inputFormat(0), outputFormat(1)  // Default: HAP input, PED output
{
	// Initialize MAF (Minor Allele Frequency) array to zero
	// MAF values will be calculated when reading input files
	for (int i = 0; i < NSNPPERCHR; i++)
	{
		for (int j = 0; j < 23; j++)
		{
			MAF[i][j] = 0;
		}
	}
	
	// Initialize genome pointers to nullptr
	// Memory will be allocated when reading input files based on NbIndiv
	// and nbsnpperchr values
	for (int i = 0; i < 23; i++)
	{
		genomes[i] = nullptr;
	}
}

/**
 * @brief Read genome data from HAP format file
 * 
 * This is a wrapper method that calls the library function readgenomelocal()
 * from libfileio, passing the necessary member variables. The HAP format
 * is SNP-centric: each line represents one SNP with genotypes for all individuals.
 * 
 * @param pathfile Base path for input files (e.g., "data/chr")
 * @param chr Chromosome number (1-22)
 * @param run Run identifier (typically 100+chr, used for file naming)
 * @param step Step identifier (typically 105, used for file naming)
 * @param gentomodify Pointer to memory buffer where genome data will be stored
 * @return 0 on success, 1 on error
 */
int PhasingProgram::readgenomelocal(char pathfile[], int chr, int run, int step, unsigned char* gentomodify)
{
	// Delegate to library function, passing member variables
	return ::readgenomelocal(pathfile, chr, run, step, gentomodify, NbIndiv, nbsnpperchr, nbsnpperchrinfile);
}

/**
 * @brief Read parent information from file
 * 
 * Reads a file containing parent-child relationships. The file format is:
 * Individual_ID Parent1_ID Parent2_ID
 * 
 * Use -1 to indicate unknown parents. This information is used to determine
 * which phasing algorithm to use (with or without parents).
 * 
 * @param filename Path to the parent information file
 * @return 0 on success, 1 on error
 */
int PhasingProgram::readParentInfo(const char* filename)
{
	// Delegate to library function, storing results in member variables
	return ::readParentInfo(filename, parentinfo, &nbparentinfo, NBINDIVMAX);
}

/**
 * @brief Get parent information for a specific individual
 * 
 * Looks up the parent IDs for a given individual from the loaded parent
 * information. If the individual is not found or has no known parents,
 * parent1 and parent2 are set to -1.
 * 
 * @param indivID Individual ID to look up (must be in range [0, NbIndiv-1])
 * @param parent1 Output parameter: Parent 1 ID (or -1 if unknown)
 * @param parent2 Output parameter: Parent 2 ID (or -1 if unknown)
 * @return 0 if parents found, 1 if individual not found or no parents
 */
int PhasingProgram::getParentInfo(int indivID, int* parent1, int* parent2)
{
	// Delegate to library function, using loaded parent information
	return ::getParentInfo(indivID, parent1, parent2, parentinfo, nbparentinfo);
}

/**
 * @brief Read list of individuals to process
 * 
 * Reads a file containing a list of individual IDs that should be processed.
 * This allows selective processing of specific individuals rather than
 * processing all individuals in the dataset.
 * 
 * File format: one individual ID per line. Lines starting with '#' are
 * treated as comments. Empty lines are ignored.
 * 
 * @param filename Path to the individual list file
 * @return 0 on success, 1 on error
 */
int PhasingProgram::readListIndiv(const char* filename)
{
	// Delegate to library function, storing results in member variables
	return ::readListIndiv(filename, listIndivToProcess, &nbIndivToProcess, &useListIndiv, NbIndiv, NBINDIVMAX);
}

int PhasingProgram::findrelative(int ID,int numtrio,int IDp1loop,int IDp2loop,int lenminseg,int version,int gentostart,char pathresult[])
{
	return ::findrelative(ID, numtrio, IDp1loop, IDp2loop, lenminseg, version, gentostart, pathresult,
		NbIndiv, nbsnpperchr, nbsnpperchrinfile, MAF, genomes, 
		pihatagainstall, bestpihatagainstallID, &IDbestpihat, &IDbestpihat2, 
		bestpihat, &placefirttoconsider, seuilpihat);
}

/**
 * @brief Perform phasing without parent information
 * 
 * This method phases haplotypes for an individual when parent information
 * is not available. It uses cross-chromosome correlations and population-level
 * information to determine phase.
 * 
 * The algorithm:
 * 1. Uses shared segments with relatives identified in findrelative()
 * 2. Analyzes correlations between chromosome segments
 * 3. Uses iterative refinement to improve phasing accuracy
 * 4. Stores phased haplotypes in the genomes array
 * 
 * This is a wrapper that delegates to the library function in libphasing.
 * 
 * @param ID Individual ID to phase
 * @param chr1 First chromosome to process (typically 1)
 * @param chr2 Last chromosome to process (0 means process all chromosomes 1-22)
 * @param numtrio Trio identifier (typically same as ID)
 * @param limitnbsnp Maximum number of SNPs to consider (0 means no limit)
 * @param run Run identifier (for tracking/versioning)
 * @param IDp1loop Parent 1 ID for loop (typically same as ID when no parents)
 * @param IDp2loop Parent 2 ID for loop (typically same as ID when no parents)
 * @param lenminseg Minimum segment length to consider
 * @param version Algorithm version to use
 * @param gentostart Starting genotype index
 * @param pathresult Base path for result files
 * @return 0 on success, 1 on error
 */
int PhasingProgram::predict_phasing_without_parents(int ID,int chr1, int chr2,int numtrio,float limitnbsnp,int run,int IDp1loop,int IDp2loop,int lenminseg,int version,int gentostart,char pathresult[])
{
	// Initialize local variables that are not class members
	float pihatagainstall2[MAXPOP] = {0};  // Second PI-HAT array (initialized to zero)
	int64_t relatpihatchr[MAXCLOSERELAT][23][2];  // PI-HAT per chromosome and haplotype
	memset(relatpihatchr, 0, sizeof(relatpihatchr));  // Initialize to zero
	
	// Delegate to library function, passing all necessary member variables
	return ::predict_phasing_without_parents(ID, chr1, chr2, numtrio, limitnbsnp, run, IDp1loop, IDp2loop, 
		lenminseg, version, gentostart, pathresult, NbIndiv, nbsnpperchr, nbsnpperchrinfile, MAF, genomes, 
		genomeoffpss, seuilpihat, chrdivider, nbchrdivider, &IDjob, nbbreak, nbrelatpihat,
		pihatagainstall, bestpihatagainstallID, pihatagainstall2, tappointdec, relatpihatchr, &placefirttoconsider);
}

int PhasingProgram::predict_phasing_with_parents_providing_GT(int ID,int chr1, int chr2,int numtrio,float limitnbsnp,int run,int IDp1loop,int IDp2loop,int lenminseg,int version,int gentostart,char pathresult[])
{
	// Initialize local variables that are not class members
	float pihatagainstall2[MAXPOP] = {0};  // Second PI-HAT array (initialized to zero)
	int64_t relatpihatchr[MAXCLOSERELAT][23][2];  // PI-HAT per chromosome and haplotype
	memset(relatpihatchr, 0, sizeof(relatpihatchr));  // Initialize to zero
	
	return ::predict_phasing_with_parents_providing_GT(ID, chr1, chr2, numtrio, limitnbsnp, run, IDp1loop, IDp2loop, 
		lenminseg, version, gentostart, pathresult, NbIndiv, nbsnpperchr, nbsnpperchrinfile, MAF, genomes, 
		genomeoffpss, seuilpihat, chrdivider, nbchrdivider, &IDjob, nbbreak, nbrelatpihat,
		pihatagainstall, bestpihatagainstallID, pihatagainstall2, tappointdec, relatpihatchr, &placefirttoconsider);
}

/**
 * @brief Read a PED file for a specific chromosome
 * 
 * Reads genotype data from a PED format file for a specific chromosome.
 * PED format is individual-centric: each line represents one individual
 * with genotypes for all SNPs.
 * 
 * This is a wrapper method that calls the library function readpedfile()
 * from libfileio, passing the necessary member variables.
 * 
 * @param filename Path to the PED file (should include chromosome number, e.g., "data/chr1.ped")
 * @param chr Chromosome number (1-22)
 * @return 0 on success, 1 on error
 */
int PhasingProgram::readpedfile_local(const char* filename, int chr)
{
	// Delegate to library function, passing member variables
	// The function will update nbsnpperchrinfile[chr] and NbIndiv
	return ::readpedfile(filename, chr, genomes, nbsnpperchr, &nbsnpperchrinfile[chr], &NbIndiv, NBINDIVMAX);
}

/**
 * @brief Write phased output to file
 * 
 * Writes the phased haplotypes for a given chromosome to an output file.
 * The output format is determined by the outputFormat member variable:
 * - outputFormat == 1: PED format (default)
 * - outputFormat == 0: HAP format (currently not fully implemented, writes as PED)
 * 
 * The phased data is read from the genomes array, which contains the
 * phased haplotypes after running the phasing algorithm.
 * 
 * @param chr Chromosome number (1-22)
 * @param pathoutput Base path for output files (e.g., "results/out")
 *                   Output file will be: <pathoutput><chr>.ped
 * @return 0 on success, 1 on error
 */
int PhasingProgram::writeoutput(int chr, char pathoutput[])
{
	if (outputFormat == 1)  // PED format (default)
	{
		// Write in PED format: one individual per line
		return ::writepedfile(chr, pathoutput, genomes, nbsnpperchr, nbsnpperchrinfile, NbIndiv);
	}
	else  // HAP format (not yet fully implemented, fallback to PED)
	{
		// TODO: Implement HAP format output
		// For now, write as PED format
		return ::writepedfile(chr, pathoutput, genomes, nbsnpperchr, nbsnpperchrinfile, NbIndiv);
	}
}

/**
 * @brief Main execution method - orchestrates the entire phasing process
 * 
 * This is the main method that orchestrates the entire phasing workflow:
 * 
 * 1. **Command-line argument parsing**: Parses all command-line options
 *    - Required: -NbIndiv, -PathInput, -PathOutput
 *    - Optional: -PathParentInfo, -ListIndiv, -InputFormat, -OutputFormat
 * 
 * 2. **Initialization**: 
 *    - Sets default values for chromosome SNP counts
 *    - Initializes chromosome dividers
 *    - Sets up random number generator
 *    - Initializes result arrays
 * 
 * 3. **Input file reading**:
 *    - Reads parent information file (if provided)
 *    - Reads individual list file (if provided)
 *    - Allocates memory for genome data
 *    - Reads genome data from HAP or PED files
 * 
 * 4. **Individual processing loop**:
 *    For each individual to process:
 *    - Loads chromosome segments and identifies relatives
 *    - Determines phasing method (with/without parents)
 *    - Executes appropriate phasing algorithm
 *    - Writes output files for all chromosomes
 * 
 * 5. **Output generation**:
 *    - Writes phased haplotypes to output files
 *    - Format determined by outputFormat (PED by default)
 * 
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return 0 on success, non-zero on error
 */
int PhasingProgram::run(int argc, char* argv[])
	{
		// Record start time for performance measurement
		clock_t begin0 = clock();
		
		// Enable detailed output/printing
		printdetail = 1;
		
		// Print program header
		printf("\n");
		printf("================================================================================\n");
		printf("                    ACROSS CHROMOSOMES PHASING PROGRAM\n");
		printf("================================================================================\n");
		printf("\n");
		
		// Default chromosome to process (1 = chromosome 1, typically process all)
		int chrtodo1 = 1;
		
		// Initialize path strings
		char PathInput[200];      // Base path for input files
		char PathOutput[200];     // Base path for output files
		char PathParentInfo[200] = "";  // Path to parent info file (empty if not provided)
		char PathListIndiv[200] = "";   // Path to individual list file (empty if not provided)
		
		// Set default file formats
		inputFormat = 0;  // Default: HAP format (SNP-centric)
		outputFormat = 1; // Default: PED format (individual-centric)
		
		// ========================================================================
		// STEP 1: Parse command-line arguments
		// ========================================================================
		// Iterate through all command-line arguments
		for (int input = 1; input < argc; input++)
		{
			// Parse -NbIndiv: Number of individuals in the dataset
			if (strncmp(argv[input], "-NbIndiv", strlen("-NbIndiv")) == 0 && input < argc - 1)
				NbIndiv = atoi(argv[++input]);
			
			// Parse -PathInput: Base path for input files
			else if (strncmp(argv[input], "-PathInput", strlen("-PathInput")) == 0 && input < argc - 1)
				strcpy(PathInput, argv[++input]);
			
			// Parse -PathOutput: Base path for output files
			else if (strncmp(argv[input], "-PathOutput", strlen("-PathOutput")) == 0 && input < argc - 1)
				strcpy(PathOutput, argv[++input]);
			
			// Parse -PathParentInfo: Path to parent information file
			else if (strncmp(argv[input], "-PathParentInfo", strlen("-PathParentInfo")) == 0 && input < argc - 1)
				strcpy(PathParentInfo, argv[++input]);
			
			// Parse -ListIndiv: Path to individual list file
			else if (strncmp(argv[input], "-ListIndiv", strlen("-ListIndiv")) == 0 && input < argc - 1)
				strcpy(PathListIndiv, argv[++input]);
			
			// Parse -InputFormat: Input file format (HAP or PED)
			else if (strncmp(argv[input], "-InputFormat", strlen("-InputFormat")) == 0 && input < argc - 1)
			{
				input++;
				if (strcmp(argv[input], "HAP") == 0 || strcmp(argv[input], "hap") == 0)
					inputFormat = 0;  // HAP format
				else if (strcmp(argv[input], "PED") == 0 || strcmp(argv[input], "ped") == 0)
					inputFormat = 1;  // PED format
				else
					printf("WARNING: Unknown input format '%s', using default (HAP)\n", argv[input]);
			}
			
			// Parse -OutputFormat: Output file format (HAP or PED)
			else if (strncmp(argv[input], "-OutputFormat", strlen("-OutputFormat")) == 0 && input < argc - 1)
			{
				input++;
				if (strcmp(argv[input], "HAP") == 0 || strcmp(argv[input], "hap") == 0)
					outputFormat = 0;  // HAP format (not fully implemented)
				else if (strcmp(argv[input], "PED") == 0 || strcmp(argv[input], "ped") == 0)
					outputFormat = 1;  // PED format (default)
				else
					printf("WARNING: Unknown output format '%s', using default (PED)\n", argv[input]);
			}
		}
		
		// Print configuration summary
		printf("Configuration:\n");
		printf("  - Number of individuals: %d\n", NbIndiv);
		printf("  - Input format: %s\n", inputFormat == 0 ? "HAP" : "PED");
		printf("  - Output format: %s\n", outputFormat == 0 ? "HAP" : "PED");
		printf("  - Input path: %s\n", PathInput);
		printf("  - Output path: %s\n", PathOutput);
		printf("\n");
		
		// ========================================================================
		// STEP 2: Validate parameters
		// ========================================================================
		// Check that number of individuals is valid
		if (NbIndiv == 0)
		{
			printf("\n");
			printf("================================================================================\n");
			printf("ERROR: Number of individuals is zero or undefined\n");
			printf("Please specify the number of individuals using -NbIndiv <number>\n");
			printf("================================================================================\n");
			printf("\n");
			exit(1);
		}
		if (NbIndiv > NBINDIVMAX)
		{
			printf("\n");
			printf("================================================================================\n");
			printf("ERROR: Number of individuals (%d) exceeds maximum allowed (%d)\n", NbIndiv, NBINDIVMAX);
			printf("Please reduce the number of individuals or increase NBINDIVMAX in types.h\n");
			printf("================================================================================\n");
			printf("\n");
			exit(1);
		}
		
		// ========================================================================
		// STEP 3: Initialize chromosome SNP counts
		// ========================================================================
		// Set maximum number of SNPs per chromosome (based on typical human genome)
		// Index 0 is metadata, indices 1-22 are chromosomes 1-22
		nbsnpperchr[0] = 330005;  // Metadata/total
		nbsnpperchr[1] = 26229;   // Chromosome 1
		nbsnpperchr[2] = 26210;   // Chromosome 2
		nbsnpperchr[3] = 22209;   // Chromosome 3
		nbsnpperchr[4] = 20690;   // Chromosome 4
		nbsnpperchr[5] = 19027;   // Chromosome 5
		nbsnpperchr[6] = 18418;   // Chromosome 6
		nbsnpperchr[7] = 18367;   // Chromosome 7
		nbsnpperchr[8] = 16283;   // Chromosome 8
		nbsnpperchr[9] = 14990;   // Chromosome 9
		nbsnpperchr[10] = 16494;  // Chromosome 10
		nbsnpperchr[11] = 15818;  // Chromosome 11
		nbsnpperchr[12] = 16008;  // Chromosome 12
		nbsnpperchr[13] = 11510;  // Chromosome 13
		nbsnpperchr[14] = 10804;  // Chromosome 14
		nbsnpperchr[15] = 10884;  // Chromosome 15
		nbsnpperchr[16] = 12195;  // Chromosome 16
		nbsnpperchr[17] = 11486;  // Chromosome 17
		nbsnpperchr[18] = 10222;  // Chromosome 18
		nbsnpperchr[19] = 9806;   // Chromosome 19
		nbsnpperchr[20] = 8985;   // Chromosome 20
		nbsnpperchr[21] = 5227;   // Chromosome 21
		nbsnpperchr[22] = 5882;   // Chromosome 22
		
		// Initialize chromosome dividers: each chromosome starts as one segment
		// Method index 25 is used for initial/default division
		for (int chrtemp1 = 1; chrtemp1 < 23; chrtemp1++)
		{
			chrdivider[25][chrtemp1][0].start = 0;  // Segment starts at SNP 0
			chrdivider[25][chrtemp1][0].end = nbsnpperchr[chrtemp1];  // Segment ends at last SNP
			nbchrdivider[25][chrtemp1] = 1;  // One segment per chromosome initially
		}
		
		// Initialize random number generator with fixed seed for reproducibility
		srand(0);
		
		// Initialize result matrix to zero (for tracking algorithm performance)
		for (int method = 0; method < 4; method++)
		{
			for (int res = 0; res < 6; res++)
			{
				result[method][res] = 0;
			}
		}
		// Memory allocation completed
		if (printdetail) 
		{
			printf("Initialization completed.\n");
			printf("Memory allocated: %zu bytes for genome data\n", sizeof(genomes));
		}
		for (int IDrun = 0; IDrun < 12; IDrun++)
		{
			distrigametic[IDrun] = 0;
		}
		for (int IDrun = 0; IDrun < 23; IDrun++)
		{
			distrigametickeep[IDrun] = 0;
		}
		// ========================================================================
		// STEP 4: Read optional input files
		// ========================================================================
		// Read parent information file if provided
		if (strlen(PathParentInfo) > 0)
		{
			printf("Reading parent information file: %s...", PathParentInfo);
			fflush(stdout);
			if (readParentInfo(PathParentInfo) != 0)
			{
				printf(" FAILED\n");
				printf("\n");
				printf("================================================================================\n");
				printf("ERROR: Failed to read parent info file: %s\n", PathParentInfo);
				printf("Please check that the file exists and is in the correct format.\n");
				printf("Expected format: Individual_ID Parent1_ID Parent2_ID (one per line)\n");
				printf("================================================================================\n");
				printf("\n");
				exit(1);
			}
			printf(" OK\n");
		};
		
		// Read individual list file if provided (for selective processing)
		if (strlen(PathListIndiv) > 0)
		{
			printf("Reading individual list file: %s...", PathListIndiv);
			fflush(stdout);
			if (readListIndiv(PathListIndiv) != 0)
			{
				printf(" FAILED\n");
				printf("\n");
				printf("================================================================================\n");
				printf("ERROR: Failed to read individual list file: %s\n", PathListIndiv);
				printf("Please check that the file exists and is in the correct format.\n");
				printf("Expected format: One individual ID per line\n");
				printf("================================================================================\n");
				printf("\n");
				exit(1);
			}
			printf(" OK\n");
			printf("  Will process %d individuals from the list\n", nbIndivToProcess);
		}
		else
		{
			printf("As there is not an individual list file, Processing first 10 individuals in the dataset\n");
			nbIndivToProcess = 10;
			for (int idx = 0; idx < nbIndivToProcess; idx++)
			{
				listIndivToProcess[idx] = idx;
			}
			useListIndiv = 1;
		}
		printf("\n");
		
		// ========================================================================
		// STEP 5: Allocate memory and read genome data
		// ========================================================================
		// Allocate memory for genome data arrays for each chromosome
		// Memory layout: compressed format with 2 bits per genotype (4 per byte)
		// Size = (number of bytes per individual) * (number of individuals)
		for (int chrtemp1 = 1; chrtemp1 < 23; chrtemp1++)
		{
			// Calculate bytes needed per individual: (SNPs/4) rounded up
			unsigned long long bytes_per_indiv = (nbsnpperchr[chrtemp1] / 4) + ((nbsnpperchr[chrtemp1] % 4) > 0);
			genomes[chrtemp1] = (unsigned char*)calloc((unsigned long long)(1 + bytes_per_indiv) * NbIndiv, sizeof(char));
		}
		
		// Read input files based on specified format
		printf("Reading genome data files...\n");
		int chromosomes_read = 0;
		if (inputFormat == 0)  // HAP format (SNP-centric)
		{
			printf("  Format: HAP (SNP-centric)\n");
			// HAP files: <PathInput>1.hap, <PathInput>2.hap, ..., <PathInput>22.hap
			for (int chrtemp1 = 1; chrtemp1 < 23; chrtemp1++)
			{
				printf("  [%2d/22] Reading chromosome %2d...", chrtemp1, chrtemp1);
				fflush(stdout);
				// Read HAP file for this chromosome
				// Parameters: run=100+chr, step=105 (used for file naming in some contexts)
				if (readgenomelocal(PathInput, chrtemp1, 100 + chrtemp1, 105, genomes[chrtemp1]) == 0)
				{
					// NOTE: Do NOT update nbsnpperchr[chrtemp1] here!
					// nbsnpperchr[chrtemp1] represents the allocated memory size (used for offset calculations)
					// nbsnpperchrinfile[chrtemp1] represents the actual number of SNPs read (used for loops)
					chromosomes_read++;
					printf(" OK (%d SNPs, allocated for %d)\n", nbsnpperchrinfile[chrtemp1], nbsnpperchr[chrtemp1]);
				}
				else
				{
					printf(" FAILED\n");
				}
			}
		}
		else  // PED format (individual-centric)
		{
			printf("  Format: PED (individual-centric)\n");
			// PED files: <PathInput>1.ped, <PathInput>2.ped, ..., <PathInput>22.ped
			// For PED format, we need to read each chromosome file separately
			// PathInput should be the base path (e.g., "data/chr")
			for (int chrtemp1 = 1; chrtemp1 < 23; chrtemp1++)
			{
				printf("  [%2d/22] Reading chromosome %2d...", chrtemp1, chrtemp1);
				fflush(stdout);
				// Construct filename: <PathInput><chr>.ped
				char pedfilename[300];
				strcpy(pedfilename, PathInput);
				char number[100];
				sprintf(number, "%d", chrtemp1);
				strcat(pedfilename, number);
				strcat(pedfilename, ".ped");
				
				// Read PED file for this chromosome
				if (readpedfile_local(pedfilename, chrtemp1) == 0)
				{
					// NOTE: Do NOT update nbsnpperchr[chrtemp1] here!
					// nbsnpperchr[chrtemp1] represents the allocated memory size (used for offset calculations)
					// nbsnpperchrinfile[chrtemp1] represents the actual number of SNPs read (used for loops)
					chromosomes_read++;
					printf(" OK (%d SNPs, allocated for %d)\n", nbsnpperchrinfile[chrtemp1], nbsnpperchr[chrtemp1]);
				}
				else
				{
					printf(" FAILED\n");
				}
			}
		}
		printf("  Successfully read %d/22 chromosomes\n", chromosomes_read);
		if (chromosomes_read < 22)
		{
			printf("  WARNING: Some chromosomes failed to load. Results may be incomplete.\n");
		}
		printf("\n");
	
		// ========================================================================
		// STEP 6: Process individuals (main phasing loop)
		// ========================================================================
		printf("Starting phasing process...\n");
		printf("================================================================================\n");
		int total_individuals = useListIndiv ? nbIndivToProcess : NbIndiv;
		int processed_count = 0;
		
		for (int idx = 0; idx < nbIndivToProcess; idx++)
		{
			int indiv = listIndivToProcess[idx];
			processed_count++;
			printf("\n[%d/%d] Processing individual %d\n", processed_count, total_individuals, indiv);
			printf("  Calculating PI-HAT scores and identifying relatives...");
			fflush(stdout);
			
			// Step 6a: Load segments and identify relatives
			// This calculates PI-HAT scores and identifies the most closely related individuals
			findrelative(indiv,
				indiv,      // numtrio: trio identifier (same as individual ID)
				indiv,      // IDp1loop: parent 1 ID for loop (same as ID if no parents)
				indiv,      // IDp2loop: parent 2 ID for loop (same as ID if no parents)
				0,          // lenminseg: minimum segment length (0 = no minimum)
				2,          // version: algorithm version
				0,          // gentostart: starting genotype index
				PathInput); // pathresult: base path for result files
			printf(" Done\n");
			
			// Step 6b: Determine phasing method based on parent information
			int parent1 = -1;
			int parent2 = -1;
			if (strlen(PathParentInfo) > 0)
			{
				// Look up parent information for this individual
				getParentInfo(indiv, &parent1, &parent2);
			}
			
			// Step 6c: Execute appropriate phasing algorithm
			if (parent1 == -1 || parent1 == indiv)
			{
				// Use phasing without parents if:
				// - Parents are unknown (-1)
				// - Parent IDs are invalid (same as individual ID)
				if (parent1 == -1 || parent2 == -1)
				{
					printf("  Method: Phasing without parents (parents unknown)\n");
				}
				else
				{
					printf("  Method: Phasing without parents (invalid parent IDs: %d, %d)\n",
						parent1, parent2);
				}
				printf("  Running phasing algorithm...");
				fflush(stdout);
				// Execute phasing without parent information
				predict_phasing_without_parents(indiv,
					chrtodo1,   // chr1: first chromosome (1)
					0,          // chr2: last chromosome (0 = all)
					indiv,      // numtrio: trio identifier
					0,          // limitnbsnp: max SNPs (0 = no limit)
					0,          // run: run identifier
					indiv,      // IDp1loop: parent 1 ID (same as ID when no parents)
					indiv,      // IDp2loop: parent 2 ID (same as ID when no parents)
					0,          // lenminseg: minimum segment length
					2,          // version: algorithm version
					0,          // gentostart: starting genotype index
					PathOutput); // pathresult: base path for output files
				printf(" Done\n");
				
			}
			else
			{
				// Use phasing with parents if both parents are known and valid
				printf("  Method: Phasing with parents (parent1=%d, parent2=%d)\n",
					parent1, parent2);
				printf("  Running phasing algorithm...");
				fflush(stdout);
				// Execute phasing with parent information
				predict_phasing_with_parents_providing_GT(indiv,
					chrtodo1,   // chr1: first chromosome (1)
					0,          // chr2: last chromosome (0 = all)
					indiv,      // numtrio: trio identifier
					0,          // limitnbsnp: max SNPs (0 = no limit)
					0,          // run: run identifier
					parent1,    // IDp1loop: actual parent 1 ID
					parent2,    // IDp2loop: actual parent 2 ID
					0,          // lenminseg: minimum segment length
					2,          // version: algorithm version
					0,          // gentostart: starting genotype index
					PathOutput); // pathresult: base path for output files
				printf(" Done\n");
			}
			
			// Step 6d: Write output files for all chromosomes
			printf("  Writing output files...");
			fflush(stdout);
			for (int chrtemp1 = 1; chrtemp1 < 23; chrtemp1++)
			{
				writeoutput(chrtemp1, PathOutput);
			}
			printf(" Done\n");
			printf("  Individual %d: Complete\n", indiv);
		};
		
		// ========================================================================
		// STEP 7: Program completed successfully
		// ========================================================================
		clock_t end_time = clock();
		float total_elapsed = (float)(end_time - begin0) / CLOCKS_PER_SEC;
		
		printf("\n");
		printf("================================================================================\n");
		printf("                    PHASING COMPLETED SUCCESSFULLY\n");
		printf("================================================================================\n");
		printf("\n");
		printf("Summary:\n");
		printf("  - Individuals processed: %d\n", processed_count);
		printf("  - Chromosomes processed: 22\n");
		printf("  - Output format: %s\n", outputFormat == 0 ? "HAP" : "PED");
		printf("  - Output location: %s\n", PathOutput);
		printf("  - Total execution time: %.2f seconds (%.2f minutes)\n", 
			total_elapsed, total_elapsed / 60.0);
		printf("\n");
		printf("Output files written:\n");
		for (int chr = 1; chr <= 22; chr++)
		{
			char outfile[300];
			strcpy(outfile, PathOutput);
			char number[100];
			sprintf(number, "%d", chr);
			strcat(outfile, number);
			if (outputFormat == 1)
				strcat(outfile, ".ped");
			else
				strcat(outfile, ".hap");
			printf("  - %s\n", outfile);
		}
		printf("\n");
		printf("================================================================================\n");
		printf("\n");
		
		return 0;
	}


