/**
 * @file phasing_program.h
 * @brief Header file for the PhasingProgram class
 * 
 * This file defines the PhasingProgram class, which encapsulates all the
 * functionality of the haplotype phasing program. The class contains:
 * - Data members for storing genome data, relationships, and configuration
 * - Methods for reading input files (HAP, PED, parent info, individual lists)
 * - Methods for performing phasing algorithms
 * - Methods for writing output files
 * - The main run() method that orchestrates the entire phasing process
 */

#ifndef PHASING_PROGRAM_H
#define PHASING_PROGRAM_H

#include "types.h"

/**
 * @class PhasingProgram
 * @brief Main class for haplotype phasing across chromosomes
 * 
 * This class encapsulates all functionality for performing haplotype phasing.
 * It can process individuals with or without known parent information, and
 * supports both HAP and PED file formats for input and output.
 */
class PhasingProgram
{
public:
	// ============================================================================
	// Data Members - Relationship and Phasing Information
	// ============================================================================
	
	/**
	 * @brief Array storing IDs of individuals with best PI-HAT scores
	 * 
	 * PI-HAT is a measure of relatedness between individuals. This array
	 * stores the IDs of the top 100 individuals most related to the current
	 * individual being processed.
	 */
	int bestpihatagainstallID[100];
	
	/**
	 * @brief ID of the individual with the best PI-HAT score
	 * 
	 * This is the ID of the most closely related individual found during
	 * relationship detection. Used for phasing decisions.
	 */
	int IDbestpihat;
	
	/**
	 * @brief ID of the individual with the second-best PI-HAT score
	 * 
	 * Used when multiple relatives are needed for phasing decisions.
	 */
	int IDbestpihat2;
	
	/**
	 * @brief Array storing PI-HAT scores against all individuals
	 * 
	 * PI-HAT (π-hat) is a measure of genetic relatedness. Values range
	 * from 0 (unrelated) to 1 (identical). This array stores the PI-HAT
	 * score for each individual in the population against the current
	 * individual being processed.
	 */
	float pihatagainstall[MAXPOP];
	
	/**
	 * @brief Number of chromosome dividers for each method and chromosome
	 * 
	 * Chromosome dividers split chromosomes into segments for analysis.
	 * First dimension: method index (0-50)
	 * Second dimension: chromosome number (0-22, where 0 is metadata)
	 */
	int nbchrdivider[51][23];
	
	/**
	 * @brief Job identifier for parallel processing
	 * 
	 * Used to track which job/thread is processing which individual
	 * in parallel execution scenarios.
	 */
	int IDjob;
	
	/**
	 * @brief Index of first relative to consider in phasing
	 * 
	 * Used to optimize phasing by starting from the most relevant
	 * relative rather than checking all individuals.
	 */
	int placefirttoconsider;
	
	/**
	 * @brief Chromosome divider information
	 * 
	 * Stores information about how chromosomes are divided into segments
	 * for analysis. First dimension: method index (0-49)
	 * Second dimension: chromosome number (0-22)
	 * Third dimension: divider index (0-19, up to 20 dividers per chromosome)
	 */
	typechrdivider chrdivider[50][23][20];
	
	/**
	 * @brief PI-HAT thresholds for phasing decisions
	 * 
	 * Array of three threshold values used to determine which relatives
	 * to use for phasing. Different thresholds may be used for different
	 * relationship levels or algorithm stages.
	 */
	float seuilpihat[3];
	
	/**
	 * @brief Results matrix for different methods and metrics
	 * 
	 * Stores results for different phasing methods (first dimension: 0-3)
	 * and different evaluation metrics (second dimension: 0-5).
	 * Used for comparing algorithm performance.
	 */
	int result[4][6];
	
	/**
	 * @brief Best PI-HAT scores for different relationship categories
	 * 
	 * Stores the best PI-HAT scores found for different categories
	 * of relationships (e.g., parent-child, siblings, etc.).
	 */
	float bestpihat[7];
	
	/**
	 * @brief Flag to enable detailed output/printing
	 * 
	 * When set to 1, the program prints detailed information during
	 * execution. When 0, only essential messages are printed.
	 */
	int printdetail;
	
	// ============================================================================
	// Data Members - Genome Data and SNPs
	// ============================================================================
	
	/**
	 * @brief Number of SNPs per chromosome (maximum capacity)
	 * 
	 * Array storing the maximum number of SNPs that can be stored
	 * for each chromosome. Index 0 is metadata, indices 1-22 are
	 * chromosomes 1-22.
	 */
	int nbsnpperchr[23];
	
	/**
	 * @brief Number of SNPs per chromosome in input files
	 * 
	 * Array storing the actual number of SNPs read from input files
	 * for each chromosome. This may be less than nbsnpperchr if
	 * the input files contain fewer SNPs.
	 */
	int nbsnpperchrinfile[23];
	
	/**
	 * @brief Minor Allele Frequency (MAF) for each SNP
	 * 
	 * Two-dimensional array storing the MAF for each SNP.
	 * First dimension: SNP index (0 to NSNPPERCHR-1)
	 * Second dimension: chromosome number (0-22)
	 * 
	 * MAF values are used to weight phasing decisions, as rare variants
	 * provide more information for phasing than common variants.
	 */
	int MAF[NSNPPERCHR][23];
	
	/**
	 * @brief Genome data for phased haplotypes
	 * 
	 * Three-dimensional array storing phased genotype data.
	 * First dimension: haplotype copy (0-2, typically 0 for unphased)
	 * Second dimension: SNP index
	 * Third dimension: chromosome number
	 * 
	 * Used to store intermediate phasing results and final phased data.
	 */
	int genomeoffpss[3][NSNPPERCHR][23];
	
	/**
	 * @brief Number of individuals in the dataset
	 * 
	 * Total number of individuals to process. Must be set before
	 * reading genome data.
	 */
	int NbIndiv;
	
	/**
	 * @brief Pointers to genome data arrays for each chromosome
	 * 
	 * Array of pointers, one for each chromosome (0-22).
	 * Each pointer points to a compressed genome data array where
	 * genotypes are stored as 2-bit values (4 genotypes per byte).
	 * 
	 * Memory layout: [individual][SNP/4] where each byte contains
	 * 4 genotypes encoded as 2 bits each.
	 */
	unsigned char* genomes[23];
	
	/**
	 * @brief Decision points for chromosome division
	 * 
	 * Array storing decision points used when dividing chromosomes
	 * into segments for analysis. Each point contains correlation
	 * information and group sizes.
	 */
	pointdecision tappointdec[MAXNBDIVISOR * 22];
	
	// ============================================================================
	// Data Members - Parent and Individual Selection
	// ============================================================================
	
	/**
	 * @brief Parent information for individuals
	 * 
	 * Array storing parent-child relationships. Each entry contains
	 * an individual ID and the IDs of their parents (or -1 if unknown).
	 */
	structparentinfo parentinfo[NBINDIVMAX];
	
	/**
	 * @brief Number of parent information entries loaded
	 * 
	 * Count of how many parent information entries have been loaded
	 * from the parent info file. May be less than NbIndiv if not
	 * all individuals have known parents.
	 */
	int nbparentinfo;
	
	/**
	 * @brief List of individual IDs to process
	 * 
	 * Array storing the IDs of individuals that should be processed.
	 * Used when -ListIndiv option is provided to process only a
	 * subset of individuals.
	 */
	int listIndivToProcess[NBINDIVMAX];
	
	/**
	 * @brief Number of individuals in the processing list
	 * 
	 * Count of how many individual IDs are in listIndivToProcess.
	 * If 0, all individuals are processed.
	 */
	int nbIndivToProcess;
	
	/**
	 * @brief Flag indicating whether to use the individual list
	 * 
	 * When 1, only individuals in listIndivToProcess are processed.
	 * When 0, all individuals are processed.
	 */
	int useListIndiv;
	
	// ============================================================================
	// Data Members - Algorithm Parameters and State
	// ============================================================================
	
	/**
	 * @brief Number of breakpoints in chromosome segments
	 * 
	 * Tracks the number of breakpoints (recombination sites) identified
	 * in chromosome segments during phasing analysis.
	 */
	int nbbreak;
	
	/**
	 * @brief Distribution of gametic states
	 * 
	 * Array storing counts of different gametic states (0-11).
	 * Used for tracking inheritance patterns and phasing decisions.
	 */
	int distrigametic[12];
	
	/**
	 * @brief Distribution of gametic states to keep
	 * 
	 * Array storing which gametic states should be retained for
	 * each chromosome during phasing.
	 */
	int distrigametickeep[23];
	
	/**
	 * @brief Array marking which individuals are relatives
	 * 
	 * Used to track which individuals in the population are identified
	 * as relatives of the current individual being processed.
	 */
	int checkrelat[MAXPOP];
	
	/**
	 * @brief Start positions of chromosome segments
	 * 
	 * Array storing the start positions of segments for each chromosome
	 * (indices 0-23, where 0 is metadata and 1-22 are chromosomes).
	 */
	int segstart[24];
	
	/**
	 * @brief Number of relatives with PI-HAT information
	 * 
	 * Count of how many relatives have been identified and have
	 * PI-HAT scores calculated.
	 */
	int nbrelatpihat;
	
	// ============================================================================
	// Data Members - File Format Configuration
	// ============================================================================
	
	/**
	 * @brief Input file format flag
	 * 
	 * 0 = HAP format (SNP-centric, one SNP per line)
	 * 1 = PED format (individual-centric, one individual per line)
	 */
	int inputFormat;
	
	/**
	 * @brief Output file format flag
	 * 
	 * 0 = HAP format (not fully implemented, writes as PED)
	 * 1 = PED format (default, individual-centric output)
	 */
	int outputFormat;
	
	// ============================================================================
	// Constructor
	// ============================================================================
	
	/**
	 * @brief Constructor for PhasingProgram
	 * 
	 * Initializes all member variables to default values:
	 * - Sets input format to HAP (0)
	 * - Sets output format to PED (1)
	 * - Initializes arrays to zero/null
	 * - Allocates memory for genome pointers (set to nullptr)
	 */
	PhasingProgram();
	
	// ============================================================================
	// Input/Output Methods
	// ============================================================================
	
	/**
	 * @brief Read genome data from HAP format file
	 * 
	 * Reads genotype data from a HAP format file for a specific chromosome.
	 * HAP format is SNP-centric: each line represents one SNP with genotypes
	 * for all individuals.
	 * 
	 * @param pathfile Base path for input files (e.g., "data/chr")
	 * @param chr Chromosome number (1-22)
	 * @param run Run identifier (typically 100+chr for file naming)
	 * @param step Step identifier (typically 105 for file naming)
	 * @param gentomodify Pointer to memory buffer where data will be stored
	 * @return 0 on success, 1 on error
	 */
	int readgenomelocal(char pathfile[], int chr, int run, int step, unsigned char* gentomodify);
	
	/**
	 * @brief Read genome data from PED format file
	 * 
	 * Reads genotype data from a PED format file for a specific chromosome.
	 * PED format is individual-centric: each line represents one individual
	 * with genotypes for all SNPs.
	 * 
	 * @param filename Full path to the PED file (e.g., "data/chr1.ped")
	 * @param chr Chromosome number (1-22)
	 * @return 0 on success, 1 on error
	 */
	int readpedfile_local(const char* filename, int chr);
	
	/**
	 * @brief Read parent information from file
	 * 
	 * Reads a file containing parent-child relationships. File format:
	 * Individual_ID Parent1_ID Parent2_ID
	 * Use -1 for unknown parents.
	 * 
	 * @param filename Path to parent information file
	 * @return 0 on success, 1 on error
	 */
	int readParentInfo(const char* filename);
	
	/**
	 * @brief Get parent information for an individual
	 * 
	 * Looks up the parent IDs for a given individual from the loaded
	 * parent information.
	 * 
	 * @param indivID Individual ID to look up
	 * @param parent1 Output: Parent 1 ID (or -1 if unknown)
	 * @param parent2 Output: Parent 2 ID (or -1 if unknown)
	 * @return 0 if found, 1 if not found
	 */
	int getParentInfo(int indivID, int* parent1, int* parent2);
	
	/**
	 * @brief Read list of individuals to process
	 * 
	 * Reads a file containing a list of individual IDs that should be
	 * processed. File format: one ID per line.
	 * 
	 * @param filename Path to individual list file
	 * @return 0 on success, 1 on error
	 */
	int readListIndiv(const char* filename);
	
	/**
	 * @brief Write phased output to file
	 * 
	 * Writes the phased haplotypes for a chromosome to an output file.
	 * The format is determined by outputFormat member variable.
	 * 
	 * @param chr Chromosome number (1-22)
	 * @param pathoutput Base path for output files
	 * @return 0 on success, 1 on error
	 */
	int writeoutput(int chr, char pathoutput[]);
	
	// ============================================================================
	// Phasing Algorithm Methods
	// ============================================================================
	
	/**
	 * @brief Find and identify relatives
	 * 
	 * Finds relatives for an individual by analyzing chromosome segments
	 * and calculating PI-HAT (relatedness coefficient). This is a preparatory step before phasing.
	 * 
	 * @param ID Individual ID to process
	 * @param numtrio Trio identifier (typically same as ID)
	 * @param IDp1loop Parent 1 ID for loop processing
	 * @param IDp2loop Parent 2 ID for loop processing
	 * @param lenminseg Minimum segment length to consider
	 * @param version Algorithm version to use
	 * @param gentostart Starting genotype index
	 * @param pathresult Base path for result files
	 * @return 0 on success, 1 on error
	 */
	int findrelative(int ID, int numtrio, int IDp1loop, int IDp2loop, int lenminseg, int version, int gentostart, char pathresult[]);
	
	/**
	 * @brief Perform phasing without parent information
	 * 
	 * Phases haplotypes for an individual when parent information is not
	 * available. Uses cross-chromosome correlations and population-level
	 * information to determine phase.
	 * 
	 * @param ID Individual ID to phase
	 * @param chr1 First chromosome to process (typically 1)
	 * @param chr2 Last chromosome to process (0 for all)
	 * @param numtrio Trio identifier
	 * @param limitnbsnp Maximum number of SNPs to consider
	 * @param run Run identifier
	 * @param IDp1loop Parent 1 ID for loop (typically same as ID)
	 * @param IDp2loop Parent 2 ID for loop (typically same as ID)
	 * @param lenminseg Minimum segment length
	 * @param version Algorithm version
	 * @param gentostart Starting genotype index
	 * @param pathresult Base path for result files
	 * @return 0 on success, 1 on error
	 */
	int predict_phasing_without_parents(int ID, int chr1, int chr2, int numtrio, float limitnbsnp, int run, int IDp1loop, int IDp2loop, int lenminseg, int version, int gentostart, char pathresult[]);
	
	/**
	 * @brief Perform phasing with parent information
	 * 
	 * Phases haplotypes for an individual when parent information is
	 * available. Uses Mendelian inheritance rules combined with
	 * cross-chromosome information for more accurate phasing.
	 * 
	 * @param ID Individual ID to phase
	 * @param chr1 First chromosome to process (typically 1)
	 * @param chr2 Last chromosome to process (0 for all)
	 * @param numtrio Trio identifier
	 * @param limitnbsnp Maximum number of SNPs to consider
	 * @param run Run identifier
	 * @param IDp1loop Parent 1 ID (actual parent, not same as ID)
	 * @param IDp2loop Parent 2 ID (actual parent, not same as ID)
	 * @param lenminseg Minimum segment length
	 * @param version Algorithm version
	 * @param gentostart Starting genotype index
	 * @param pathresult Base path for result files
	 * @return 0 on success, 1 on error
	 */
	int predict_phasing_with_parents_providing_GT(int ID, int chr1, int chr2, int numtrio, float limitnbsnp, int run, int IDp1loop, int IDp2loop, int lenminseg, int version, int gentostart, char pathresult[]);
	
	// ============================================================================
	// Main Execution Method
	// ============================================================================
	
	/**
	 * @brief Main execution method
	 * 
	 * This is the main method that orchestrates the entire phasing process:
	 * 1. Parses command-line arguments
	 * 2. Reads input files (genome data, parent info, individual lists)
	 * 3. For each individual to process:
	 *    - Finds relatives and calculates PI-HAT scores
	 *    - Determines phasing method (with/without parents)
	 *    - Executes phasing algorithm
	 *    - Writes output files
	 * 
	 * @param argc Number of command-line arguments
	 * @param argv Array of command-line argument strings
	 * @return 0 on success, non-zero on error
	 */
	int run(int argc, char* argv[]);
};

#endif

