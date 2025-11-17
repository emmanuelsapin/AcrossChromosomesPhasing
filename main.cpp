/**
 * @file main.cpp
 * @brief Main entry point for the Across Chromosomes Phasing program
 * 
 * This file provides the main() function and wrapper functions that maintain
 * backward compatibility with the original C-style API. All functions delegate
 * to a global PhasingProgram instance.
 * 
 * The wrapper functions allow existing code that calls these functions directly
 * to continue working while the implementation has been refactored into an
 * object-oriented structure.
 */

#include "phasing_program.h"

/**
 * @brief Global instance of the PhasingProgram class
 * 
 * This global instance is used by all wrapper functions to maintain
 * backward compatibility with the original C-style API. All wrapper
 * functions delegate their calls to methods of this instance.
 */
PhasingProgram g_program;

/**
 * @brief Wrapper function to read genome data from HAP format files
 * 
 * This function maintains backward compatibility with code that calls
 * readgenomelocal() directly. It delegates to the PhasingProgram::readgenomelocal()
 * method.
 * 
 * @param pathfile Base path for input HAP files (e.g., "data/chr")
 * @param chr Chromosome number (1-22)
 * @param run Run identifier (used for file naming, typically 100+chr)
 * @param step Step identifier (used for file naming, typically 105)
 * @param gentomodify Pointer to memory buffer where genome data will be stored
 * @return 0 on success, 1 on error
 */
int readgenomelocal(char pathfile[], int chr, int run, int step, unsigned char* gentomodify)
{
	return g_program.readgenomelocal(pathfile, chr, run, step, gentomodify);
}

/**
 * @brief Wrapper function to read parent information from a file
 * 
 * Reads a file containing parent-child relationships. The file format is:
 * Individual_ID Parent1_ID Parent2_ID
 * 
 * @param filename Path to the parent information file
 * @return 0 on success, 1 on error
 */
int readParentInfo(const char* filename)
{
	return g_program.readParentInfo(filename);
}

/**
 * @brief Wrapper function to retrieve parent information for an individual
 * 
 * Looks up the parent IDs for a given individual from the loaded parent
 * information. If the individual is not found or has no known parents,
 * parent1 and parent2 are set to -1.
 * 
 * @param indivID Individual ID to look up
 * @param parent1 Output parameter: Parent 1 ID (or -1 if unknown)
 * @param parent2 Output parameter: Parent 2 ID (or -1 if unknown)
 * @return 0 if parents found, 1 if individual not found or no parents
 */
int getParentInfo(int indivID, int* parent1, int* parent2)
{
	return g_program.getParentInfo(indivID, parent1, parent2);
}

/**
 * @brief Wrapper function to read a list of individuals to process
 * 
 * Reads a file containing a list of individual IDs that should be processed.
 * The file format is one ID per line. This allows selective processing
 * of specific individuals rather than processing all individuals.
 * 
 * @param filename Path to the individual list file
 * @return 0 on success, 1 on error
 */
int readListIndiv(const char* filename)
{
	return g_program.readListIndiv(filename);
}

/**
 * @brief Wrapper function to load chromosome segments for phasing
 * 
 * This function loads and analyzes chromosome segments for a given individual.
 * It identifies shared segments with relatives and prepares data structures
 * needed for the phasing algorithm.
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
int findrelative(int ID, int numtrio, int IDp1loop, int IDp2loop, int lenminseg, int version, int gentostart, char pathresult[])
{
	return g_program.findrelative(ID, numtrio, IDp1loop, IDp2loop, lenminseg, version, gentostart, pathresult);
}

/**
 * @brief Wrapper function for phasing without parent information
 * 
 * Performs haplotype phasing for an individual when parent information
 * is not available. This method uses cross-chromosome correlations and
 * population-level information to determine phase.
 * 
 * @param ID Individual ID to phase
 * @param chr1 First chromosome to process (typically 1)
 * @param chr2 Last chromosome to process (typically 0 for all)
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
int predict_phasing_without_parents(int ID, int chr1, int chr2, int numtrio, float limitnbsnp, int run, int IDp1loop, int IDp2loop, int lenminseg, int version, int gentostart, char pathresult[])
{
	return g_program.predict_phasing_without_parents(ID, chr1, chr2, numtrio, limitnbsnp, run, IDp1loop, IDp2loop, lenminseg, version, gentostart, pathresult);
}

/**
 * @brief Wrapper function for phasing with parent information
 * 
 * Performs haplotype phasing for an individual when parent information
 * is available. This method uses Mendelian inheritance rules combined
 * with cross-chromosome information for more accurate phasing.
 * 
 * @param ID Individual ID to phase
 * @param chr1 First chromosome to process (typically 1)
 * @param chr2 Last chromosome to process (typically 0 for all)
 * @param numtrio Trio identifier
 * @param limitnbsnp Maximum number of SNPs to consider
 * @param run Run identifier
 * @param IDp1loop Parent 1 ID (actual parent ID, not same as ID)
 * @param IDp2loop Parent 2 ID (actual parent ID, not same as ID)
 * @param lenminseg Minimum segment length
 * @param version Algorithm version
 * @param gentostart Starting genotype index
 * @param pathresult Base path for result files
 * @return 0 on success, 1 on error
 */
int predict_phasing_with_parents_providing_GT(int ID, int chr1, int chr2, int numtrio, float limitnbsnp, int run, int IDp1loop, int IDp2loop, int lenminseg, int version, int gentostart, char pathresult[])
{
	return g_program.predict_phasing_with_parents_providing_GT(ID, chr1, chr2, numtrio, limitnbsnp, run, IDp1loop, IDp2loop, lenminseg, version, gentostart, pathresult);
}

/**
 * @brief Wrapper function to write phased output to files
 * 
 * Writes the phased haplotypes for a given chromosome to an output file.
 * The output format is determined by the outputFormat member variable
 * (PED format by default).
 * 
 * @param chr Chromosome number (1-22)
 * @param pathoutput Base path for output files
 * @return 0 on success, 1 on error
 */
int writeoutput(int chr, char pathoutput[])
{
	return g_program.writeoutput(chr, pathoutput);
}

/**
 * @brief Main entry point of the program
 * 
 * This is the main() function that serves as the entry point for the program.
 * It delegates all processing to the PhasingProgram::run() method, which
 * handles command-line argument parsing and program execution.
 * 
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return 0 on successful completion, non-zero on error
 */
int main(int argc, char* argv[])
{
	return g_program.run(argc, argv);
}

