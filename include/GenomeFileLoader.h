/**
 * @file GenomeFileLoader.h
 * @brief Factory pattern for loading genome data from various formats
 */

#ifndef GENOME_FILE_LOADER_H
#define GENOME_FILE_LOADER_H

#include "utils/readinteger.h"

namespace PhasingEngine {

/**
 * @class GenomeFileLoader
 * @brief Factory pattern for loading genome data from various formats
 */
class GenomeFileLoader {
public:
    enum class FileFormat {
        HAP_FORMAT,
        PED_FORMAT,
        VCF_FORMAT
    };
    
    static bool loadGenome(const char* filePath, int chromosome, 
                          unsigned char* genomeBuffer, int numberOfIndividuals,
                          int* snpCountPerChr, int* snpCountInFile,
                          FileFormat format = FileFormat::HAP_FORMAT);
    
    static bool validateFileFormat(const char* filePath, FileFormat format);
    
private:
    static bool loadHAPFormat(const char* filePath, int chromosome,
                             unsigned char* genomeBuffer, int numberOfIndividuals,
                             int* snpCountPerChr, int* snpCountInFile);
};

}

#endif // GENOME_FILE_LOADER_H

