/**
 * @file PhasingProgram.h
 * @brief Haplotype Phasing Engine - Professional Implementation
 * @details This module provides a comprehensive object-oriented framework for 
 *          genetic haplotype phasing using relative-based inference algorithms.
 * @version 2.0
 * @date 2024
 */

#ifndef PHASING_PROGRAM_H
#define PHASING_PROGRAM_H

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <memory>
#include <string>
#include <stdexcept>
#include <cfloat>
#include <climits>
#include <omp.h>
#include <algorithm>
#include <functional>

#include "readinteger.h"
#include "readreal.h"
#include "readnegativereal.h"

/**
 * @namespace PhasingEngine
 * @brief Main namespace for haplotype phasing functionality
 */
namespace PhasingEngine {

/**
 * @namespace Constants
 * @brief System-wide constants and configuration parameters
 */
namespace Constants {
    constexpr int NSNPPERCHR = 60000;
    constexpr int MAXPOP = 435188;
    constexpr int MAXCLOSERELAT = 6000;
    constexpr int MAXCLOSERELATTEMP = 450000;
    constexpr int MAXNBDIVISOR = 25;
    constexpr int NBINDIVMAX = 100000;
    constexpr int NBINDIV = 100000;
    constexpr int NBINDIVEA = 1;
    constexpr int MAXGEN = 1;
    constexpr int MAXBREAK = 1000;
    constexpr int NUM_CHROMOSOMES = 23;
    constexpr float DEFAULT_PIHAT_THRESHOLD = 0.33f;
    constexpr float PIHAT_NORMALIZATION_FACTOR = 330005.0f;
}

/**
 * @namespace ErrorCodes
 * @brief Error code enumeration for exception handling
 */
namespace ErrorCodes {
    enum class PhasingError {
        SUCCESS = 0,
        FILE_NOT_FOUND = 1,
        INVALID_INPUT = 2,
        MEMORY_ALLOCATION_FAILED = 3,
        INVALID_CHROMOSOME = 4,
        INVALID_INDIVIDUAL = 5,
        INVALID_SNP_INDEX = 6,
        DATA_CORRUPTION = 7
    };
}

/**
 * @class PhasingException
 * @brief Custom exception class for phasing operations
 */
class PhasingException : public std::runtime_error {
private:
    ErrorCodes::PhasingError errorCode;
    
public:
    PhasingException(const std::string& message, ErrorCodes::PhasingError code = ErrorCodes::PhasingError::INVALID_INPUT)
        : std::runtime_error(message), errorCode(code) {}
    
    ErrorCodes::PhasingError getErrorCode() const { return errorCode; }
};

/**
 * @class IGenomeDataAccessor
 * @brief Abstract interface for genome data access operations
 */
class IGenomeDataAccessor {
public:
    virtual ~IGenomeDataAccessor() = default;
    virtual int getGenotype(int chromosome, int individual, int snpIndex) const = 0;
    virtual void setGenotype(int chromosome, int individual, int snpIndex, int genotype) = 0;
    virtual int getSNPCount(int chromosome) const = 0;
    virtual bool isValidChromosome(int chromosome) const = 0;
    virtual bool isValidIndividual(int individual) const = 0;
};

/**
 * @class IRelativeFinder
 * @brief Abstract interface for relative identification algorithms
 */
class IRelativeFinder {
public:
    virtual ~IRelativeFinder() = default;
    virtual void computePIHATMatrix() = 0;
    virtual void identifyBestRelatives(float threshold) = 0;
    virtual int getBestRelativeID(int rank) const = 0;
    virtual float getPIHATValue(int individualID) const = 0;
};

/**
 * @class IPhasingStrategy
 * @brief Strategy pattern interface for different phasing algorithms
 */
class IPhasingStrategy {
public:
    virtual ~IPhasingStrategy() = default;
    virtual int executePhasing(int individualID, int chromosomeStart, int chromosomeEnd) = 0;
    virtual std::string getStrategyName() const = 0;
};

/**
 * @struct ChromosomeDivider
 * @brief Represents a chromosomal segment with phasing information
 */
struct ChromosomeDivider {
    int start;          ///< Starting SNP index
    int end;            ///< Ending SNP index
    int phasing;        ///< Phasing orientation (1 or 2)
    int nbright;        ///< Number of correctly phased positions
    int nbwrong;        ///< Number of incorrectly phased positions
    int segment;        ///< Segment identifier
    
    ChromosomeDivider() : start(0), end(0), phasing(0), nbright(0), nbwrong(0), segment(0) {}
    
    /**
     * @brief Calculate phasing accuracy percentage
     * @return Accuracy as percentage (0-100)
     */
    double calculateAccuracy() const {
        int total = nbright + nbwrong;
        return total > 0 ? (100.0 * nbright) / total : 0.0;
    }
    
    /**
     * @brief Check if segment is valid
     * @return true if start < end
     */
    bool isValid() const { return start < end && start >= 0; }
    
    /**
     * @brief Get segment length in SNPs
     * @return Number of SNPs in segment
     */
    int getLength() const { return end - start; }
};

/**
 * @class GenomeDataManager
 * @brief Comprehensive genome data management with validation and caching
 * @implements IGenomeDataAccessor
 */
class GenomeDataManager : public IGenomeDataAccessor {
private:
    unsigned char* genomes[Constants::NUM_CHROMOSOMES];
    int snpCountPerChromosome[Constants::NUM_CHROMOSOMES];
    int snpCountInFile[Constants::NUM_CHROMOSOMES];
    int numberOfIndividuals;
    int minorAlleleFrequency[Constants::NSNPPERCHR][Constants::NUM_CHROMOSOMES];
    int genomeOffspringData[3][Constants::NSNPPERCHR][Constants::NUM_CHROMOSOMES];
    bool isInitialized;
    
    void validateChromosomeIndex(int chromosome) const;
    void validateIndividualIndex(int individual) const;
    void validateSNPIndex(int chromosome, int snpIndex) const;
    size_t calculateGenomeBufferSize(int chromosome) const;
    
public:
    GenomeDataManager();
    virtual ~GenomeDataManager();
    
    // IGenomeDataAccessor interface
    virtual int getGenotype(int chromosome, int individual, int snpIndex) const override;
    virtual void setGenotype(int chromosome, int individual, int snpIndex, int genotype) override;
    virtual int getSNPCount(int chromosome) const override;
    virtual bool isValidChromosome(int chromosome) const override;
    virtual bool isValidIndividual(int individual) const override;
    
    // Extended interface
    bool initializeChromosome(int chromosome, int snpCount, int nbIndiv);
    bool loadFromFile(const char* pathfile, int chromosome, int nbIndiv);
    void setGenomeBuffer(int chromosome, unsigned char* buffer);
    unsigned char* getGenomeBuffer(int chromosome) const;
    
    int getSNPCountPerChr(int chromosome) const;
    void setSNPCountPerChr(int chromosome, int count);
    void setSNPCountInFile(int chromosome, int count);
    int* getSNPCountPerChrArray() { return snpCountPerChromosome; }
    int* getSNPCountInFileArray() { return snpCountInFile; }
    
    int getMAF(int snpIndex, int chromosome) const;
    void setMAF(int snpIndex, int chromosome, int value);
    void computeMAFForChromosome(int chromosome);
    
    int getGenomeOffspring(int index, int snpIndex, int chromosome) const;
    void setGenomeOffspring(int index, int snpIndex, int chromosome, int value);
    
    void setNumberOfIndividuals(int count);
    int getNumberOfIndividuals() const { return numberOfIndividuals; }
    
    bool isDataInitialized() const { return isInitialized; }
    void reset();
};

/**
 * @class RelativeIdentificationEngine
 * @brief Advanced relative identification using PIHAT computation
 * @implements IRelativeFinder
 */
class RelativeIdentificationEngine : public IRelativeFinder {
private:
    float pihatMatrix[Constants::MAXPOP];
    float pihatMatrixSecondary[Constants::MAXPOP];
    int bestRelativeIDs[100];
    int primaryBestRelativeID;
    int secondaryBestRelativeID;
    int firstConsiderationIndex;
    bool isComputed;
    
    void sortRelativesByPIHAT();
    void updateBestRelativeRankings();
    
public:
    RelativeIdentificationEngine();
    virtual ~RelativeIdentificationEngine() = default;
    
    // IRelativeFinder interface
    virtual void computePIHATMatrix() override;
    virtual void identifyBestRelatives(float threshold) override;
    virtual int getBestRelativeID(int rank) const override;
    virtual float getPIHATValue(int individualID) const override;
    
    // Extended interface
    void reset();
    void accumulatePIHAT(int relativeID, float contribution);
    void setPIHAT2(int relativeID, float value);
    float getPIHAT2(int relativeID) const;
    
    int getPrimaryBestRelativeID() const { return primaryBestRelativeID; }
    int getSecondaryBestRelativeID() const { return secondaryBestRelativeID; }
    int getFirstConsiderationIndex() const { return firstConsiderationIndex; }
    
    int getRelativeCountAboveThreshold(float threshold) const;
    bool hasComputed() const { return isComputed; }
};

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

/**
 * @class PhasingAlgorithmEngine
 * @brief Core phasing algorithm implementation with multiple strategies
 */
class PhasingAlgorithmEngine {
private:
    GenomeDataManager* genomeDataManager;
    RelativeIdentificationEngine* relativeEngine;
    ChromosomeDivider chromosomeDividers[50][Constants::NUM_CHROMOSOMES][20];
    int chromosomeDividerCounts[51][Constants::NUM_CHROMOSOMES];
    float pihatThresholds[3];
    bool verboseOutput;
    int jobIdentifier;
    int relativeCount;
    int breakpointCount;
    
    std::unique_ptr<IPhasingStrategy> currentStrategy;
    
    void initializeChromosomeDividers();
    void loadGenomeOffspringData(int individualID, int parent1ID, int parent2ID);
    void processPhasingCorrections(int individualID, int relativeID);
    void mergeChromosomeWindows(int breakpointIndex);
    void applyPhasingToGenome(int individualID, int breakpointIndex);
    void computeCorrelationMatrix(int chromosome1, int chromosome2, 
                                 double* correlationMatrix) const;
    void optimizeWindowMerging(int breakpointIndex);
    
public:
    PhasingAlgorithmEngine(GenomeDataManager* gdm, RelativeIdentificationEngine* rie);
    virtual ~PhasingAlgorithmEngine();
    
    bool loadSegment(int individualID, int trioNumber, int parent1ID, int parent2ID,
                    int minimumSegmentLength, int algorithmVersion, int generationStart,
                    const char* resultPath);
    
    int executePhasingAlgorithm(int individualID, int chromosomeStart, int chromosomeEnd,
                               int trioNumber, float snpLimit, int runNumber,
                               int parent1ID, int parent2ID, int minimumSegmentLength,
                               int algorithmVersion, int generationStart, const char* resultPath);
    
    void setPhasingStrategy(std::unique_ptr<IPhasingStrategy> strategy);
    void setVerboseOutput(bool enabled) { verboseOutput = enabled; }
    bool isVerboseOutput() const { return verboseOutput; }
    
    ChromosomeDivider* getChromosomeDivider(int sizeIndex, int chromosome, int window);
    int getChromosomeDividerCount(int sizeIndex, int chromosome) const;
    void setChromosomeDividerCount(int sizeIndex, int chromosome, int count);
    
    void setRelativeCount(int count) { relativeCount = count; }
    int getRelativeCount() const { return relativeCount; }
    void setBreakpointCount(int count) { breakpointCount = count; }
    int getBreakpointCount() const { return breakpointCount; }
    
    float getPIHATThreshold(int index) const;
    void setPIHATThreshold(int index, float value);
};

/**
 * @class OutputFileWriter
 * @brief Professional output generation with multiple format support
 */
class OutputFileWriter {
private:
    GenomeDataManager* genomeDataManager;
    std::string outputDirectory;
    bool validateOutputPath(const char* path) const;
    void writePEDHeader(FILE* file, int chromosome) const;
    void writePEDRow(FILE* file, int individualID, int chromosome) const;
    
public:
    explicit OutputFileWriter(GenomeDataManager* gdm);
    virtual ~OutputFileWriter() = default;
    
    bool writeOutput(int chromosome, const char* outputPath);
    bool writeOutputBatch(const std::vector<int>& chromosomes, const char* outputPath);
    void setOutputDirectory(const std::string& directory) { outputDirectory = directory; }
    std::string getOutputDirectory() const { return outputDirectory; }
};

/**
 * @class ConfigurationManager
 * @brief Centralized configuration management
 */
class ConfigurationManager {
private:
    int numberOfIndividuals;
    std::string inputPath;
    std::string outputPath;
    std::string mafFilePath;
    bool verboseMode;
    int algorithmVersion;
    float pihatThreshold;
    
public:
    ConfigurationManager();
    
    bool parseCommandLineArguments(int argc, char* argv[]);
    bool validateConfiguration() const;
    
    int getNumberOfIndividuals() const { return numberOfIndividuals; }
    void setNumberOfIndividuals(int count) { numberOfIndividuals = count; }
    
    std::string getInputPath() const { return inputPath; }
    void setInputPath(const std::string& path) { inputPath = path; }
    
    std::string getOutputPath() const { return outputPath; }
    void setOutputPath(const std::string& path) { outputPath = path; }
    
    std::string getMAFFilePath() const { return mafFilePath; }
    void setMAFFilePath(const std::string& path) { mafFilePath = path; }
    
    bool isVerboseMode() const { return verboseMode; }
    void setVerboseMode(bool enabled) { verboseMode = enabled; }
    
    int getAlgorithmVersion() const { return algorithmVersion; }
    void setAlgorithmVersion(int version) { algorithmVersion = version; }
    
    float getPIHATThreshold() const { return pihatThreshold; }
    void setPIHATThreshold(float threshold) { pihatThreshold = threshold; }
};

/**
 * @class HaplotypePhasingProgram
 * @brief Main application class orchestrating the phasing pipeline
 */
class HaplotypePhasingProgram {
private:
    std::unique_ptr<GenomeDataManager> genomeDataManager;
    std::unique_ptr<RelativeIdentificationEngine> relativeEngine;
    std::unique_ptr<PhasingAlgorithmEngine> phasingEngine;
    std::unique_ptr<OutputFileWriter> outputWriter;
    std::unique_ptr<ConfigurationManager> configuration;
    
    void initializeSNPCounts();
    void initializeChromosomeDividers();
    bool validateInputFiles() const;
    void logExecutionStatistics(clock_t startTime, clock_t endTime) const;
    
public:
    HaplotypePhasingProgram();
    virtual ~HaplotypePhasingProgram();
    
    bool initialize(int argc, char* argv[]);
    bool execute();
    void shutdown();
    
    ConfigurationManager* getConfiguration() const { return configuration.get(); }
    GenomeDataManager* getGenomeDataManager() const { return genomeDataManager.get(); }
};

} // namespace PhasingEngine

#endif // PHASING_PROGRAM_H
