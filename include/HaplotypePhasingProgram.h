/**
 * @file HaplotypePhasingProgram.h
 * @brief Main application class orchestrating the phasing pipeline
 */

#ifndef HAPLOTYPE_PHASING_PROGRAM_H
#define HAPLOTYPE_PHASING_PROGRAM_H

#include <memory>
#include "ConfigurationManager.h"

namespace PhasingEngine {

class GenomeDataManager;
class RelativeIdentificationEngine;
class PhasingAlgorithmEngine;
class OutputFileWriter;

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
    
private:
    bool processIndividuals();
};

}

#endif // HAPLOTYPE_PHASING_PROGRAM_H


