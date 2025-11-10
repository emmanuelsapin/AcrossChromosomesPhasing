/**
 * @file HaplotypePhasingProgram.cpp
 * @brief Implementation of HaplotypePhasingProgram
 */

#include "../include/HaplotypePhasingProgram.h"
#include "../include/GenomeDataManager.h"
#include "../include/RelativeIdentificationEngine.h"
#include "../include/PhasingAlgorithmEngine.h"
#include "../include/OutputFileWriter.h"
#include "../include/ConfigurationManager.h"
#include "../include/GenomeFileLoader.h"
#include "../include/Constants.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <memory>

using namespace PhasingEngine;
using namespace PhasingEngine::Constants;

HaplotypePhasingProgram::HaplotypePhasingProgram() {
    genomeDataManager = std::make_unique<GenomeDataManager>();
    relativeEngine = std::make_unique<RelativeIdentificationEngine>();
    phasingEngine = std::make_unique<PhasingAlgorithmEngine>(
        genomeDataManager.get(), relativeEngine.get());
    outputWriter = std::make_unique<OutputFileWriter>(genomeDataManager.get());
    configuration = std::make_unique<ConfigurationManager>();
}

HaplotypePhasingProgram::~HaplotypePhasingProgram() {
    shutdown();
}

bool HaplotypePhasingProgram::initialize(int argc, char* argv[]) {
    if(!configuration->parseCommandLineArguments(argc, argv)) {
        return false;
    }
    
    genomeDataManager->setNumberOfIndividuals(configuration->getNumberOfIndividuals());
    phasingEngine->setVerboseOutput(configuration->isVerboseMode());
    
    return true;
}

void HaplotypePhasingProgram::initializeSNPCounts() {
    const int snpCounts[] = {
        330005, 26229, 26210, 22209, 20690, 19027, 18418, 18367, 16283, 14990,
        16494, 15818, 16008, 11510, 10804, 10884, 12195, 11486, 10222, 9806, 
        8985, 5227, 5882
    };
    
    for(int i = 0; i < NUM_CHROMOSOMES - 1; i++) {
        genomeDataManager->setSNPCountPerChr(i + 1, snpCounts[i]);
    }
}

bool HaplotypePhasingProgram::execute() {
    clock_t startTime = clock();
    
    if(!validateInputFiles()) {
        return false;
    }
    
    initializeSNPCounts();
    
    if(configuration->isVerboseMode()) {
        printf("Loading genome data from: %s\n", configuration->getInputPath().c_str());
    }
    
    for(int chr = 1; chr < NUM_CHROMOSOMES; chr++) {
        int snpCount = genomeDataManager->getSNPCountPerChr(chr);
        size_t bufferSize = ((size_t)snpCount / 4 + ((snpCount % 4) > 0 ? 1 : 0)) * 
                           configuration->getNumberOfIndividuals();
        unsigned char* buffer = (unsigned char*)calloc(bufferSize, sizeof(unsigned char));
        
        if(buffer == nullptr) {
            printf("Error: Memory allocation failed for chromosome %d\n", chr);
            return false;
        }
        
        genomeDataManager->setGenomeBuffer(chr, buffer);
        
        if(!GenomeFileLoader::loadGenome(configuration->getInputPath().c_str(), chr,
                                        buffer, configuration->getNumberOfIndividuals(),
                                        genomeDataManager->getSNPCountPerChrArray(),
                                        genomeDataManager->getSNPCountInFileArray())) {
            free(buffer);
            genomeDataManager->setGenomeBuffer(chr, nullptr);
            return false;
        }
    }
    
    clock_t loadTime = clock();
    if(configuration->isVerboseMode()) {
        float elapsed = (float)(loadTime - startTime) / CLOCKS_PER_SEC;
        printf("Genome loading completed in %.2f seconds\n", elapsed);
    }
    
    if(!processIndividuals()) {
        return false;
    }
    
    clock_t endTime = clock();
    logExecutionStatistics(startTime, endTime);
    
    return true;
}

bool HaplotypePhasingProgram::processIndividuals() {
    int numberOfIndividuals = configuration->getNumberOfIndividuals();
    
    for(int individual = 0; individual < numberOfIndividuals; individual++) {
        if(configuration->isVerboseMode()) {
            printf("Processing individual %d of %d\n", individual + 1, numberOfIndividuals);
        }
        
        if(!phasingEngine->loadSegment(individual, individual, individual, individual,
                                      0, configuration->getAlgorithmVersion(), 0,
                                      configuration->getInputPath().c_str())) {
            printf("Warning: Failed to load segment for individual %d\n", individual);
            continue;
        }
        
        int relativeCount = relativeEngine->getRelativeCountAboveThreshold(0.1f);
        phasingEngine->setRelativeCount(relativeCount);
        phasingEngine->setBreakpointCount(1000);
        
        phasingEngine->executePhasingAlgorithm(individual, 1, 0, individual, 0.0f, 0,
                                              individual, individual, 0,
                                              configuration->getAlgorithmVersion(), 0,
                                              configuration->getOutputPath().c_str());
        
        for(int chr = 1; chr < NUM_CHROMOSOMES; chr++) {
            outputWriter->writeOutput(chr, configuration->getOutputPath().c_str());
        }
    }
    
    return true;
}

bool HaplotypePhasingProgram::validateInputFiles() const {
    FILE* testFile = fopen(configuration->getInputPath().c_str(), "r");
    if(testFile == nullptr) {
        printf("Error: Cannot access input path: %s\n", configuration->getInputPath().c_str());
        return false;
    }
    fclose(testFile);
    return true;
}

void HaplotypePhasingProgram::logExecutionStatistics(clock_t startTime, clock_t endTime) const {
    if(configuration->isVerboseMode()) {
        float totalTime = (float)(endTime - startTime) / CLOCKS_PER_SEC;
        printf("\n=== Execution Statistics ===\n");
        printf("Total execution time: %.2f seconds\n", totalTime);
        printf("Number of individuals processed: %d\n", configuration->getNumberOfIndividuals());
        printf("Output directory: %s\n", configuration->getOutputPath().c_str());
        printf("===========================\n");
    }
}

void HaplotypePhasingProgram::shutdown() {
    if(genomeDataManager) {
        genomeDataManager->reset();
    }
}

void HaplotypePhasingProgram::initializeChromosomeDividers() {
    // Implementation if needed
}


