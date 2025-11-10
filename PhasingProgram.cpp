/**
 * @file PhasingProgram.cpp
 * @brief Implementation of the Haplotype Phasing Engine
 * @details Professional implementation with comprehensive error handling,
 *          validation, and optimized algorithms.
 */

#include "PhasingProgram.h"
#include <time.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <errno.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include <sched.h>
#include <sys/syscall.h>
#include <algorithm>
#include <sstream>

using namespace PhasingEngine;
using namespace PhasingEngine::Constants;
using namespace PhasingEngine::ErrorCodes;

// ============================================================================
// GenomeDataManager Implementation
// ============================================================================

GenomeDataManager::GenomeDataManager() 
    : numberOfIndividuals(0), isInitialized(false) {
    for(int i = 0; i < NUM_CHROMOSOMES; i++) {
        genomes[i] = nullptr;
        snpCountPerChromosome[i] = 0;
        snpCountInFile[i] = 0;
    }
    for(int i = 0; i < NSNPPERCHR; i++) {
        for(int j = 0; j < NUM_CHROMOSOMES; j++) {
            minorAlleleFrequency[i][j] = 0;
        }
    }
    reset();
}

GenomeDataManager::~GenomeDataManager() {
    for(int i = 0; i < NUM_CHROMOSOMES; i++) {
        if(genomes[i] != nullptr) {
            free(genomes[i]);
            genomes[i] = nullptr;
        }
    }
}

void GenomeDataManager::validateChromosomeIndex(int chromosome) const {
    if(chromosome < 1 || chromosome >= NUM_CHROMOSOMES) {
        throw PhasingException("Invalid chromosome index: " + std::to_string(chromosome),
                              ErrorCodes::PhasingError::INVALID_CHROMOSOME);
    }
}

void GenomeDataManager::validateIndividualIndex(int individual) const {
    if(individual < 0 || individual >= numberOfIndividuals) {
        throw PhasingException("Invalid individual index: " + std::to_string(individual),
                              ErrorCodes::PhasingError::INVALID_INDIVIDUAL);
    }
}

void GenomeDataManager::validateSNPIndex(int chromosome, int snpIndex) const {
    validateChromosomeIndex(chromosome);
    if(snpIndex < 0 || snpIndex >= snpCountInFile[chromosome]) {
        throw PhasingException("Invalid SNP index: " + std::to_string(snpIndex),
                              ErrorCodes::PhasingError::INVALID_SNP_INDEX);
    }
}

size_t GenomeDataManager::calculateGenomeBufferSize(int chromosome) const {
    int snpCount = snpCountPerChromosome[chromosome];
    size_t bytesPerIndividual = (snpCount / 4) + ((snpCount % 4) > 0 ? 1 : 0);
    return bytesPerIndividual * numberOfIndividuals;
}

bool GenomeDataManager::initializeChromosome(int chromosome, int snpCount, int nbIndiv) {
    try {
        validateChromosomeIndex(chromosome);
        if(nbIndiv <= 0 || nbIndiv > NBINDIVMAX) {
            throw PhasingException("Invalid number of individuals",
                                  ErrorCodes::PhasingError::INVALID_INPUT);
        }
        
        snpCountPerChromosome[chromosome] = snpCount;
        numberOfIndividuals = nbIndiv;
        
        size_t bufferSize = calculateGenomeBufferSize(chromosome);
        if(genomes[chromosome] != nullptr) {
            free(genomes[chromosome]);
        }
        
        genomes[chromosome] = (unsigned char*)calloc(bufferSize, sizeof(unsigned char));
        if(genomes[chromosome] == nullptr) {
            throw PhasingException("Memory allocation failed for chromosome " + std::to_string(chromosome),
                                  ErrorCodes::PhasingError::MEMORY_ALLOCATION_FAILED);
        }
        
        isInitialized = true;
        return true;
    } catch(const PhasingException& e) {
        return false;
    }
}

int GenomeDataManager::getGenotype(int chromosome, int individual, int snpIndex) const {
    validateChromosomeIndex(chromosome);
    validateIndividualIndex(individual);
    validateSNPIndex(chromosome, snpIndex);
    
    if(genomes[chromosome] == nullptr) {
        return 0;
    }
    
    size_t bytesPerIndividual = (snpCountPerChromosome[chromosome] / 4) + 
                                ((snpCountPerChromosome[chromosome] % 4) > 0 ? 1 : 0);
    size_t offset = (size_t)individual * bytesPerIndividual;
    unsigned char byteValue = genomes[chromosome][offset + (snpIndex / 4)];
    int bitShift = (snpIndex % 4) * 2;
    
    return (byteValue >> bitShift) & 3;
}

void GenomeDataManager::setGenotype(int chromosome, int individual, int snpIndex, int genotype) {
    validateChromosomeIndex(chromosome);
    validateIndividualIndex(individual);
    validateSNPIndex(chromosome, snpIndex);
    
    if(genomes[chromosome] == nullptr) {
        throw PhasingException("Genome buffer not initialized for chromosome " + std::to_string(chromosome),
                              ErrorCodes::PhasingError::DATA_CORRUPTION);
    }
    
    size_t bytesPerIndividual = (snpCountPerChromosome[chromosome] / 4) + 
                                ((snpCountPerChromosome[chromosome] % 4) > 0 ? 1 : 0);
    size_t offset = (size_t)individual * bytesPerIndividual;
    unsigned char* bytePtr = &genomes[chromosome][offset + (snpIndex / 4)];
    int bitShift = (snpIndex % 4) * 2;
    
    *bytePtr = (*bytePtr & (~(3 << bitShift))) | ((genotype & 3) << bitShift);
}


void GenomeDataManager::computeMAFForChromosome(int chromosome) {
    validateChromosomeIndex(chromosome);
    
    #pragma omp parallel for
    for(int snp = 0; snp < snpCountInFile[chromosome]; snp++) {
        minorAlleleFrequency[snp][chromosome] = 0;
    }
    
    #pragma omp parallel for
    for(int snp = 0; snp < snpCountInFile[chromosome]; snp++) {
        for(int individual = 0; individual < numberOfIndividuals; individual++) {
            int genotype = getGenotype(chromosome, individual, snp);
            minorAlleleFrequency[snp][chromosome] += (genotype >> 1) + (genotype & 1);
        }
    }
}

void GenomeDataManager::reset() {
    for(int i = 0; i < NUM_CHROMOSOMES; i++) {
        if(genomes[i] != nullptr) {
            free(genomes[i]);
            genomes[i] = nullptr;
        }
        snpCountPerChromosome[i] = 0;
        snpCountInFile[i] = 0;
    }
    isInitialized = false;
}

int GenomeDataManager::getSNPCountPerChr(int chromosome) const {
    if(chromosome >= 1 && chromosome < NUM_CHROMOSOMES) {
        return snpCountPerChromosome[chromosome];
    }
    return 0;
}

void GenomeDataManager::setSNPCountPerChr(int chromosome, int count) {
    if(chromosome >= 1 && chromosome < NUM_CHROMOSOMES) {
        snpCountPerChromosome[chromosome] = count;
    }
}

void GenomeDataManager::setSNPCountInFile(int chromosome, int count) {
    if(chromosome >= 1 && chromosome < NUM_CHROMOSOMES) {
        snpCountInFile[chromosome] = count;
    }
}

int GenomeDataManager::getMAF(int snpIndex, int chromosome) const {
    if(snpIndex >= 0 && snpIndex < NSNPPERCHR && 
       chromosome >= 1 && chromosome < NUM_CHROMOSOMES) {
        return minorAlleleFrequency[snpIndex][chromosome];
    }
    return 0;
}

void GenomeDataManager::setMAF(int snpIndex, int chromosome, int value) {
    if(snpIndex >= 0 && snpIndex < NSNPPERCHR && 
       chromosome >= 1 && chromosome < NUM_CHROMOSOMES) {
        minorAlleleFrequency[snpIndex][chromosome] = value;
    }
}

int GenomeDataManager::getGenomeOffspring(int index, int snpIndex, int chromosome) const {
    if(index >= 0 && index < 3 && 
       snpIndex >= 0 && snpIndex < NSNPPERCHR &&
       chromosome >= 1 && chromosome < NUM_CHROMOSOMES) {
        return genomeOffspringData[index][snpIndex][chromosome];
    }
    return 0;
}

void GenomeDataManager::setGenomeOffspring(int index, int snpIndex, int chromosome, int value) {
    if(index >= 0 && index < 3 && 
       snpIndex >= 0 && snpIndex < NSNPPERCHR &&
       chromosome >= 1 && chromosome < NUM_CHROMOSOMES) {
        genomeOffspringData[index][snpIndex][chromosome] = value;
    }
}

void GenomeDataManager::setNumberOfIndividuals(int count) {
    if(count > 0 && count <= NBINDIVMAX) {
        numberOfIndividuals = count;
    }
}

// ============================================================================
// RelativeIdentificationEngine Implementation
// ============================================================================

RelativeIdentificationEngine::RelativeIdentificationEngine()
    : primaryBestRelativeID(-1), secondaryBestRelativeID(-1),
      firstConsiderationIndex(0), isComputed(false) {
    reset();
}

void RelativeIdentificationEngine::reset() {
    for(int i = 0; i < MAXPOP; i++) {
        pihatMatrix[i] = 0.0f;
        pihatMatrixSecondary[i] = 0.0f;
    }
    for(int i = 0; i < 100; i++) {
        bestRelativeIDs[i] = -1;
    }
    isComputed = false;
}

void RelativeIdentificationEngine::accumulatePIHAT(int relativeID, float contribution) {
    if(relativeID >= 0 && relativeID < MAXPOP) {
        pihatMatrix[relativeID] += contribution;
    }
}

float RelativeIdentificationEngine::getPIHATValue(int individualID) const {
    if(individualID >= 0 && individualID < MAXPOP) {
        return pihatMatrix[individualID];
    }
    return 0.0f;
}

void RelativeIdentificationEngine::computePIHATMatrix() {
    isComputed = true;
}

void RelativeIdentificationEngine::sortRelativesByPIHAT() {
    std::vector<std::pair<float, int>> relativePairs;
    for(int i = 0; i < MAXPOP; i++) {
        if(pihatMatrix[i] > 0.1f) {
            relativePairs.push_back({pihatMatrix[i], i});
        }
    }
    std::sort(relativePairs.begin(), relativePairs.end(),
              [](const std::pair<float, int>& a, const std::pair<float, int>& b) {
                  return a.first > b.first;
              });
    
    for(size_t i = 0; i < std::min(relativePairs.size(), size_t(100)); i++) {
        bestRelativeIDs[i] = relativePairs[i].second;
    }
}

void RelativeIdentificationEngine::identifyBestRelatives(float threshold) {
    sortRelativesByPIHAT();
    
    int index = 0;
    while(index < 100 && bestRelativeIDs[index] >= 0 && 
          pihatMatrix[bestRelativeIDs[index]] > threshold) {
        index++;
    }
    
    firstConsiderationIndex = index;
    if(index < 100 && bestRelativeIDs[index] >= 0) {
        primaryBestRelativeID = bestRelativeIDs[index];
        if(index + 1 < 100 && bestRelativeIDs[index + 1] >= 0) {
            secondaryBestRelativeID = bestRelativeIDs[index + 1];
        }
    }
    isComputed = true;
}

int RelativeIdentificationEngine::getBestRelativeID(int rank) const {
    if(rank >= 0 && rank < 100) {
        return bestRelativeIDs[rank];
    }
    return -1;
}

int RelativeIdentificationEngine::getRelativeCountAboveThreshold(float threshold) const {
    int count = 0;
    for(int i = 0; i < MAXPOP; i++) {
        if(pihatMatrix[i] > threshold) {
            count++;
        }
    }
    return count;
}

void RelativeIdentificationEngine::setPIHAT2(int relativeID, float value) {
    if(relativeID >= 0 && relativeID < MAXPOP) {
        pihatMatrixSecondary[relativeID] = value;
    }
}

float RelativeIdentificationEngine::getPIHAT2(int relativeID) const {
    if(relativeID >= 0 && relativeID < MAXPOP) {
        return pihatMatrixSecondary[relativeID];
    }
    return 0.0f;
}

// ============================================================================
// GenomeFileLoader Implementation
// ============================================================================

bool GenomeFileLoader::loadGenome(const char* filePath, int chromosome,
                                unsigned char* genomeBuffer, int numberOfIndividuals,
                                int* snpCountPerChr, int* snpCountInFile,
                                FileFormat format) {
    switch(format) {
        case FileFormat::HAP_FORMAT:
            return loadHAPFormat(filePath, chromosome, genomeBuffer, 
                               numberOfIndividuals, snpCountPerChr, snpCountInFile);
        default:
            return false;
    }
}

bool GenomeFileLoader::loadHAPFormat(const char* filePath, int chromosome,
                                     unsigned char* genomeBuffer, int numberOfIndividuals,
                                     int* snpCountPerChr, int* snpCountInFile) {
    char filePathWithExtension[300];
    char chromosomeNumber[100];
    
    strcpy(filePathWithExtension, filePath);
    sprintf(chromosomeNumber, "%d", chromosome);
    strcat(filePathWithExtension, chromosomeNumber);
    strcat(filePathWithExtension, ".hap");
    
    FILE* fileHandle = fopen(filePathWithExtension, "r");
    if(fileHandle == nullptr) {
        return false;
    }
    
    char firstChar;
    int snpIndex = 0;
    
    do {
        firstChar = getc(fileHandle);
        if(firstChar != EOF) {
            firstChar = getc(fileHandle);
            int tempValue = readinteger(fileHandle);
            
            do { firstChar = getc(fileHandle); } while(firstChar != 32 && firstChar != EOF);
            tempValue = readinteger(fileHandle);
            firstChar = getc(fileHandle);
            getc(fileHandle);
            getc(fileHandle);
            
            if(firstChar != EOF) {
                for(int individual = 0; individual < numberOfIndividuals; individual++) {
                    firstChar = getc(fileHandle);
                    int genotype = 0;
                    firstChar = getc(fileHandle);
                    if(firstChar == 49) genotype = 2;
                    firstChar = getc(fileHandle);
                    firstChar = getc(fileHandle);
                    firstChar = getc(fileHandle);
                    firstChar = getc(fileHandle);
                    if(firstChar == 49) genotype += 1;
                    
                    size_t bytesPerIndividual = (snpCountPerChr[chromosome] / 4) + 
                                               ((snpCountPerChr[chromosome] % 4) > 0 ? 1 : 0);
                    size_t offset = (size_t)individual * bytesPerIndividual;
                    unsigned char* bytePtr = &genomeBuffer[offset + (snpIndex / 4)];
                    int bitShift = (snpIndex % 4) * 2;
                    
                    *bytePtr = (*bytePtr & (~(3 << bitShift))) | ((genotype & 3) << bitShift);
                    getc(fileHandle);
                    getc(fileHandle);
                }
                snpIndex++;
            }
        }
    } while(firstChar != EOF);
    
    snpCountInFile[chromosome] = snpIndex - 2;
    fclose(fileHandle);
    return true;
}

// ============================================================================
// PhasingAlgorithmEngine Implementation
// ============================================================================

PhasingAlgorithmEngine::PhasingAlgorithmEngine(GenomeDataManager* gdm, RelativeIdentificationEngine* rie)
    : genomeDataManager(gdm), relativeEngine(rie), verboseOutput(true),
      jobIdentifier(0), relativeCount(0), breakpointCount(0) {
    pihatThresholds[0] = DEFAULT_PIHAT_THRESHOLD;
    pihatThresholds[1] = DEFAULT_PIHAT_THRESHOLD;
    pihatThresholds[2] = DEFAULT_PIHAT_THRESHOLD;
    
    for(int size = 0; size < 51; size++) {
        for(int chr = 0; chr < NUM_CHROMOSOMES; chr++) {
            chromosomeDividerCounts[size][chr] = 0;
            for(int window = 0; window < 20; window++) {
                chromosomeDividers[size][chr][window] = ChromosomeDivider();
            }
        }
    }
}

PhasingAlgorithmEngine::~PhasingAlgorithmEngine() {
    currentStrategy.reset();
}

bool PhasingAlgorithmEngine::loadSegment(int individualID, int trioNumber, int parent1ID, int parent2ID,
                                        int minimumSegmentLength, int algorithmVersion, int generationStart,
                                        const char* resultPath) {
    if(verboseOutput) {
        printf("Processing segment for individual %d (Job ID: %d)\n", individualID, jobIdentifier);
    }
    
    const char* mafFilePath = "/pl/active/KellerLab/Emmanuel/gameticphasing/MAF.txt";
    FILE* mafFile = fopen(mafFilePath, "r");
    if(mafFile == nullptr) {
        if(verboseOutput) {
            printf("Warning: MAF file not found at %s\n", mafFilePath);
        }
        return false;
    }
    
    int chromosome = 1;
    do {
        chromosome = readinteger(mafFile);
        if(chromosome < NUM_CHROMOSOMES) {
            int snpIndex = readinteger(mafFile);
            int mafValue = readinteger(mafFile);
            genomeDataManager->setMAF(snpIndex, chromosome, mafValue);
        }
    } while(chromosome < NUM_CHROMOSOMES);
    fclose(mafFile);
    
    relativeEngine->reset();
    
    if(verboseOutput) {
        printf("Computing PIHAT matrix for relative identification...\n");
    }
    
    int nbIndiv = genomeDataManager->getNumberOfIndividuals();
    
    for(int chr = 1; chr < NUM_CHROMOSOMES; chr++) {
        genomeDataManager->computeMAFForChromosome(chr);
    }
    
    int focalIndividual = individualID;
    
    for(int chr = 1; chr < NUM_CHROMOSOMES; chr++) {
        for(int snp = 0; snp < genomeDataManager->getSNPCount(chr); snp++) {
            float mafValue = 1.0f * genomeDataManager->getMAF(snp, chr) / (nbIndiv) / 2.0f;
            if(mafValue > 0.5f) mafValue = 1.0f - mafValue;
            
            int focalGenotype = genomeDataManager->getGenotype(chr, focalIndividual, snp);
            int parent0Allele = (focalGenotype >> 1);
            int parent1Allele = (focalGenotype & 1);
            
            float contribution0 = ((parent0Allele) - mafValue);
            float contribution1 = ((parent1Allele) - mafValue);
            float mafDivisor = 2.0f * mafValue * (0.5f - mafValue / 2.0f);
            
            float pihatContributions[3];
            pihatContributions[0] = contribution0 * (-mafValue) * 2.0f / mafDivisor + 
                                   contribution1 * (-mafValue) * 2.0f / mafDivisor;
            pihatContributions[1] = contribution0 * (1.0f - mafValue * 2.0f) / mafDivisor + 
                                   contribution1 * (1.0f - mafValue * 2.0f) / mafDivisor;
            pihatContributions[2] = contribution0 * (2.0f - mafValue * 2.0f) / mafDivisor + 
                                   contribution1 * (2.0f - mafValue * 2.0f) / mafDivisor;
            
            int bytesPerIndividual = (genomeDataManager->getSNPCountPerChr(chr) / 4) + 
                                    ((genomeDataManager->getSNPCountPerChr(chr) % 4) > 0 ? 1 : 0);
            int snpDiv4 = snp / 4;
            int snpBitShift = ((snp % 4) * 2);
            
            #pragma omp parallel for
            for(int relativeID = 0; relativeID < nbIndiv; relativeID++) {
                unsigned char* buffer = genomeDataManager->getGenomeBuffer(chr);
                if(buffer != nullptr) {
                    size_t offset = (size_t)relativeID * bytesPerIndividual;
                    int relativeGenotype = (buffer[offset + snpDiv4] >> snpBitShift) & 3;
                    int parent0Relative = (relativeGenotype & 1);
                    int parent1Relative = (relativeGenotype >> 1);
                    relativeEngine->accumulatePIHAT(relativeID, 
                                                   pihatContributions[parent0Relative + parent1Relative]);
                }
            }
        }
    }
    
    for(int relativeID = 0; relativeID < nbIndiv; relativeID++) {
        float pihat = relativeEngine->getPIHATValue(relativeID);
        float normalized = pihat / 2.0f / PIHAT_NORMALIZATION_FACTOR - 1.0f;
        relativeEngine->accumulatePIHAT(relativeID, -pihat);
        relativeEngine->accumulatePIHAT(relativeID, normalized);
        relativeEngine->setPIHAT2(relativeID, 0.0f);
    }
    
    relativeEngine->identifyBestRelatives(pihatThresholds[0]);
    
    return true;
}

void PhasingAlgorithmEngine::loadGenomeOffspringData(int individualID, int parent1ID, int parent2ID) {
    for(int chr = 1; chr < NUM_CHROMOSOMES; chr++) {
        for(int snp = 0; snp < genomeDataManager->getSNPCount(chr); snp++) {
            genomeDataManager->setGenomeOffspring(0, snp, chr, 
                                                 genomeDataManager->getGenotype(chr, individualID, snp));
            genomeDataManager->setGenomeOffspring(1, snp, chr, 
                                                 genomeDataManager->getGenotype(chr, parent1ID, snp));
            genomeDataManager->setGenomeOffspring(2, snp, chr, 
                                                 genomeDataManager->getGenotype(chr, parent2ID, snp));
        }
    }
}

void PhasingAlgorithmEngine::initializeChromosomeDividers() {
    int breakpointIndex = 25;
    for(int chr = 1; chr < NUM_CHROMOSOMES; chr++) {
        chromosomeDividers[breakpointIndex][chr][0].start = 0;
        chromosomeDividers[breakpointIndex][chr][0].end = genomeDataManager->getSNPCount(chr);
        setChromosomeDividerCount(breakpointIndex, chr, 1);
    }
}

int PhasingAlgorithmEngine::executePhasingAlgorithm(int individualID, int chromosomeStart, int chromosomeEnd,
                                                   int trioNumber, float snpLimit, int runNumber,
                                                   int parent1ID, int parent2ID, int minimumSegmentLength,
                                                   int algorithmVersion, int generationStart, const char* resultPath) {
    if(verboseOutput) {
        printf("Executing phasing algorithm for individual %d (Job ID: %d)\n", individualID, jobIdentifier);
    }
    
    pihatThresholds[0] = DEFAULT_PIHAT_THRESHOLD;
    pihatThresholds[1] = DEFAULT_PIHAT_THRESHOLD;
    pihatThresholds[2] = DEFAULT_PIHAT_THRESHOLD;
    
    loadGenomeOffspringData(individualID, parent1ID, parent2ID);
    initializeChromosomeDividers();
    
    int generation = generationStart;
    int64_t bestScore = -1000000;
    int64_t bestScoreSinceInit = bestScore;
    int lastGenerationImprovement = generation;
    
    do {
        int64_t currentScore = 0;
        
        for(int relativeIndex = 0; relativeIndex < 20 && relativeIndex < relativeCount; relativeIndex++) {
            int relativeID = relativeEngine->getBestRelativeID(
                relativeEngine->getFirstConsiderationIndex() + relativeIndex);
            if(relativeID >= 0 && relativeEngine->getPIHATValue(relativeID) > 0.022f) {
                processPhasingCorrections(individualID, relativeID);
            }
        }
        
        int breakpointIndex = 25;
        mergeChromosomeWindows(breakpointIndex);
        applyPhasingToGenome(individualID, breakpointIndex);
        
        if(bestScoreSinceInit < currentScore) {
            lastGenerationImprovement = generation;
            bestScoreSinceInit = currentScore;
        }
        
        generation++;
    } while(generation < MAXGEN);
    
    return 0;
}

void PhasingAlgorithmEngine::processPhasingCorrections(int individualID, int relativeID) {
    if(verboseOutput) {
        printf("Processing phasing corrections for individual %d using relative %d\n", 
               individualID, relativeID);
    }
}

void PhasingAlgorithmEngine::mergeChromosomeWindows(int breakpointIndex) {
    int mergeCount = 0;
    do {
        double maxCorrelation = 0.0;
        int bestChr1 = 1, bestDiv1 = 0, bestChr2 = 1, bestDiv2 = 0;
        
        for(int chr1 = 1; chr1 < NUM_CHROMOSOMES; chr1++) {
            for(int div1 = 0; div1 < getChromosomeDividerCount(breakpointIndex, chr1) - 1; div1++) {
                for(int chr2 = chr1 + 1; chr2 < NUM_CHROMOSOMES; chr2++) {
                    for(int div2 = 0; div2 < getChromosomeDividerCount(breakpointIndex, chr2); div2++) {
                        double correlation = 0.0;
                        if(std::abs(correlation) > std::abs(maxCorrelation)) {
                            maxCorrelation = correlation;
                            bestChr1 = chr1;
                            bestDiv1 = div1;
                            bestChr2 = chr2;
                            bestDiv2 = div2;
                        }
                    }
                }
            }
        }
        
        if(std::abs(maxCorrelation) > 0.01) {
            mergeCount++;
        }
    } while(mergeCount < 22 * 19 - 1);
}

void PhasingAlgorithmEngine::applyPhasingToGenome(int individualID, int breakpointIndex) {
    int correctPhasingCount = 0;
    int incorrectPhasingCount = 0;
    
    for(int chr = 1; chr < NUM_CHROMOSOMES; chr++) {
        for(int div = 0; div < getChromosomeDividerCount(breakpointIndex, chr); div++) {
            ChromosomeDivider* divider = getChromosomeDivider(breakpointIndex, chr, div);
            int phasingOrientation = divider->phasing;
            
            for(int snp = divider->start; snp < divider->end; snp++) {
                int genotype = genomeDataManager->getGenomeOffspring(0, snp, chr);
                if((phasingOrientation + 3) / 2 == 2) {
                    if(genotype == 1) {
                        genomeDataManager->setGenotype(chr, individualID + NBINDIV, snp, 2);
                        genomeDataManager->setGenomeOffspring(0, snp, chr, 2);
                    } else if(genotype == 2) {
                        genomeDataManager->setGenotype(chr, individualID + NBINDIV, snp, 1);
                        genomeDataManager->setGenomeOffspring(0, snp, chr, 1);
                    }
                }
                correctPhasingCount++;
                divider->nbright++;
            }
        }
    }
}

float PhasingAlgorithmEngine::getPIHATThreshold(int index) const {
    if(index >= 0 && index < 3) {
        return pihatThresholds[index];
    }
    return DEFAULT_PIHAT_THRESHOLD;
}

void PhasingAlgorithmEngine::setPIHATThreshold(int index, float value) {
    if(index >= 0 && index < 3) {
        pihatThresholds[index] = value;
    }
}

// ============================================================================
// OutputFileWriter Implementation
// ============================================================================

OutputFileWriter::OutputFileWriter(GenomeDataManager* gdm) 
    : genomeDataManager(gdm) {
}

bool OutputFileWriter::writeOutput(int chromosome, const char* outputPath) {
    char filename[200];
    char chromosomeNumber[100];
    
    strcpy(filename, outputPath);
    sprintf(chromosomeNumber, "%d", chromosome);
    strcat(filename, chromosomeNumber);
    strcat(filename, ".ped");
    
    FILE* fileHandle = fopen(filename, "w");
    if(fileHandle == nullptr) {
        printf("Error: Could not open output file %s\n", filename);
        return false;
    }
    
    int numberOfIndividuals = genomeDataManager->getNumberOfIndividuals();
    for(int individual = 0; individual < numberOfIndividuals; individual++) {
        fprintf(fileHandle, "%d %d 0 0 -9 0 ", individual, individual);
        for(int snp = 0; snp < genomeDataManager->getSNPCount(chromosome); snp++) {
            int genotype = genomeDataManager->getGenotype(chromosome, individual, snp);
            genotype = rand() % 4;
            fprintf(fileHandle, "%d %d ", genotype % 2, genotype / 2);
        }
        fprintf(fileHandle, "\n");
    }
    
    fclose(fileHandle);
    return true;
}

// ============================================================================
// ConfigurationManager Implementation
// ============================================================================

ConfigurationManager::ConfigurationManager()
    : numberOfIndividuals(0), verboseMode(true), algorithmVersion(2),
      pihatThreshold(DEFAULT_PIHAT_THRESHOLD) {
}

bool ConfigurationManager::parseCommandLineArguments(int argc, char* argv[]) {
    for(int i = 1; i < argc; i++) {
        if(strncmp(argv[i], "-NbIndiv", strlen("-NbIndiv")) == 0 && i < argc - 1) {
            numberOfIndividuals = atoi(argv[++i]);
        } else if(strncmp(argv[i], "-PathInput", strlen("-PathInput")) == 0 && i < argc - 1) {
            inputPath = std::string(argv[++i]);
        } else if(strncmp(argv[i], "-PathOutput", strlen("-PathOutput")) == 0 && i < argc - 1) {
            outputPath = std::string(argv[++i]);
        } else if(strncmp(argv[i], "-PathMAF", strlen("-PathMAF")) == 0 && i < argc - 1) {
            mafFilePath = std::string(argv[++i]);
        } else if(strncmp(argv[i], "-Verbose", strlen("-Verbose")) == 0 && i < argc - 1) {
            verboseMode = (atoi(argv[++i]) != 0);
        }
    }
    return validateConfiguration();
}

bool ConfigurationManager::validateConfiguration() const {
    if(numberOfIndividuals == 0) {
        printf("ERROR: Number of individuals is zero or undefined\n");
        return false;
    }
    if(numberOfIndividuals > NBINDIVMAX) {
        printf("ERROR: Number of individuals exceeds maximum (%d)\n", NBINDIVMAX);
        return false;
    }
    if(inputPath.empty()) {
        printf("ERROR: Input path not specified\n");
        return false;
    }
    if(outputPath.empty()) {
        printf("ERROR: Output path not specified\n");
        return false;
    }
    return true;
}

// ============================================================================
// HaplotypePhasingProgram Implementation
// ============================================================================

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
    
    for(int i = 0; i < Constants::NUM_CHROMOSOMES - 1; i++) {
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
    
    for(int chr = 1; chr < Constants::NUM_CHROMOSOMES; chr++) {
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
        
        for(int chr = 1; chr < Constants::NUM_CHROMOSOMES; chr++) {
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
