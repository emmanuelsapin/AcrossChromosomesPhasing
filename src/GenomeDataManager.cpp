/**
 * @file GenomeDataManager.cpp
 * @brief Implementation of GenomeDataManager
 */

#include "../include/GenomeDataManager.h"
#include "../include/GenomeFileLoader.h"
#include "../include/Exceptions.h"
#include "../include/ErrorCodes.h"
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <string>

using namespace PhasingEngine;
using namespace PhasingEngine::Constants;
using namespace PhasingEngine::ErrorCodes;

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

bool GenomeDataManager::isValidChromosome(int chromosome) const {
    return chromosome >= 1 && chromosome < NUM_CHROMOSOMES;
}

bool GenomeDataManager::isValidIndividual(int individual) const {
    return individual >= 0 && individual < numberOfIndividuals;
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

bool GenomeDataManager::loadFromFile(const char* pathfile, int chromosome, int nbIndiv) {
    return GenomeFileLoader::loadGenome(pathfile, chromosome,
                                       genomes[chromosome], nbIndiv,
                                       snpCountPerChromosome, snpCountInFile);
}

