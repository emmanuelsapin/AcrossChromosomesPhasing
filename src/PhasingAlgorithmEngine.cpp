/**
 * @file PhasingAlgorithmEngine.cpp
 * @brief Implementation of PhasingAlgorithmEngine
 */

#include "../include/PhasingAlgorithmEngine.h"
#include "../include/GenomeDataManager.h"
#include "../include/RelativeIdentificationEngine.h"
#include "../include/ChromosomeDivider.h"
#include "../include/Constants.h"
#include "../include/utils/readinteger.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

using namespace PhasingEngine;
using namespace PhasingEngine::Constants;

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

ChromosomeDivider* PhasingAlgorithmEngine::getChromosomeDivider(int sizeIndex, int chromosome, int window) {
    if(sizeIndex >= 0 && sizeIndex < 50 && 
       chromosome >= 0 && chromosome < NUM_CHROMOSOMES &&
       window >= 0 && window < 20) {
        return &chromosomeDividers[sizeIndex][chromosome][window];
    }
    return nullptr;
}

int PhasingAlgorithmEngine::getChromosomeDividerCount(int sizeIndex, int chromosome) const {
    if(sizeIndex >= 0 && sizeIndex < 51 && 
       chromosome >= 0 && chromosome < NUM_CHROMOSOMES) {
        return chromosomeDividerCounts[sizeIndex][chromosome];
    }
    return 0;
}

void PhasingAlgorithmEngine::setChromosomeDividerCount(int sizeIndex, int chromosome, int count) {
    if(sizeIndex >= 0 && sizeIndex < 51 && 
       chromosome >= 0 && chromosome < NUM_CHROMOSOMES) {
        chromosomeDividerCounts[sizeIndex][chromosome] = count;
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

void PhasingAlgorithmEngine::setPhasingStrategy(std::unique_ptr<IPhasingStrategy> strategy) {
    currentStrategy = std::move(strategy);
}

void PhasingAlgorithmEngine::computeCorrelationMatrix(int chromosome1, int chromosome2, 
                                                     double* correlationMatrix) const {
    // Implementation for correlation computation
}

void PhasingAlgorithmEngine::optimizeWindowMerging(int breakpointIndex) {
    // Implementation for window merging optimization
}

