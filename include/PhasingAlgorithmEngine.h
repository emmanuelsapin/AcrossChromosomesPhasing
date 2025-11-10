/**
 * @file PhasingAlgorithmEngine.h
 * @brief Core phasing algorithm implementation with multiple strategies
 */

#ifndef PHASING_ALGORITHM_ENGINE_H
#define PHASING_ALGORITHM_ENGINE_H

#include "Interfaces.h"
#include "ChromosomeDivider.h"
#include "Constants.h"
#include <memory>

namespace PhasingEngine {

class GenomeDataManager;
class RelativeIdentificationEngine;

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

}

#endif // PHASING_ALGORITHM_ENGINE_H


