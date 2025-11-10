/**
 * @file RelativeIdentificationEngine.h
 * @brief Advanced relative identification using PIHAT computation
 */

#ifndef RELATIVE_IDENTIFICATION_ENGINE_H
#define RELATIVE_IDENTIFICATION_ENGINE_H

#include "Interfaces.h"
#include "Constants.h"

namespace PhasingEngine {

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

}

#endif // RELATIVE_IDENTIFICATION_ENGINE_H


