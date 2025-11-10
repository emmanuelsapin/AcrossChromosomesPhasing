/**
 * @file RelativeIdentificationEngine.cpp
 * @brief Implementation of RelativeIdentificationEngine
 */

#include "../include/RelativeIdentificationEngine.h"
#include "../include/Constants.h"
#include <algorithm>
#include <vector>
#include <utility>

using namespace PhasingEngine;
using namespace PhasingEngine::Constants;

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

void RelativeIdentificationEngine::updateBestRelativeRankings() {
    // Implementation for updating rankings
}


