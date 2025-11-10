/**
 * @file Interfaces.h
 * @brief Abstract interfaces for phasing engine components
 */

#ifndef INTERFACES_H
#define INTERFACES_H

#include <string>

namespace PhasingEngine {

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

}

#endif // INTERFACES_H


