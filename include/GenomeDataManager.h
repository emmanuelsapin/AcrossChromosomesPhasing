/**
 * @file GenomeDataManager.h
 * @brief Comprehensive genome data management with validation and caching
 */

#ifndef GENOME_DATA_MANAGER_H
#define GENOME_DATA_MANAGER_H

#include "Interfaces.h"
#include "Constants.h"
#include <cstddef>

namespace PhasingEngine {

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

}

#endif // GENOME_DATA_MANAGER_H


