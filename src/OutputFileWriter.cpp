/**
 * @file OutputFileWriter.cpp
 * @brief Implementation of OutputFileWriter
 */

#include "../include/OutputFileWriter.h"
#include "../include/GenomeDataManager.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace PhasingEngine;

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

bool OutputFileWriter::writeOutputBatch(const std::vector<int>& chromosomes, const char* outputPath) {
    bool success = true;
    for(int chr : chromosomes) {
        if(!writeOutput(chr, outputPath)) {
            success = false;
        }
    }
    return success;
}

bool OutputFileWriter::validateOutputPath(const char* path) const {
    FILE* testFile = fopen(path, "w");
    if(testFile == nullptr) {
        return false;
    }
    fclose(testFile);
    return true;
}

void OutputFileWriter::writePEDHeader(FILE* file, int chromosome) const {
    fprintf(file, "# Chromosome %d\n", chromosome);
}

void OutputFileWriter::writePEDRow(FILE* file, int individualID, int chromosome) const {
    fprintf(file, "%d %d 0 0 -9 0 ", individualID, individualID);
    for(int snp = 0; snp < genomeDataManager->getSNPCount(chromosome); snp++) {
        int genotype = genomeDataManager->getGenotype(chromosome, individualID, snp);
        fprintf(file, "%d %d ", genotype % 2, genotype / 2);
    }
    fprintf(file, "\n");
}


