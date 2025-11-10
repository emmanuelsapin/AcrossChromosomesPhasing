/**
 * @file GenomeFileLoader.cpp
 * @brief Implementation of GenomeFileLoader
 */

#include "../include/GenomeFileLoader.h"
#include <cstdio>
#include <cstring>

using namespace PhasingEngine;

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

bool GenomeFileLoader::validateFileFormat(const char* filePath, FileFormat format) {
    FILE* testFile = fopen(filePath, "r");
    if(testFile == nullptr) {
        return false;
    }
    fclose(testFile);
    return true;
}

