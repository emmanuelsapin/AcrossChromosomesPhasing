/**
 * @file ConfigurationManager.cpp
 * @brief Implementation of ConfigurationManager
 */

#include "../include/ConfigurationManager.h"
#include "../include/Constants.h"
#include <cstring>
#include <cstdio>

using namespace PhasingEngine;
using namespace PhasingEngine::Constants;

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


