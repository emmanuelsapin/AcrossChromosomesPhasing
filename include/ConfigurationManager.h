/**
 * @file ConfigurationManager.h
 * @brief Centralized configuration management
 */

#ifndef CONFIGURATION_MANAGER_H
#define CONFIGURATION_MANAGER_H

#include <string>

namespace PhasingEngine {

/**
 * @class ConfigurationManager
 * @brief Centralized configuration management
 */
class ConfigurationManager {
private:
    int numberOfIndividuals;
    std::string inputPath;
    std::string outputPath;
    std::string mafFilePath;
    bool verboseMode;
    int algorithmVersion;
    float pihatThreshold;
    
public:
    ConfigurationManager();
    
    bool parseCommandLineArguments(int argc, char* argv[]);
    bool validateConfiguration() const;
    
    int getNumberOfIndividuals() const { return numberOfIndividuals; }
    void setNumberOfIndividuals(int count) { numberOfIndividuals = count; }
    
    std::string getInputPath() const { return inputPath; }
    void setInputPath(const std::string& path) { inputPath = path; }
    
    std::string getOutputPath() const { return outputPath; }
    void setOutputPath(const std::string& path) { outputPath = path; }
    
    std::string getMAFFilePath() const { return mafFilePath; }
    void setMAFFilePath(const std::string& path) { mafFilePath = path; }
    
    bool isVerboseMode() const { return verboseMode; }
    void setVerboseMode(bool enabled) { verboseMode = enabled; }
    
    int getAlgorithmVersion() const { return algorithmVersion; }
    void setAlgorithmVersion(int version) { algorithmVersion = version; }
    
    float getPIHATThreshold() const { return pihatThreshold; }
    void setPIHATThreshold(float threshold) { pihatThreshold = threshold; }
};

}

#endif // CONFIGURATION_MANAGER_H


