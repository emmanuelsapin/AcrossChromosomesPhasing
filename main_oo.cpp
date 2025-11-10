/**
 * @file main_oo.cpp
 * @brief Main entry point for the Haplotype Phasing Program
 * @details Professional implementation with comprehensive error handling
 *          and execution statistics.
 */

#include "include/PhasingProgram.h"
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <exception>

using namespace PhasingEngine;

/**
 * @brief Main application entry point
 * @param argc Number of command line arguments
 * @param argv Command line arguments array
 * @return Exit code (0 on success, non-zero on error)
 */
int main(int argc, char* argv[]) {
    clock_t programStartTime = clock();
    
    try {
        HaplotypePhasingProgram phasingApplication;
        
        std::cout << "===========================================" << std::endl;
        std::cout << "  Haplotype Phasing Engine v2.0" << std::endl;
        std::cout << "  Professional Implementation" << std::endl;
        std::cout << "===========================================" << std::endl;
        std::cout << std::endl;
        
        if(!phasingApplication.initialize(argc, argv)) {
            std::cerr << "ERROR: Failed to initialize phasing application" << std::endl;
            std::cerr << "Please check your command line arguments:" << std::endl;
            std::cerr << "  -NbIndiv <number>     : Number of individuals" << std::endl;
            std::cerr << "  -PathInput <path>     : Input file path prefix" << std::endl;
            std::cerr << "  -PathOutput <path>    : Output file path prefix" << std::endl;
            std::cerr << "  -PathMAF <path>       : MAF file path (optional)" << std::endl;
            std::cerr << "  -Verbose <0|1>        : Enable verbose output" << std::endl;
            return 1;
        }
        
        ConfigurationManager* config = phasingApplication.getConfiguration();
        std::cout << "Configuration loaded successfully:" << std::endl;
        std::cout << "  Individuals: " << config->getNumberOfIndividuals() << std::endl;
        std::cout << "  Input path: " << config->getInputPath() << std::endl;
        std::cout << "  Output path: " << config->getOutputPath() << std::endl;
        std::cout << "  Verbose mode: " << (config->isVerboseMode() ? "Enabled" : "Disabled") << std::endl;
        std::cout << std::endl;
        
        clock_t executionStartTime = clock();
        float initializationTime = (float)(executionStartTime - programStartTime) / CLOCKS_PER_SEC;
        std::cout << "Initialization completed in " << initializationTime << " seconds" << std::endl;
        std::cout << std::endl;
        
        std::cout << "Starting phasing pipeline execution..." << std::endl;
        std::cout << "-------------------------------------------" << std::endl;
        
        if(!phasingApplication.execute()) {
            std::cerr << "ERROR: Phasing execution failed" << std::endl;
            return 1;
        }
        
        clock_t programEndTime = clock();
        float totalExecutionTime = (float)(programEndTime - programStartTime) / CLOCKS_PER_SEC;
        float phasingTime = (float)(programEndTime - executionStartTime) / CLOCKS_PER_SEC;
        
        std::cout << std::endl;
        std::cout << "===========================================" << std::endl;
        std::cout << "Execution Summary:" << std::endl;
        std::cout << "  Total time: " << totalExecutionTime << " seconds" << std::endl;
        std::cout << "  Phasing time: " << phasingTime << " seconds" << std::endl;
        std::cout << "  Status: SUCCESS" << std::endl;
        std::cout << "===========================================" << std::endl;
        
        phasingApplication.shutdown();
        
        return 0;
        
    } catch(const PhasingException& e) {
        std::cerr << "Phasing Exception: " << e.what() << std::endl;
        std::cerr << "Error Code: " << static_cast<int>(e.getErrorCode()) << std::endl;
        return static_cast<int>(e.getErrorCode());
        
    } catch(const std::exception& e) {
        std::cerr << "Standard Exception: " << e.what() << std::endl;
        return 1;
        
    } catch(...) {
        std::cerr << "Unknown exception occurred" << std::endl;
        return 1;
    }
}
