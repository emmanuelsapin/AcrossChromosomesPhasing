/**
 * @file PhasingProgram.h
 * @brief Main header file including all dependencies
 */

#ifndef PHASING_PROGRAM_H
#define PHASING_PROGRAM_H

// Core dependencies
#include "Constants.h"
#include "ErrorCodes.h"
#include "Exceptions.h"
#include "Interfaces.h"
#include "ChromosomeDivider.h"

// Utility functions
#include "utils/FileIOUtils.h"

// Data management
#include "GenomeDataManager.h"
#include "RelativeIdentificationEngine.h"

// I/O operations
#include "GenomeFileLoader.h"
#include "OutputFileWriter.h"

// Configuration
#include "ConfigurationManager.h"

// Main program
#include "HaplotypePhasingProgram.h"

#endif // PHASING_PROGRAM_H

