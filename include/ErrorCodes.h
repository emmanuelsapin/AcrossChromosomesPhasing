/**
 * @file ErrorCodes.h
 * @brief Error code enumeration for exception handling
 */

#ifndef ERROR_CODES_H
#define ERROR_CODES_H

namespace PhasingEngine {
namespace ErrorCodes {
    enum class PhasingError {
        SUCCESS = 0,
        FILE_NOT_FOUND = 1,
        INVALID_INPUT = 2,
        MEMORY_ALLOCATION_FAILED = 3,
        INVALID_CHROMOSOME = 4,
        INVALID_INDIVIDUAL = 5,
        INVALID_SNP_INDEX = 6,
        DATA_CORRUPTION = 7
    };
}
}

#endif // ERROR_CODES_H


