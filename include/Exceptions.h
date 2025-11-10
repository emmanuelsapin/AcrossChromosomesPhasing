/**
 * @file Exceptions.h
 * @brief Custom exception classes for phasing operations
 */

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include "ErrorCodes.h"
#include <stdexcept>
#include <string>

namespace PhasingEngine {

/**
 * @class PhasingException
 * @brief Custom exception class for phasing operations
 */
class PhasingException : public std::runtime_error {
private:
    ErrorCodes::PhasingError errorCode;
    
public:
    PhasingException(const std::string& message, 
                    ErrorCodes::PhasingError code = ErrorCodes::PhasingError::INVALID_INPUT)
        : std::runtime_error(message), errorCode(code) {}
    
    ErrorCodes::PhasingError getErrorCode() const { return errorCode; }
};

}

#endif // EXCEPTIONS_H


