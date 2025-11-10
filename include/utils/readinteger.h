/**
 * @file readinteger.h
 * @brief Utility function to read positive integers from files
 * @details Reads an integer from a file stream, skipping non-numeric characters.
 *          Returns INT32_MAX if EOF is reached.
 */

#ifndef READ_INTEGER_H
#define READ_INTEGER_H

#include <inttypes.h>
#include <stdio.h>

#define INT int32_t
#define MAXINT INT32_MAX

/**
 * @brief Read a positive integer from a file
 * @param file File pointer to read from
 * @return The integer read, or INT32_MAX if EOF is reached
 */
INT readinteger(FILE * file)
{	
    INT ID = 0;
    char carac;
    do 
    {	
        carac = getc(file);
        if (carac != EOF && carac > 47 && carac < 58) 
            ID = ID * 10 + carac - 48;
    } while (carac != EOF && carac > 47 && carac < 58);
    
    if (carac == EOF) 
        ID = MAXINT;
    
    return(ID); 	
}

#endif // READ_INTEGER_H


