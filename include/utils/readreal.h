/**
 * @file readreal.h
 * @brief Utility function to read positive real numbers from files
 * @details Reads a floating point number from a file stream, supporting decimal notation.
 *          Returns FLT_MAX if EOF is reached.
 */

#ifndef READ_REAL_H
#define READ_REAL_H

#include <inttypes.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>

#define INT int32_t
#define MAXINT INT32_MAX

/**
 * @brief Read a positive real number from a file
 * @param file File pointer to read from
 * @return The floating point number read, or FLT_MAX if EOF is reached
 */
float readreal(FILE * file)
{	
    float ID = 0;
    char carac;
    do 
    {	
        carac = getc(file);
        if (carac != EOF && carac > 47 && carac < 58) 
            ID = ID * 10 + carac - 48;
    } while (carac != EOF && carac > 47 && carac < 58);
    
    if (carac == EOF) 
        ID = FLT_MAX;
    else
    {	
        if (carac == '.')
        {	
            float dec = 10;
            do 
            {	
                carac = getc(file);
                if (carac != EOF && carac > 47 && carac < 58) 
                    ID = ID + (carac - 48) / dec;
                dec = dec * 10;
            } while (carac != EOF && carac > 47 && carac < 58);
        }
    }
    
    if (carac == EOF) 
        ID = FLT_MAX;
    
    return(ID); 	
}

#endif // READ_REAL_H


