#ifndef READINTEGER_H
#define READINTEGER_H

#include <inttypes.h>
#include <stdio.h>

#define INT int32_t
#define MAXINT INT32_MAX

//
//  Read a positive integer from the file
//  Returns INT32_MAX if EOF
//	
INT readinteger(FILE * file);

#endif