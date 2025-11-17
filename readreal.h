#ifndef READREAL_H
#define READREAL_H

#include <inttypes.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>

#define INT int32_t
#define MAXINT INT32_MAX

//
//  Read a positive real number from the file
//  Returns FLT_MAX if EOF
//	
float readreal(FILE * file);

#endif