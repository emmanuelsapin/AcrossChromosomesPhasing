#ifndef READNEGATIVEREAL_H
#define READNEGATIVEREAL_H

#include <inttypes.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <math.h>   

//
//  Read a real number (possibly negative) from the file
//  Returns FLT_MAX if EOF
//	
double readnegativereal(FILE * file);

#endif