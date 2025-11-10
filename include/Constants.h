/**
 * @file Constants.h
 * @brief System-wide constants and configuration parameters
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace PhasingEngine {
namespace Constants {
    constexpr int NSNPPERCHR = 60000;
    constexpr int MAXPOP = 435188;
    constexpr int MAXCLOSERELAT = 6000;
    constexpr int MAXCLOSERELATTEMP = 450000;
    constexpr int MAXNBDIVISOR = 25;
    constexpr int NBINDIVMAX = 100000;
    constexpr int NBINDIV = 100000;
    constexpr int NBINDIVEA = 1;
    constexpr int MAXGEN = 1;
    constexpr int MAXBREAK = 1000;
    constexpr int NUM_CHROMOSOMES = 23;
    constexpr float DEFAULT_PIHAT_THRESHOLD = 0.33f;
    constexpr float PIHAT_NORMALIZATION_FACTOR = 330005.0f;
}
}

#endif // CONSTANTS_H


