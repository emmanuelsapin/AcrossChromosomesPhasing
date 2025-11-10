/**
 * @file ChromosomeDivider.h
 * @brief Chromosomal segment representation with phasing information
 */

#ifndef CHROMOSOME_DIVIDER_H
#define CHROMOSOME_DIVIDER_H

namespace PhasingEngine {

/**
 * @struct ChromosomeDivider
 * @brief Represents a chromosomal segment with phasing information
 */
struct ChromosomeDivider {
    int start;          ///< Starting SNP index
    int end;            ///< Ending SNP index
    int phasing;        ///< Phasing orientation (1 or 2)
    int nbright;        ///< Number of correctly phased positions
    int nbwrong;        ///< Number of incorrectly phased positions
    int segment;        ///< Segment identifier
    
    ChromosomeDivider() : start(0), end(0), phasing(0), nbright(0), nbwrong(0), segment(0) {}
    
    /**
     * @brief Calculate phasing accuracy percentage
     * @return Accuracy as percentage (0-100)
     */
    double calculateAccuracy() const {
        int total = nbright + nbwrong;
        return total > 0 ? (100.0 * nbright) / total : 0.0;
    }
    
    /**
     * @brief Check if segment is valid
     * @return true if start < end
     */
    bool isValid() const { return start < end && start >= 0; }
    
    /**
     * @brief Get segment length in SNPs
     * @return Number of SNPs in segment
     */
    int getLength() const { return end - start; }
};

}

#endif // CHROMOSOME_DIVIDER_H


