/**
 * @file OutputFileWriter.h
 * @brief Professional output generation with multiple format support
 */

#ifndef OUTPUT_FILE_WRITER_H
#define OUTPUT_FILE_WRITER_H

#include <string>
#include <vector>

namespace PhasingEngine {

class GenomeDataManager;

/**
 * @class OutputFileWriter
 * @brief Professional output generation with multiple format support
 */
class OutputFileWriter {
private:
    GenomeDataManager* genomeDataManager;
    std::string outputDirectory;
    bool validateOutputPath(const char* path) const;
    void writePEDHeader(FILE* file, int chromosome) const;
    void writePEDRow(FILE* file, int individualID, int chromosome) const;
    
public:
    explicit OutputFileWriter(GenomeDataManager* gdm);
    virtual ~OutputFileWriter() = default;
    
    bool writeOutput(int chromosome, const char* outputPath);
    bool writeOutputBatch(const std::vector<int>& chromosomes, const char* outputPath);
    void setOutputDirectory(const std::string& directory) { outputDirectory = directory; }
    std::string getOutputDirectory() const { return outputDirectory; }
};

}

#endif // OUTPUT_FILE_WRITER_H


