#ifndef APPCONFIG_H
#define APPCONFIG_H

#include <string>

class AppConfig {
public:
    std::string inputFileName;
    std::string outputFilePath;
    std::string inputSPECalibrationPath;
    int eventNumber;
    int triggerType;

    AppConfig();

    static AppConfig parseArgs(int argc, char** argv);
};

#endif // APPCONFIG_H
