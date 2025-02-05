#include "AppConfig.h"
#include "DataProcessor.h"
#include <chrono>
#include <iostream>

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();

    AppConfig config = AppConfig::parseArgs(argc, argv);

    DataProcessor processor(config);
    processor.dailyCheck();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;
    return 0;
}
