#include "AppConfig.h"
#include "DataProcessor.h"
#include <chrono>
#include <iostream>

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();

    AppConfig config = AppConfig::parseArgs(argc, argv);

    std::vector<std::string> pmts_1ton = {"adc_b1_ch1",  "adc_b1_ch2",  "adc_b1_ch3",  "adc_b1_ch4",
             "adc_b1_ch5",  "adc_b1_ch6",  "adc_b1_ch7",  "adc_b1_ch8",
             "adc_b1_ch9",  "adc_b1_ch10", "adc_b1_ch11", "adc_b1_ch12",
             "adc_b1_ch13", "adc_b1_ch14", "adc_b1_ch15", "adc_b2_ch0",
             "adc_b2_ch1",  "adc_b2_ch2",  "adc_b2_ch3",  "adc_b2_ch4",
             "adc_b2_ch5",  "adc_b2_ch6",  "adc_b2_ch7",  "adc_b2_ch8",
             "adc_b2_ch9",  "adc_b2_ch10", "adc_b2_ch11", "adc_b2_ch12",
             "adc_b2_ch13", "adc_b2_ch14", "adc_b3_ch0",  "adc_b3_ch1",
             "adc_b3_ch2",  "adc_b3_ch3",  "adc_b3_ch4",  "adc_b3_ch5",
             "adc_b3_ch6",  "adc_b3_ch7",  "adc_b3_ch8",  "adc_b3_ch9",
             "adc_b3_ch10", "adc_b3_ch11", "adc_b3_ch12", "adc_b3_ch13",
             "adc_b3_ch14", "adc_b3_ch15", "adc_b4_ch0",  "adc_b4_ch1",
             "adc_b4_ch2",  "adc_b4_ch3",  "adc_b4_ch4",  "adc_b4_ch5",
             "adc_b4_ch6",  "adc_b4_ch7",  "adc_b4_ch8",  "adc_b4_ch9",
             "adc_b4_ch10", "adc_b4_ch11"};
    std::vector<std::string> triggers_1ton = {"adc_b5_ch33", "adc_b5_ch34", "adc_b5_ch35"};

    DataProcessor processor(config);
    processor.setPMTs(pmts_1ton);
    processor.setTriggers(triggers_1ton);
    processor.run();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;
    return 0;
}
