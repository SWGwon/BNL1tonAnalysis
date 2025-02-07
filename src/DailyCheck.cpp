#include "AppConfig.h"
#include "DataProcessor.h"
#include <chrono>
#include <iostream>

int main(int argc, char** argv) {
    auto start = std::chrono::high_resolution_clock::now();

    AppConfig config = AppConfig::parseArgs(argc, argv);
    std::vector<std::string> pmts_all_1ton = {"adc_b1_ch0",  "adc_b1_ch1",  "adc_b1_ch2",  "adc_b1_ch3",
                "adc_b1_ch4",  "adc_b1_ch5",  "adc_b1_ch6",  "adc_b1_ch7",
                "adc_b1_ch8",  "adc_b1_ch9",  "adc_b1_ch10", "adc_b1_ch11",
                "adc_b1_ch12", "adc_b1_ch13", "adc_b1_ch14", "adc_b1_ch15",
                "adc_b2_ch0",  "adc_b2_ch1",  "adc_b2_ch2",  "adc_b2_ch3",
                "adc_b2_ch4",  "adc_b2_ch5",  "adc_b2_ch6",  "adc_b2_ch7",
                "adc_b2_ch8",  "adc_b2_ch9",  "adc_b2_ch10", "adc_b2_ch11",
                "adc_b2_ch12", "adc_b2_ch13", "adc_b2_ch14", "adc_b2_ch15",
                "adc_b3_ch0",  "adc_b3_ch1",  "adc_b3_ch2",  "adc_b3_ch3",
                "adc_b3_ch4",  "adc_b3_ch5",  "adc_b3_ch6",  "adc_b3_ch7",
                "adc_b3_ch8",  "adc_b3_ch9",  "adc_b3_ch10", "adc_b3_ch11",
                "adc_b3_ch12", "adc_b3_ch13", "adc_b3_ch14", "adc_b3_ch15",
                "adc_b4_ch0",  "adc_b4_ch1",  "adc_b4_ch2",  "adc_b4_ch3",
                "adc_b4_ch4",  "adc_b4_ch5",  "adc_b4_ch6",  "adc_b4_ch7",
                "adc_b4_ch8",  "adc_b4_ch9",  "adc_b4_ch10", "adc_b4_ch11",
                "adc_b4_ch12", "adc_b4_ch13", "adc_b4_ch14", "adc_b4_ch15",
                "adc_b5_ch0",  "adc_b5_ch1",  "adc_b5_ch2",  "adc_b5_ch3",
                "adc_b5_ch4",  "adc_b5_ch5",  "adc_b5_ch6",  "adc_b5_ch7",
                "adc_b5_ch8",  "adc_b5_ch9",  "adc_b5_ch10", "adc_b5_ch11",
                "adc_b5_ch12", "adc_b5_ch13", "adc_b5_ch14", "adc_b5_ch15",
                "adc_b5_ch16", "adc_b5_ch17", "adc_b5_ch18", "adc_b5_ch19",
                "adc_b5_ch20", "adc_b5_ch21", "adc_b5_ch22", "adc_b5_ch23",
                "adc_b5_ch24", "adc_b5_ch25", "adc_b5_ch26", "adc_b5_ch27",
                "adc_b5_ch28", "adc_b5_ch29", "adc_b5_ch30", "adc_b5_ch31",
                "adc_b5_ch32"};

    DataProcessor processor(config);
    processor.setPMTs(pmts_all_1ton);
    processor.dailyCheck();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;
    return 0;
}
