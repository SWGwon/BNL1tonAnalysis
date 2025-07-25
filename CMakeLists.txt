cmake_minimum_required(VERSION 3.10)

# You can turn this on by running: cmake -DENABLE_ASAN=ON ..
option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
project(Process1tonData CXX)

if(ENABLE_ASAN)
  message(STATUS "AddressSanitizer (ASan) is being force-enabled.")
  # We know AppleClang supports ASan, so just add the flags directly.
  add_compile_options(-fsanitize=address -g -O1 -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address)
endif()

set(CMAKE_CXX_STANDARD 17)

include_directories(${PROJECT_SOURCE_DIR}/eigen-3.4.0)
include_directories(${PROJECT_SOURCE_DIR}/include)

find_package(ROOT REQUIRED COMPONENTS RIO Net Hist Tree TreePlayer MathCore Thread)
find_package(Threads REQUIRED)

add_executable(SelectCrossingMuon src/SelectCrossingMuon.cpp)
target_include_directories(SelectCrossingMuon PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(SelectCrossingMuon PUBLIC ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::MathCore ROOT::Thread Threads::Threads)

add_executable(MakeResultRootFile src/MakeResultRootFile.cpp)
target_include_directories(MakeResultRootFile PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(MakeResultRootFile PUBLIC ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::MathCore ROOT::Thread Threads::Threads)
target_link_libraries(MakeResultRootFile PUBLIC Threads::Threads)

add_executable(MakeCrossingMuonRootFile src/MakeCrossingMuonRootFile.cpp)
target_include_directories(MakeCrossingMuonRootFile PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(MakeCrossingMuonRootFile PUBLIC ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::MathCore ROOT::Thread Threads::Threads)
target_link_libraries(MakeCrossingMuonRootFile PUBLIC Threads::Threads)

add_executable(ProcessMCFile src/ProcessMCFile.cpp)
target_include_directories(ProcessMCFile PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ProcessMCFile PUBLIC ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::MathCore ROOT::Thread Threads::Threads)
target_link_libraries(ProcessMCFile PUBLIC Threads::Threads)

add_executable(ProcessRawRoot src/ProcessRawRoot.cpp src/Waveform.cpp
    src/DataProcessor.cpp src/AppConfig.cpp)
target_include_directories(ProcessRawRoot PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ProcessRawRoot PUBLIC ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::MathCore ROOT::Thread Threads::Threads)

add_executable(Process1tonData src/Process1tonData.cpp src/Waveform.cpp
    src/DataProcessor1ton.cpp src/AppConfig.cpp)
target_include_directories(Process1tonData PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(Process1tonData PUBLIC ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::MathCore ROOT::Thread Threads::Threads)

add_executable(DailyCheck src/DailyCheck.cpp src/Waveform.cpp
    src/DataProcessor.cpp src/AppConfig.cpp)
target_include_directories(DailyCheck PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(DailyCheck PUBLIC ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::MathCore ROOT::Thread Threads::Threads)

add_executable(DailyCheck30T src/DailyCheck30T.cpp src/Waveform.cpp
    src/DataProcessor.cpp src/AppConfig.cpp)
target_include_directories(DailyCheck30T PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(DailyCheck30T PUBLIC ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::MathCore ROOT::Thread Threads::Threads)

add_executable(ConvertMCtoDataFormat src/ConvertMCtoDataFormat.cpp
    src/AppConfig.cpp)
target_include_directories(ConvertMCtoDataFormat PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(ConvertMCtoDataFormat PUBLIC ROOT::RIO ROOT::Net ROOT::Hist ROOT::Tree ROOT::TreePlayer ROOT::MathCore ROOT::Thread Threads::Threads)

#find_package(OSQP REQUIRED)
#find_package(OsqpEigen REQUIRED)
#
#include_directories(/usr/local/include/qpOASES)
#link_directories(/usr/local/lib)
#
#add_executable(BNL1tonAnalysis src/main.cpp)
#
#target_include_directories(BNL1tonAnalysis PUBLIC
#    ${ROOT_INCLUDE_DIRS}
#    ${PROJECT_SOURCE_DIR}/eigen-3.4.0
#    ${OSQP_INCLUDE_DIRS}
#    /usr/local/include/qpOASES
#)
#
#target_link_libraries(BNL1tonAnalysis PUBLIC
#    ${ROOT_LIBRARIES}
#    osqp::osqp
#    OsqpEigen::OsqpEigen
#    qpOASES
#)
