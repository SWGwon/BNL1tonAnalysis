#include "Waveform.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

// Constructor: Reads waveform data from a UShort_t pointer and exits if a zero
// value is encountered.
Waveform::Waveform(const unsigned short *data, int maxSampleSize) {
    for (int i = 0; i < maxSampleSize; ++i) {
        if (data[i] != 0)
            samples_.push_back((double)data[i]);
        else
            break;
    }
}

// Constructor: Copies an int vector directly.
Waveform::Waveform(const std::vector<double> &samples) : samples_(samples) {}

// Getter: Returns waveform samples.
const std::vector<double> &Waveform::getSamples() const { return samples_; }

// Getter: Returns the PE-corrected waveform.
const std::vector<double> &Waveform::getAmpPE() const {
    if (!isAmpPeUpToDate_) {
        std::cerr << "Warning (Waveform::getAmpPE): PE waveform cache is not up-to-date. Call setAmpPE() to calculate it." << std::endl;
    }
    return amp_pe_;
}

void Waveform::setAmpPE(double spe_mean) {
    // getAmpMV() ensures that amp_mv_ is calculated and up-to-date
    const auto& ampMV = getAmpMV();

    amp_pe_.clear();
    if (spe_mean == 0.0) {
        std::cerr << "Warning (Waveform::setAmpPE): spe_mean is zero. PE waveform will be empty." << std::endl;
    } else {
        if (!ampMV.empty()) {
            amp_pe_.reserve(ampMV.size());
            for (const auto& value : ampMV) {
                amp_pe_.push_back(value / spe_mean);
            }
        }
    }
    
    // Mark the PE cache as up-to-date and store the spe_mean value used.
    isAmpPeUpToDate_ = true;
    cached_spe_mean_ = spe_mean;
}

// Static helper: Checks if there is any value below the threshold.
bool Waveform::hasValueLessThan(const std::vector<double> &vec, int threshold) {
    return std::any_of(vec.begin(), vec.end(),
                       [threshold](double value) { return value < threshold; });
}

int Waveform::computeFlatBaseline(int start, int end) const {
    if (start < 0 || end >= static_cast<int>(samples_.size()) || start > end) {
        std::cerr << "Error (Waveform::computeFlatBaseline): Invalid range. Start: " << start
                  << ", End: " << end << ", samples_.size(): " << samples_.size() << std::endl;
        // Or handle as per your project's error strategy, e.g., return a specific error code or 0
        // throw std::invalid_argument("Invalid range for baseline computation"); // Original behavior
        return 0; // Example: return 0 on invalid range if throw is not desired here
    }

    std::vector<double> subVector_double;
    // Reserve space if the range is valid and non-empty
    if (end >= start) {
        subVector_double.reserve(end - start + 1);
    }

    for (int i = start; i <= end; ++i) {
        double current_sample_val = samples_[i]; // Assuming samples_ stores or can be treated as double
        if (std::isnan(current_sample_val)) {
            std::cerr << "Debug (Waveform::computeFlatBaseline): samples_[" << i << "] is NaN. Value: "
                      << current_sample_val << std::endl;
        }
        subVector_double.push_back(current_sample_val);
    }

    if (subVector_double.empty()) {
         std::cerr << "Warning (Waveform::computeFlatBaseline): subVector for baseline computation is empty for range ["
                   << start << ", " << end << "]. Returning 0." << std::endl;
        return 0;
    }

    // std::sort behavior with NaNs: NaNs are typically moved to the end of the sorted range.
    // This affects median calculation if NaNs are present.
    std::sort(subVector_double.begin(), subVector_double.end());

    int n = subVector_double.size();
    double median_val_double;

    if (n % 2 == 0) {
        double val1 = subVector_double[n / 2 - 1];
        double val2 = subVector_double[n / 2];
        if (std::isnan(val1)) {
            std::cerr << "Debug (Waveform::computeFlatBaseline): Median component subVector_double[" << (n / 2 - 1) << "] (value: " << val1 << ") is NaN." << std::endl;
        }
        if (std::isnan(val2)) {
            std::cerr << "Debug (Waveform::computeFlatBaseline): Median component subVector_double[" << (n / 2) << "] (value: " << val2 << ") is NaN." << std::endl;
        }
        median_val_double = (val1 + val2) / 2.0;
    } else {
        median_val_double = subVector_double[n / 2];
        if (std::isnan(median_val_double)) {
            std::cerr << "Debug (Waveform::computeFlatBaseline): Median value subVector_double[" << (n / 2) << "] (value: " << median_val_double << ") is NaN." << std::endl;
        }
    }

    if (std::isnan(median_val_double)) {
        std::cerr << "Warning (Waveform::computeFlatBaseline): Calculated median (as double) is NaN for range ["
                  << start << ", " << end << "]." << std::endl;
        std::cerr << "  Returning 0 as int to avoid Undefined Behavior from NaN-to-int conversion." << std::endl;
        return 0; // Return a defined int value (e.g., 0) if the double median is NaN
    }
    
    // Convert non-NaN double median to int. Note potential precision loss.
    return static_cast<int>(median_val_double);
}

void Waveform::subtractFlatBaseline(int start, int end) {
    int baseline_int = computeFlatBaseline(start, end);
    double baseline = static_cast<double>(baseline_int); // Work with double for consistency

    int index = 0;
    for (double &sample : samples_) { // Ensure sample is a reference to modify original
        double original_sample_for_log = sample; // For logging original value
        bool was_nan_initially = std::isnan(sample);

        if (was_nan_initially) {
            std::cerr << "Debug (Waveform::subtractFlatBaseline): samples_[" << index << "] is NaN (" << sample
                      << ") before operations." << std::endl;
        }

        sample -= baseline;
        sample *= -1.0; // 파형을 양수화 (use -1.0 for clarity with doubles)

        bool is_nan_after_ops = std::isnan(sample);

        if (is_nan_after_ops) {
            if (!was_nan_initially) {
                // This would be unusual if 'baseline' is not NaN and 'original_sample_for_log' was finite.
                std::cerr << "Debug (Waveform::subtractFlatBaseline): samples_[" << index << "] BECAME NaN after operations." << std::endl;
            } else {
                std::cerr << "Debug (Waveform::subtractFlatBaseline): samples_[" << index << "] REMAINED NaN after operations." << std::endl;
            }
            std::cerr << "  Details: original_sample = " << original_sample_for_log
                      << ", baseline_subtracted = " << baseline
                      << ", final_sample = " << sample << std::endl;
        }
        index++;
    }
    this->isBaselineSubstracted_ = true;
    isAmpMvUpToDate_ = false;
    isAmpPeUpToDate_ = false;
}

const std::vector<double>& Waveform::getAmpMV() const {
    // 데이터가 최신이 아니면 다시 계산
    if (!isAmpMvUpToDate_) {
        // --- 이 부분이 예전 setAmpMV의 역할 ---
        amp_mv_.clear();
        if (!samples_.empty()) {
            amp_mv_.reserve(samples_.size());
            const double scale = 2000.0 / (std::pow(2, 14) - 1.0);
            for (const auto& sample : samples_) {
                amp_mv_.push_back(sample * scale / 50.0);
            }
        }
        // ------------------------------------
        isAmpMvUpToDate_ = true; // 계산했으므로 최신 상태로 플래그 변경
    }
    return this->amp_mv_;
}

bool Waveform::hasPeakAboveThreshold(double threshold) {
    const auto& ampMV = getAmpMV();
    if (ampMV.empty()) {
        return false;
    }
    double maxValue = *std::max_element(ampMV.begin(), ampMV.end());

    // Print the found maximum value.
    //std::cout << "Maximum peak value found: " << maxValue << " mV." << std::endl;

    return maxValue > threshold;
}

double Waveform::getCharge(int integrateStart, int integrateEnd) const {
    const auto& ampMV = getAmpMV();
    if (integrateStart < 0 || integrateEnd >= static_cast<int>(ampMV.size()) || integrateStart > integrateEnd) {
        std::cerr << "Error (Waveform::getPE): Invalid integration range." << std::endl;
        return 0.0;
    }

    double charge = 0.0;
    for (int i = integrateStart; i <= integrateEnd; ++i) {
        charge += ampMV[i];
    }

    return charge;
}
double Waveform::getPE(int integrateStart, int integrateEnd, double spe_mean, int sampleToNS) const {
    if (spe_mean == 0.0) {
        std::cerr << "Error (Waveform::getPE): spe_mean cannot be zero." << std::endl;
        return 0.0;
    }
    double charge = this->getCharge(integrateStart, integrateEnd);
    return charge * sampleToNS / spe_mean;
}

// Performs daisy chain trigger delay correction based on the channel name
// (ch_name). Correction is applied by appropriately slicing the amp_pe_ vector.
void Waveform::correctDaisyChainTrgDelay(const std::string &ch_name, int dt_ns,
                                         int sampleToNS) {
    int dS = dt_ns / sampleToNS;
    int start_offset = 0;
    int end_offset = 0;

    if (ch_name.find("_b1") != std::string::npos) {
        start_offset = dS * 3;
    } else if (ch_name.find("_b2") != std::string::npos) {
        start_offset = dS * 2;
        end_offset = dS;
    } else if (ch_name.find("_b3") != std::string::npos) {
        start_offset = dS;
        end_offset = dS * 2;
    } else if (ch_name.find("_b4") != std::string::npos || ch_name.find("_b5") !=
            std::string::npos) {
        end_offset = dS * 3;
    } else {
        std::cerr << "ERROR in correctDaisyChainTrgDelay: invalid boardId in channel name " <<
            ch_name << std::endl;
        return;
    }
    // 벡터를 슬라이싱하는 람다 함수
    auto slice_vector = [&](std::vector<double>& vec) {
        if (vec.size() <= static_cast<size_t>(start_offset + end_offset)) {
            std::cerr << "Not enough samples for correction for " << ch_name << std::endl;
            vec.clear(); // 에러 발생 시 벡터를 비움
            return;
        }
        // 새 벡터를 만드는 대신 기존 벡터를 수정하여 효율성 증대
        vec.erase(vec.begin(), vec.begin() + start_offset);
        vec.erase(vec.end() - end_offset, vec.end());
    };

    slice_vector(amp_pe_);
    slice_vector(samples_);
    isAmpMvUpToDate_ = false;
    isAmpPeUpToDate_ = false;

    /*
    if (amp_pe_.empty()) {
        std::cerr << "amp_pe_ is empty, cannot correct delay." << std::endl;
        return;
    }
    std::vector<double> corrected;
    std::vector<double> corrected_samples;
    if (ch_name.find("_b1") != std::string::npos) {
        if (amp_pe_.size() > static_cast<size_t>(dS * 3))
            corrected =
                std::vector<double>(amp_pe_.begin() + dS * 3, amp_pe_.end());
        else {
            std::cerr << "Not enough samples for _b1 correction" << std::endl;
            return;
        }
        if (samples_.size() > static_cast<size_t>(dS * 3))
            corrected_samples =
                std::vector<double>(samples_.begin() + dS * 3, samples_.end());
        else {
            std::cerr << "Not enough samples for _b1 correction" << std::endl;
            return;
        }
    } else if (ch_name.find("_b2") != std::string::npos) {
        if (amp_pe_.size() > static_cast<size_t>(dS * 3))
            corrected = std::vector<double>(amp_pe_.begin() + dS * 2,
                                            amp_pe_.end() - dS);
        else {
            std::cerr << "Not enough samples for _b2 correction" << std::endl;
            return;
        }
        if (samples_.size() > static_cast<size_t>(dS * 3))
            corrected_samples = std::vector<double>(samples_.begin() + dS * 2,
                                            samples_.end() - dS);
        else {
            std::cerr << "Not enough samples for _b2 correction" << std::endl;
            return;
        }
    } else if (ch_name.find("_b3") != std::string::npos) {
        if (amp_pe_.size() > static_cast<size_t>(dS * 3))
            corrected = std::vector<double>(amp_pe_.begin() + dS,
                                            amp_pe_.end() - dS * 2);
        else {
            std::cerr << "Not enough samples for _b3 correction" << std::endl;
            return;
        }
        if (samples_.size() > static_cast<size_t>(dS * 3))
            corrected_samples = std::vector<double>(samples_.begin() + dS,
                                            samples_.end() - dS * 2);
        else {
            std::cerr << "Not enough samples for _b3 correction" << std::endl;
            return;
        }
    } else if (ch_name.find("_b4") != std::string::npos ||
               ch_name.find("_b5") != std::string::npos) {
        if (amp_pe_.size() > static_cast<size_t>(dS * 3))
            corrected =
                std::vector<double>(amp_pe_.begin(), amp_pe_.end() - dS * 3);
        else {
            std::cerr << "Not enough samples for _b4/_b5 correction"
                      << std::endl;
            return;
        }
        if (samples_.size() > static_cast<size_t>(dS * 3))
            corrected_samples =
                std::vector<double>(samples_.begin(), samples_.end() - dS * 3);
        else {
            std::cerr << "Not enough samples for _b4/_b5 correction"
                      << std::endl;
            return;
        }
    } else {
        std::cerr << "ERROR in correctDaisyChainTrgDelay: invalid boardId in "
                     "channel name "
                  << ch_name << std::endl;
        return;
    }
    amp_pe_ = corrected;
    samples_ = corrected_samples;
    */
}

// Converts amp_pe_ data into a TGraph and returns it.
TGraph *Waveform::drawPEAsGraph(const std::string &title) const {
    // Use the getter to ensure data consistency
    const auto& amp_pe = getAmpPE();
    int n = amp_pe.size();
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = i;
    }
    TGraph *graph = new TGraph(n, x.data(), amp_pe.data());
    graph->SetTitle(title.c_str());
    return graph;
}

TGraph *Waveform::drawMVAsGraph(const std::string &title) const {
    int n = samples_.size();
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = i;
    }
    TGraph *graph = new TGraph(n, x.data(), samples_.data());
    graph->SetTitle(title.c_str());
    return graph;
}
