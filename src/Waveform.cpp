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
const std::vector<double> &Waveform::getAmpPE() const { return amp_pe_; }

// Static helper: Checks if there is any value below the threshold.
bool Waveform::hasValueLessThan(const std::vector<double> &vec, int threshold) {
    return std::any_of(vec.begin(), vec.end(),
                       [threshold](double value) { return value < threshold; });
}

// Internal helper: Calculates the median value from a given range of samples
// and uses it as the baseline.
int Waveform::computeFlatBaseline(int start, int end) const {
    if (start < 0 || end >= static_cast<int>(samples_.size()) || start > end) {
        throw std::invalid_argument("Invalid range for baseline computation");
    }
    std::vector<int> subVector(samples_.begin() + start,
                               samples_.begin() + end + 1);
    std::sort(subVector.begin(), subVector.end());
    int n = subVector.size();
    if (n % 2 == 0) {
        return (subVector[n / 2 - 1] + subVector[n / 2]) / 2.;
    } else {
        return subVector[n / 2];
    }
}

// Baseline correction: Subtracts the baseline value from each sample and
// multiplies the result by -1.
void Waveform::subtractFlatBaseline(int start, int end) {
    int baseline = computeFlatBaseline(start, end);
    for (auto &sample : samples_) {
        sample -= baseline;
        sample *= -1; // 파형을 양수화
    }
}

// Converts waveform samples to mV units.
std::vector<double> Waveform::getAmpMV(const std::string & /*ch_name*/) const {
    std::vector<double> ampMV;
    ampMV.reserve(samples_.size());
    const double scale = 2000.0 / (std::pow(2, 14) - 1);
    for (auto sample : samples_) {
        ampMV.push_back(sample * scale);
    }
    return ampMV;
}

// Converts mV data to PE units.
// Conversion formula: amp_pe = amp_mV / factor / spe_mean
// factor: 50 ohm
void Waveform::setAmpPE(double spe_mean, double factor) {
    std::vector<double> ampMV = getAmpMV();
    amp_pe_.clear();
    amp_pe_.reserve(ampMV.size());
    for (auto value : ampMV) {
        double temp = value / factor / spe_mean;
        amp_pe_.push_back(temp);
    }
}

// Integrates a specific range of the amp_pe_ vector and returns the PE value.
//double Waveform::getPE(int integrateStart, int integrateEnd,
//                       int sampleToNS) const {
//    if (integrateStart < 0 ||
//        integrateEnd >= static_cast<int>(amp_pe_.size()) ||
//        integrateStart > integrateEnd) {
//        std::cerr << "Error: Invalid range for getPE. Start: " << integrateStart
//                  << ", End: " << integrateEnd << std::endl;
//        return 0;
//    }
//    double total = 0;
//    for (int i = integrateStart; i <= integrateEnd; ++i) {
//        total += amp_pe_[i] * sampleToNS;
//    }
//    return total;
//}
double Waveform::getPE(int integrateStart, int integrateEnd,
                       int sampleToNS) const {
    // 범위 유효성 검사
    if (integrateStart < 0 ||
        integrateEnd >= static_cast<int>(amp_pe_.size()) ||
        integrateStart > integrateEnd) {
        std::cerr << "Error (Waveform::getPE): Invalid range. Start: " << integrateStart
                  << ", End: " << integrateEnd << ", amp_pe_.size(): " << amp_pe_.size() << std::endl;
        return 0.0; // 함수 반환 타입이 double이므로 0.0 반환
    }

    double total = 0.0; // double 타입 명시를 위해 0.0으로 초기화
    bool total_became_nan_logged = false; // total이 NaN이 되었을 때 로그를 한 번만 출력하기 위한 플래그

    for (int i = integrateStart; i <= integrateEnd; ++i) {
        double current_amp_pe = amp_pe_[i]; // 현재 샘플 값

        // 선택적 로그: 입력 데이터(amp_pe_[i]) 자체가 NaN인 경우
        if (std::isnan(current_amp_pe)) {
            std::cerr << "Debug (Waveform::getPE): Input amp_pe_[" << i << "] is NaN." << std::endl;
        }

        double term = current_amp_pe * sampleToNS; // 현재 항 계산

        // 선택적 로그: 계산된 항(term)이 NaN인 경우 (예: infinity * 0)
        // (입력 current_amp_pe가 NaN이 아니었을 때만 이 로그가 특별히 의미 있을 수 있음)
        if (std::isnan(term) && !std::isnan(current_amp_pe)) {
             std::cerr << "Debug (Waveform::getPE): Calculated term (amp_pe_[" << i << "] * sampleToNS) is NaN at index " << i << "." << std::endl;
             std::cerr << "  Details: amp_pe_[" << i << "] = " << current_amp_pe
                       << ", sampleToNS = " << sampleToNS << std::endl;
        }
        
        double total_before_addition = total; // 더하기 전의 total 값 저장
        total += term; // total에 현재 항 더하기

        // total이 NaN이 되었는지 확인하고, 아직 로그를 남기지 않았다면 상세 정보 출력
        if (!total_became_nan_logged && std::isnan(total)) {
            std::cerr << "Debug (Waveform::getPE): 'total' became NaN at loop index i = " << i << "." << std::endl;
            std::cerr << "  Function arguments: integrateStart = " << integrateStart
                      << ", integrateEnd = " << integrateEnd
                      << ", sampleToNS = " << sampleToNS << std::endl;
            std::cerr << "  Iteration details (i = " << i << "):" << std::endl;
            std::cerr << "    amp_pe_[" << i << "] = " << current_amp_pe << std::endl;
            std::cerr << "    term (amp_pe_[" << i << "] * sampleToNS) = " << term << std::endl;
            std::cerr << "    total before this addition = " << total_before_addition << std::endl;
            std::cerr << "    total after this addition = " << total << " (is NaN)" << std::endl;
            total_became_nan_logged = true; // 로그 남겼음을 표시

            // NaN이 발생한 이후의 계산은 의미가 없을 수 있으므로, 여기서 루프를 중단할 수 있습니다.
            // 예를 들어: break;
        }
    }
    return total;
}

// Performs daisy chain trigger delay correction based on the channel name
// (ch_name). Correction is applied by appropriately slicing the amp_pe_ vector.
void Waveform::correctDaisyChainTrgDelay(const std::string &ch_name, int dt_ns,
                                         int sampleToNS) {
    int dS = dt_ns / sampleToNS;
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
}

// Converts amp_pe_ data into a TGraph and returns it.
TGraph *Waveform::drawPEAsGraph(const std::string &title) const {
    int n = amp_pe_.size();
    std::vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = i;
    }
    TGraph *graph = new TGraph(n, x.data(), amp_pe_.data());
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
