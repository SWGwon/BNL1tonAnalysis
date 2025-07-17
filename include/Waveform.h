#ifndef WAVEFORM_H
#define WAVEFORM_H

#include <vector>
#include <string>
#include <TGraph.h>

// The Waveform class stores raw waveform data and provides processing methods
// such as baseline correction, mV conversion, photoelectron (PE) calculation,
// and trigger delay correction.
class Waveform {
public:
    // Default constants (modifiable if needed)
    static const int DEFAULT_MAX_SAMPLE_SIZE = 3000;
    static const int DEFAULT_SAMPLE_TO_NS = 2; // Conversion factor: samples to nanoseconds
    static const int DEFAULT_DT_NS = 48;       // Used for daisy chain correction
    static const int BASELINE_END_BIN = 100;

    // Constructor: Creates a waveform from a ROOT UShort_t array
    // (e.g., TTree branch data). Stops if a zero value is encountered.
    Waveform(const unsigned short* data, int maxSampleSize = DEFAULT_MAX_SAMPLE_SIZE);

    // Constructor: Creates a waveform from an int vector.
    Waveform(const std::vector<double>& samples);
    
    // Constructor: Creates a waveform from an int vector.
    Waveform(const std::vector<double>& samples, std::string ch_name, double spe_mean);

    // Getter: Returns the current waveform samples (either raw or baseline-corrected).
    const std::vector<double>& getSamples() const;

    // Computes a flat baseline (default: median of the range 0~100),
    // subtracts it from each sample, and multiplies the result by -1 to make the waveform positive.
    void subtractFlatBaseline(int start = 0, int end = BASELINE_END_BIN);

    // Converts the current waveform samples to mV units.
    const std::vector<double>& getAmpMV() const;

    // Converts the waveform from mV units to photoelectron (PE) units
    // using spe_mean (Single Photon Response) correction.
    // Calculation formula: amp_pe = amp_mV / factor / spe_mean (default factor = 50)
    void setAmpPE(double spe_mean);
    // Getter: Returns the PE-converted data (after calling setAmpPE()).
    const std::vector<double>& getAmpPE() const;

    // Integrates a specific range (integrateStart ~ integrateEnd, inclusive) of the PE waveform (amp_pe_).
    // Each sample is multiplied by sampleToNS (default 2) before returning the result.
    double getPE(int integrateStart, int integrateEnd, double spe_mean, int sampleToNS = DEFAULT_SAMPLE_TO_NS) const;
    double getCharge(int integrateStart, int integrateEnd) const;

    // Performs daisy chain trigger delay correction.
    // The waveform is appropriately sliced based on the channel name (ch_name).
    void correctDaisyChainTrgDelay(const std::string& ch_name, int dt_ns = DEFAULT_DT_NS, int sampleToNS = DEFAULT_SAMPLE_TO_NS);

    // Converts the amp_pe_ vector into a TGraph for visualization.
    TGraph* drawPEAsGraph(const std::string& title) const;

    TGraph* drawMVAsGraph(const std::string& title) const;

    // Static helper: Checks if the given vector contains any value smaller than the threshold.
    static bool hasValueLessThan(const std::vector<double>& vec, int threshold);
    // Computes the baseline as the median of a given range [start, end] (inclusive).
    int computeFlatBaseline(int start, int end) const;

    void setName(std::string ch_name) {this->ch_name_ = ch_name;};
    std::string getName() {return this->ch_name_;};

    bool hasPeakAboveThreshold(double threshold);

private:
    //raw digitized waveform
    std::vector<double> samples_;

    mutable std::vector<double> amp_mv_;
    mutable std::vector<double> amp_pe_;

    mutable bool isAmpMvUpToDate_ = false;
    mutable bool isAmpPeUpToDate_ = false;
    mutable double cached_spe_mean_ = 0.0;

    bool isBaselineSubstracted_ = false;

    std::string ch_name_;
};

#endif // WAVEFORM_H
