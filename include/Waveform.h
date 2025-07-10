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

    // Constructor: Creates a waveform from a ROOT UShort_t array
    // (e.g., TTree branch data). Stops if a zero value is encountered.
    Waveform(const unsigned short* data, int maxSampleSize = DEFAULT_MAX_SAMPLE_SIZE);

    // Constructor: Creates a waveform from an int vector.
    Waveform(const std::vector<double>& samples);

    // Getter: Returns the current waveform samples (either raw or baseline-corrected).
    const std::vector<double>& getSamples() const;

    // Getter: Returns the PE-converted data (after calling setAmpPE()).
    const std::vector<double>& getAmpPE() const;

    // Computes a flat baseline (default: median of the range 0~100),
    // subtracts it from each sample, and multiplies the result by -1 to make the waveform positive.
    void subtractFlatBaseline(int start = 0, int end = 100);

    // Converts the current waveform samples to mV units.
    // Conversion formula: mV = sample * (2000 / (2^14 - 1)) / 50
    std::vector<double> getAmpMV(const std::string& /*ch_name*/ = "") const;

    // Converts the waveform from mV units to photoelectron (PE) units
    // using spe_mean (Single Photon Response) correction.
    // Calculation formula: amp_pe = amp_mV / factor / spe_mean (default factor = 50)
    void setAmpPE(double spe_mean);

    // Integrates a specific range (integrateStart ~ integrateEnd, inclusive) of the PE waveform (amp_pe_).
    // Each sample is multiplied by sampleToNS (default 2) before returning the result.
    double getPE(int integrateStart, int integrateEnd, int sampleToNS = DEFAULT_SAMPLE_TO_NS) const;

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

private:
    std::string ch_name_;
    std::vector<double> samples_;         // Raw or baseline-corrected waveform samples
    std::vector<double> amp_pe_;       // PE-corrected waveform (after calling setAmpPE())
};

#endif // WAVEFORM_H
