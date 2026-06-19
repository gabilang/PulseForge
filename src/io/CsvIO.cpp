#include "bsp/io/CsvIO.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace bsp {
namespace io {

std::vector<double> loadColumn(const std::string& path, std::size_t column) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open file: " + path);

    std::vector<double> values;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        // Split on comma, tab, or whitespace.
        for (char& c : line) {
            if (c == ',' || c == '\t' || c == ';') c = ' ';
        }
        std::istringstream ss(line);
        std::string tok;
        std::size_t col = 0;
        bool got = false;
        double val = 0.0;
        while (ss >> tok) {
            if (col == column) {
                try {
                    std::size_t pos = 0;
                    val = std::stod(tok, &pos);
                    if (pos == tok.size()) got = true;
                } catch (...) {
                    got = false;  // header / non-numeric -> skip line
                }
                break;
            }
            ++col;
        }
        if (got) values.push_back(val);
    }
    return values;
}

Signal loadSignal(const std::string& path, double fs, std::size_t column, const std::string& label) {
    return Signal(loadColumn(path, column), fs, label);
}

void saveSignal(const std::string& path, const Signal& sig) {
    // One sample per line (the time axis is implicit from the sample rate).
    // This is the canonical format that loadColumn() reads back with col 0.
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write file: " + path);
    out << "value\n";
    for (std::size_t i = 0; i < sig.size(); ++i)
        out << sig.samples[i] << "\n";
}

void saveColumns(const std::string& path, const std::vector<std::string>& headers,
                 const std::vector<std::vector<double>>& columns) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write file: " + path);
    for (std::size_t i = 0; i < headers.size(); ++i) {
        out << headers[i];
        if (i + 1 < headers.size()) out << ",";
    }
    out << "\n";
    std::size_t rows = columns.empty() ? 0 : columns[0].size();
    for (std::size_t r = 0; r < rows; ++r) {
        for (std::size_t c = 0; c < columns.size(); ++c) {
            out << (r < columns[c].size() ? columns[c][r] : 0.0);
            if (c + 1 < columns.size()) out << ",";
        }
        out << "\n";
    }
}

} // namespace io
} // namespace bsp
