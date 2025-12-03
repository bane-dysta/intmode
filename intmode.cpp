#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <regex>
#include <algorithm>
#include <iomanip>
#include <cmath>

// ==================== 数据结构定义 ====================

// 振动模式数据
struct VibrationalMode {
    int mode_number;
    double frequency;      // cm^-1
    double huang_rhys;     // Huang-Rhys factor
    double lambda;         // Reorganization energy contribution (cm^-1)
};

// 内坐标对振动模式的贡献
struct CoordinateContribution {
    std::string name;           // 坐标名称，如 R1, A1, D1
    std::string definition;     // 坐标定义，如 R(1,2), A(2,1,3)
    double weight;              // 相对权重(%)
};

// 单个振动模式的内坐标贡献
struct ModeCoordinateData {
    int mode_number;
    std::vector<CoordinateContribution> contributions;
};

// ==================== 工具函数 ====================

std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string& str) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

bool isNumber(const std::string& str) {
    if (str.empty()) return false;
    size_t start = 0;
    if (str[0] == '-' || str[0] == '+') start = 1;
    if (start >= str.length()) return false;
    
    bool hasDot = false;
    for (size_t i = start; i < str.length(); i++) {
        if (str[i] == '.') {
            if (hasDot) return false;
            hasDot = true;
        } else if (!isdigit(str[i])) {
            return false;
        }
    }
    return true;
}

// 获取坐标类型（R, A, D等）
std::string getCoordinateType(const std::string& coordName) {
    if (coordName.empty()) return "";
    std::string type;
    for (char c : coordName) {
        if (isalpha(c)) {
            type += c;
        } else {
            break;
        }
    }
    return type;
}

// ==================== 重组能分析类 ====================

class ReorganizationEnergyAnalyzer {
private:
    std::vector<VibrationalMode> modes_;
    std::vector<ModeCoordinateData> coordinate_data_;
    std::string filename_;
    
    // 解析频率
    void parseFrequencies(std::ifstream& file) {
        std::string line;
        std::regex freq_pattern(R"(Frequencies\s+--\s+(.+))");
        int mode_count = 0;
        
        file.clear();
        file.seekg(0);
        
        while (std::getline(file, line)) {
            std::smatch match;
            if (std::regex_search(line, match, freq_pattern)) {
                std::istringstream iss(match[1].str());
                double freq;
                while (iss >> freq) {
                    VibrationalMode mode;
                    mode.mode_number = ++mode_count;
                    mode.frequency = freq;
                    mode.huang_rhys = 0.0;
                    mode.lambda = 0.0;
                    modes_.push_back(mode);
                }
            }
        }
    }
    
    // 解析Huang-Rhys因子
    void parseHuangRhys(std::ifstream& file) {
        std::string line;
        std::regex hr_pattern(R"(Mode\s+num\.\s+(\d+)\s+-\s+Factor:\s+([0-9.]+)D([+-]?\d+))");
        bool in_hr_section = false;
        int matched_count = 0;
        
        file.clear();
        file.seekg(0);
        
        while (std::getline(file, line)) {
            if (line.find("Huang-Rhys") != std::string::npos) {
                in_hr_section = true;
                continue;
            }
            
            if (in_hr_section) {
                std::smatch match;
                if (std::regex_search(line, match, hr_pattern)) {
                    int mode_num = std::stoi(match[1].str());
                    double mantissa = std::stod(match[2].str());
                    int exponent = std::stoi(match[3].str());
                    double huang_rhys = mantissa * std::pow(10.0, exponent);
                    
                    matched_count++;
                    
                    if (mode_num > 0 && mode_num <= static_cast<int>(modes_.size())) {
                        modes_[mode_num - 1].huang_rhys = huang_rhys;
                    }
                }
                
                if (matched_count > 10 && line.find("===") == std::string::npos && 
                    line.find("Mode") == std::string::npos && line.length() < 10) {
                    break;
                }
            }
        }
    }
    
    // 计算重组能贡献
    void calculateLambda() {
        for (auto& mode : modes_) {
            mode.lambda = mode.huang_rhys * mode.frequency;
        }
    }
    
    // 解析内坐标贡献
    CoordinateContribution parseCoordinateLine(const std::string& line) {
        CoordinateContribution contrib;
        
        if (line.find('!') == std::string::npos) return contrib;
        if (line.find("----") != std::string::npos) return contrib;
        if (line.find("Name") != std::string::npos) return contrib;
        
        size_t firstExclaim = line.find('!');
        size_t lastExclaim = line.rfind('!');
        
        if (firstExclaim == std::string::npos || lastExclaim == std::string::npos ||
            firstExclaim == lastExclaim) {
            return contrib;
        }
        
        std::string dataStr = line.substr(firstExclaim + 1, lastExclaim - firstExclaim - 1);
        std::vector<std::string> tokens = split(dataStr);
        
        if (tokens.size() < 4) return contrib;
        
        try {
            contrib.name = tokens[0];
            contrib.definition = tokens[1];
            
            std::vector<double> numbers;
            for (size_t i = 2; i < tokens.size(); i++) {
                if (isNumber(tokens[i])) {
                    numbers.push_back(std::stod(tokens[i]));
                }
            }
            
            if (numbers.size() >= 2) {
                contrib.weight = numbers[numbers.size() - 1];  // 最后一个数字是权重
            } else {
                contrib.name.clear();
            }
            
        } catch (const std::exception& e) {
            contrib.name.clear();
        }
        
        return contrib;
    }
    
    // 解析所有内坐标数据
    void parseCoordinateContributions(std::ifstream& file) {
        std::string line;
        
        file.clear();
        file.seekg(0);
        
        ModeCoordinateData* currentMode = nullptr;
        bool in_mode_data = false;
        
        while (std::getline(file, line)) {
            // 检测新的Normal Mode
            if (line.find("Normal Mode") != std::string::npos && 
                line.find("!") != std::string::npos) {
                size_t pos = line.find("Normal Mode");
                std::string afterMode = line.substr(pos + 11);
                std::istringstream iss(afterMode);
                int modeNum;
                if (iss >> modeNum) {
                    coordinate_data_.push_back(ModeCoordinateData());
                    currentMode = &coordinate_data_.back();
                    currentMode->mode_number = modeNum;
                    in_mode_data = true;
                }
                continue;
            }
            
            // 检测表头行
            if (currentMode != nullptr && in_mode_data && 
                line.find("Name") != std::string::npos && 
                line.find("Definition") != std::string::npos) {
                if (line.find("Relative Weight") == std::string::npos) {
                    in_mode_data = false;
                    currentMode = nullptr;
                }
                continue;
            }
            
            // 检测结束标记
            if (line.find("****") != std::string::npos || 
                line.find("GradGrad") != std::string::npos ||
                line.find("Optimized Parameters") != std::string::npos) {
                in_mode_data = false;
                currentMode = nullptr;
                continue;
            }
            
            // 解析数据行
            if (currentMode != nullptr && in_mode_data) {
                CoordinateContribution contrib = parseCoordinateLine(line);
                if (!contrib.name.empty()) {
                    currentMode->contributions.push_back(contrib);
                }
            }
        }
    }
    
public:
    ReorganizationEnergyAnalyzer(const std::string& filename) 
        : filename_(filename) {}
    
    bool loadData() {
        std::ifstream file(filename_);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file " << filename_ << std::endl;
            return false;
        }
        
        std::cout << "Processing file: " << filename_ << "\n\n";
        
        std::cout << "Parsing frequencies..." << std::endl;
        parseFrequencies(file);
        std::cout << "  Found " << modes_.size() << " vibrational modes\n";
        
        std::cout << "Parsing Huang-Rhys factors..." << std::endl;
        parseHuangRhys(file);
        int hr_count = 0;
        for (const auto& mode : modes_) {
            if (mode.huang_rhys > 0.0) hr_count++;
        }
        std::cout << "  Found " << hr_count << " Huang-Rhys factors\n";
        
        std::cout << "Calculating reorganization energies..." << std::endl;
        calculateLambda();
        
        std::cout << "Parsing internal coordinate contributions..." << std::endl;
        parseCoordinateContributions(file);
        std::cout << "  Found coordinate data for " << coordinate_data_.size() << " modes\n";
        if (coordinate_data_.empty()) {
            std::cout << "  Warning: No internal coordinate contribution data detected.\n";
        }
        
        file.close();
        return !modes_.empty();
    }
    
    // 计算内坐标对重组能的贡献
    void calculateCoordinateReorgContributions() {
        // 计算总重组能
        double total_lambda = 0.0;
        for (const auto& mode : modes_) {
            total_lambda += mode.lambda;
        }
        
        if (total_lambda == 0.0) {
            std::cerr << "Error: Total reorganization energy is zero!\n";
            return;
        }
        
        std::cout << "\n" << std::string(100, '=') << "\n";
        std::cout << "Internal Coordinate Contribution to Reorganization Energy\n";
        std::cout << std::string(100, '=') << "\n\n";
        
        // 存储每个内坐标的贡献
        std::map<std::string, double> coord_contributions;
        std::map<std::string, std::string> coord_definitions;
        
        // 对每个模式，计算各内坐标的贡献
        for (const auto& coord_mode : coordinate_data_) {
            int mode_num = coord_mode.mode_number;
            
            // 找到对应的振动模式
            if (mode_num < 1 || mode_num > static_cast<int>(modes_.size())) {
                continue;
            }
            
            const auto& mode = modes_[mode_num - 1];
            double lambda_i = mode.lambda;
            
            // 对该模式中的每个内坐标
            for (const auto& contrib : coord_mode.contributions) {
                double ci = contrib.weight / 100.0;  // 转换为小数
                double contribution = ci * lambda_i / total_lambda * 100.0;  // 转换为百分比
                
                coord_contributions[contrib.name] += contribution;
                coord_definitions[contrib.name] = contrib.definition;
            }
        }
        
        // 按坐标类型分组
        std::map<std::string, std::vector<std::pair<std::string, double>>> type_contributions;
        for (const auto& pair : coord_contributions) {
            std::string type = getCoordinateType(pair.first);
            type_contributions[type].push_back(pair);
        }
        
        // 打印详细贡献（按类型）
        std::cout << "Detailed Contributions by Coordinate Type:\n";
        std::cout << std::string(100, '-') << "\n";
        
        // 标准类型顺序：R, A, D
        std::vector<std::string> standard_types = {"R", "A", "D"};
        
        for (const auto& type : standard_types) {
            if (type_contributions.find(type) == type_contributions.end()) {
                continue;
            }
            
            auto& coords = type_contributions[type];
            
            // 按贡献大小排序
            std::sort(coords.begin(), coords.end(),
                     [](const auto& a, const auto& b) { return a.second > b.second; });
            
            std::cout << "\n" << type << " (Bond Lengths/Angles/Dihedrals):\n";
            std::cout << std::setw(8) << "Coord"
                      << std::setw(20) << "Definition "
                      << std::setw(15) << "Contribution(%)"
                      << "\n";
            std::cout << std::string(45, '-') << "\n";
            
            for (const auto& coord : coords) {
                std::cout << std::setw(8) << coord.first
                          << std::setw(20) << coord_definitions[coord.first]
                          << std::setw(15) << std::fixed << std::setprecision(4) 
                          << coord.second << "\n";
            }
        }
        
        // 计算并打印各类型的总贡献
        std::cout << "\n" << std::string(100, '=') << "\n";
        std::cout << "Summary: Total Contribution by Coordinate Type\n";
        std::cout << std::string(100, '=') << "\n\n";
        
        std::map<std::string, double> type_totals;
        for (const auto& type_pair : type_contributions) {
            double total = 0.0;
            for (const auto& coord : type_pair.second) {
                total += coord.second;
            }
            type_totals[type_pair.first] = total;
        }
        
        // 打印汇总结果
        std::cout << std::setw(20) << "Coordinate Type"
                  << std::setw(15) << "Total %"
                  << std::setw(25) << "Description"
                  << "\n";
        std::cout << std::string(60, '-') << "\n";
        
        // 标准输出顺序
        if (type_totals.find("R") != type_totals.end()) {
            std::cout << std::setw(20) << "R (Bond Lengths)"
                      << std::setw(15) << std::fixed << std::setprecision(2) 
                      << type_totals["R"] << "%"
                      << "\n";
        }
        
        if (type_totals.find("A") != type_totals.end()) {
            std::cout << std::setw(20) << "A (Bond Angles)"
                      << std::setw(15) << std::fixed << std::setprecision(2) 
                      << type_totals["A"] << "%"
                      << "\n";
        }
        
        if (type_totals.find("D") != type_totals.end()) {
            std::cout << std::setw(20) << "D (Dihedral Angles)"
                      << std::setw(15) << std::fixed << std::setprecision(2) 
                      << type_totals["D"] << "%"
                      << "\n";
        }
        
        // 打印其他类型
        for (const auto& pair : type_totals) {
            if (pair.first != "R" && pair.first != "A" && pair.first != "D") {
                std::cout << std::setw(20) << (pair.first + " (Other)")
                          << std::setw(15) << std::fixed << std::setprecision(2) 
                          << pair.second << "%"
                          << "\n";
            }
        }
        
        std::cout << std::string(60, '=') << "\n";
        
        // 验证总和
        double grand_total = 0.0;
        for (const auto& pair : type_totals) {
            grand_total += pair.second;
        }
        std::cout << "\nTotal (should be ~100%): " << std::fixed << std::setprecision(2) 
                  << grand_total << "%\n";
        
        // 保存到文件
        saveResults(coord_contributions, coord_definitions, type_totals);
        saveLowFrequencyCoordinateFile();
    }
    
    void saveResults(const std::map<std::string, double>& coord_contributions,
                    const std::map<std::string, std::string>& coord_definitions,
                    const std::map<std::string, double>& type_totals) {
        
        std::string output_file = "intmode.dat";
        std::ofstream out(output_file);
        if (!out.is_open()) {
            std::cerr << "Error: Cannot create output file " << output_file << std::endl;
            return;
        }
        
        out << "# Internal Coordinate Contributions to Reorganization Energy\n";
        out << "# Input file: " << filename_ << "\n";
        out << "#\n";
        out << "# Coordinate  Type  Definition  Contribution(%)\n";
        
        // 按类型分组输出
        std::vector<std::string> types = {"R", "A", "D"};
        for (const auto& type : types) {
            out << "#\n# " << type << " coordinates:\n";

            std::vector<std::pair<std::string, double>> type_entries;
            for (const auto& pair : coord_contributions) {
                if (getCoordinateType(pair.first) == type) {
                    type_entries.push_back(pair);
                }
            }

            std::sort(type_entries.begin(), type_entries.end(),
                      [](const auto& a, const auto& b) {
                          return a.second > b.second;
                      });

            for (const auto& entry : type_entries) {
                out << std::setw(10) << entry.first << "  "
                    << std::setw(4) << type << "  "
                    << std::setw(20) << coord_definitions.at(entry.first) << "  "
                    << std::setw(12) << std::fixed << std::setprecision(6) 
                    << entry.second << "\n";
            }
        }
        
        out << "#\n# Summary:\n";
        for (const auto& pair : type_totals) {
            out << "# " << pair.first << "_total: " 
                << std::fixed << std::setprecision(4) << pair.second << " %\n";
        }
        
        out.close();
        std::cout << "\nInternal coordinate contributions saved to: " << output_file << "\n";
    }
    
    void saveLambdaResults(const std::string& output_file = "lambda.dat") {
        std::ofstream out(output_file);
        if (!out.is_open()) {
            std::cerr << "Error: Cannot create output file " << output_file << std::endl;
            return;
        }
        
        out << "# Reorganization Energy Analysis\n";
        out << "# Input file: " << filename_ << "\n";
        out << "#\n";
        out << "# Mode  Frequency(cm-1)  Huang-Rhys  Lambda(cm-1)  %Contribution\n";
        
        double total_lambda = 0.0;
        for (const auto& mode : modes_) {
            total_lambda += mode.lambda;
        }
        
        for (const auto& mode : modes_) {
            double percent = (total_lambda > 0.0) ? (mode.lambda / total_lambda * 100.0) : 0.0;
            out << std::setw(6) << mode.mode_number << "  "
                << std::setw(14) << std::fixed << std::setprecision(6) << mode.frequency << "  "
                << std::setw(14) << std::scientific << std::setprecision(6) << mode.huang_rhys << "  "
                << std::setw(14) << std::fixed << std::setprecision(6) << mode.lambda << "  "
                << std::setw(12) << std::setprecision(4) << percent << "\n";
        }
        
        double low_freq_lambda = 0.0;
        int low_freq_count = 0;
        for (const auto& mode : modes_) {
            if (mode.frequency < 200.0) {
                low_freq_lambda += mode.lambda;
                low_freq_count++;
            }
        }
        
        out << "#\n";
        out << "# Summary:\n";
        out << "# Total lambda: " << std::fixed << std::setprecision(4) << total_lambda << " cm-1\n";
        out << "# Low freq (<200 cm-1) lambda: " << low_freq_lambda << " cm-1\n";
        out << "# Low freq percentage: " 
            << std::setprecision(2) 
            << (total_lambda > 0.0 ? low_freq_lambda / total_lambda * 100.0 : 0.0) << " %\n";
        out << "# Low freq mode count: " << low_freq_count << "\n";
        
        out.close();
        std::cout << "Reorganization energy data saved to: " << output_file << "\n";
    }
    
    void saveLowFrequencyCoordinateFile(double threshold = 208.5,
                                        const std::string& output_file = "int_lowfreq.dat") {
        if (coordinate_data_.empty()) {
            std::cout << "Skipping low-frequency coordinate file: no coordinate data.\n";
            return;
        }

        std::vector<bool> is_low(modes_.size(), false);
        double total_lambda = 0.0;
        int low_freq_mode_count = 0;

        for (size_t i = 0; i < modes_.size(); ++i) {
            if (modes_[i].frequency <= threshold) {
                is_low[i] = true;
                total_lambda += modes_[i].lambda;
                if (modes_[i].lambda > 0.0) {
                    low_freq_mode_count++;
                }
            }
        }

        if (total_lambda == 0.0) {
            std::cout << "Skipping low-frequency coordinate file: total lambda below threshold is zero.\n";
            return;
        }

        std::map<std::string, double> coord_contributions;
        std::map<std::string, std::string> coord_definitions;

        for (const auto& coord_mode : coordinate_data_) {
            int mode_index = coord_mode.mode_number - 1;
            if (mode_index < 0 || mode_index >= static_cast<int>(modes_.size())) {
                continue;
            }
            if (!is_low[mode_index]) {
                continue;
            }

            double lambda_i = modes_[mode_index].lambda;
            if (lambda_i == 0.0) {
                continue;
            }

            for (const auto& contrib : coord_mode.contributions) {
                double ci = contrib.weight / 100.0;
                double contribution = ci * lambda_i / total_lambda * 100.0;
                coord_contributions[contrib.name] += contribution;
                coord_definitions[contrib.name] = contrib.definition;
            }
        }

        if (coord_contributions.empty()) {
            std::cout << "Skipping low-frequency coordinate file: no contributions found under threshold.\n";
            return;
        }

        std::vector<std::pair<std::string, double>> sorted_contribs(coord_contributions.begin(),
                                                                    coord_contributions.end());
        std::sort(sorted_contribs.begin(), sorted_contribs.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });

        std::ofstream out(output_file);
        if (!out.is_open()) {
            std::cerr << "Error: Cannot create output file " << output_file << std::endl;
            return;
        }

        out << "# Low-frequency Internal Coordinate Contributions\n";
        out << "# Input file: " << filename_ << "\n";
        out << "# Threshold: <= " << threshold << " cm-1\n";
        out << "# Modes included: " << low_freq_mode_count << "\n";
        out << "#\n";
        out << "# Coordinate  Definition  Contribution(%)\n";

        for (const auto& entry : sorted_contribs) {
            out << std::setw(10) << entry.first << "  "
                << std::setw(20) << coord_definitions[entry.first] << "  "
                << std::setw(12) << std::fixed << std::setprecision(6) << entry.second << "\n";
        }

        out.close();
        std::cout << "Low-frequency coordinate contributions saved to: " << output_file << "\n";
    }

    bool hasCoordinateData() const {
        return !coordinate_data_.empty();
    }
    
    void printReorgEnergySummary() {
        double total_lambda = 0.0;
        for (const auto& mode : modes_) {
            total_lambda += mode.lambda;
        }
        

        std::cout << "\nTotal reorganization energy: " 
                  << std::fixed << std::setprecision(4) << total_lambda << " cm^-1\n";
    }
};

// ==================== 主函数 ====================

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <gaussian_log_file>\n";
        std::cerr << "\nThis program calculates internal coordinate contributions\n";
        std::cerr << "to reorganization energy from Gaussian output.\n";
        std::cerr << "\nThe input file must contain:\n";
        std::cerr << "  1. Frequency calculation results\n";
        std::cerr << "  2. Huang-Rhys factors\n";
        std::cerr << "  3. Internal coordinate analysis (freq=intmodes)\n";
        std::cerr << "\nExample: " << argv[0] << " molecule.log\n";
        return 1;
    }
    
    std::string input_file = argv[1];
    
    ReorganizationEnergyAnalyzer analyzer(input_file);
    
    if (!analyzer.loadData()) {
        std::cerr << "Error: Failed to load data from file\n";
        std::cerr << "Make sure the file contains frequency, Huang-Rhys, and intmodes data\n";
        return 1;
    }
    
    analyzer.printReorgEnergySummary();
    if (analyzer.hasCoordinateData()) {
        analyzer.calculateCoordinateReorgContributions();
    } else {
        std::cout << "Skipping internal coordinate contribution analysis: "
                  << "no intmode results were found in the input.\n";
    }
    analyzer.saveLambdaResults();
    
    std::cout << "\n" << std::string(100, '=') << "\n";
    std::cout << "Analysis complete!\n";
    std::cout << "Author: Bane Dysta\n";
    std::cout << "Feedback: https://github.com/bane-dysta/intmode\n";
    std::cout << "          http://bbs.keinsci.com/forum.php?mod=viewthread&tid=57134&fromuid=63020\n";
    std::cout << std::string(100, '=') << "\n\n";
    
    return 0;
}