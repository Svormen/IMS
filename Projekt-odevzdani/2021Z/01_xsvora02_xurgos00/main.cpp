// IMS project
// Authors: Peter Urgos (xurgos00), Slavomir Svorada (xsvora02)
// 12/2021

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

// Implementation from https://c.mql5.com/3/133/Tutorial_on_Burg_smethod_algorithm_recursion.pdf
// Details about the Burg's AR approximation mehtod: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1163429
#include "arburg.h"

#define ROWS 13
#define COLS 8

//#define DEBUG

void printHelp();
[[noreturn]] void invalidArguments();
void predictCPI(unsigned variant);
void doExperiment();
void parseData(std::string filename, unsigned startRow, unsigned startCol, unsigned endRow, unsigned endCol, std::vector<std::vector<double>> &v);
void writeData(unsigned row, unsigned col, double data, std::vector<std::vector<double>> &v);
std::vector<std::string> readFile(std::string filename);
bool noArgument = false;
bool experiment1 = false;
bool experiment2 = false;

int main(int argc, char *argv[])
{
    if (argc == 1)
    {
        noArgument = true;
        doExperiment();
    }
    else if (argc == 2)
    {
        if (argv[1] == std::string("-h") || argv[1] == std::string("--help"))
        {
            printHelp();
            return 0;
        }
        else if (argv[1] == std::string("-p_1") || argv[1] == std::string("--predict-1"))
        {
            if (argc != 2)
                invalidArguments();

            predictCPI(1);
            return 0;
        }
        else if (argv[1] == std::string("-p_2") || argv[1] == std::string("--predict-2"))
        {
            if (argc != 2)
                invalidArguments();

            predictCPI(2);
            return 0;
        }
        else if (argv[1] == std::string("-p_3") || argv[1] == std::string("--predict-3"))
        {
            if (argc != 2)
                invalidArguments();

            predictCPI(3);
            return 0;
        }
        else if (argv[1] == std::string("-p_4") || argv[1] == std::string("--predict-4"))
        {
            if (argc != 2)
                invalidArguments();

            predictCPI(4);
            return 0;
        }
        else if ((argv[1] == std::string("-e_1")) || (argv[1] == std::string("-e_2")))
        {
            if (argv[1] == std::string("-e_1")) {
                experiment1 = true;
            } else {
                experiment2 = true;
            }
            doExperiment();
            return 0;
        }

        invalidArguments();
    }
    else
    {
        invalidArguments();
    }

    return 0;
}

void printHelp()
{
    std::cout << "Usage: ./main [options]\n";
    std::cout << "Options:\n";
    std::cout << "    -h, --help               Print this message\n";
    std::cout << "    -e_1                     Adjust most significant factor of CPI                | experiment 1\n";
    std::cout << "    -e_2                     Adjust least significant factor of CPI               | experiment 2\n";
    std::cout << "    -p_1, --predict-1        Predict next 20 values of CPI inflation using AR(18) | experiment 3\n";
    std::cout << "    -p_2, --predict-2        Predict next 20 values of CPI inflation using AR(15) | experiment 3\n";
    std::cout << "    -p_3, --predict-3        Predict next 20 values of CPI inflation using AR(10) | experiment 3\n";
    std::cout << "    -p_4, --predict-4        Predict next 20 values of CPI inflation using AR(5)  | experiment 3\n";
}

[[noreturn]] void invalidArguments()
{
    std::cerr << "Invalid arguments. For help see './main --help'\n";
    exit(1);
}

void predictCPI(unsigned variant)
{
    // AR degree must be less than or equal to 'data.length() - 1'
    // Degree is the same as lag, which means how many previous values are taken into account in the AR model
    unsigned degree;
    if (variant == 2) degree = 15;
    else if (variant == 3) degree = 10;
    else if (variant == 4) degree = 5;
    else degree = 18;

    unsigned predict_n = 20;

    // Yearly CPI % increase from 2002 to 2020 (http://datacube.statistics.sk/#!/view/sk/VBD_SLOVSTAT/sp2043rs/v_sp2043rs_00_00_00_sk)
    std::vector<double> dataCPI = { 3.30, 8.50, 7.50, 2.70, 4.50, 2.80, 4.60, 1.60, 1, 3.90, 3.60, 1.40, -0.1, -0.3, -0.5, 1.3, 2.5, 2.7, 1.9 };
    std::vector<double> coeffs(degree);

    for (unsigned pred_count = 0; pred_count < predict_n; pred_count++)
    {
        BurgAlgorithm(coeffs, dataCPI);

        double prediction = 0;
        for (unsigned i = 0; i < degree; i++)
        {
            //std::cout << -coeffs[i] << " | ";
            prediction += -coeffs[i] * dataCPI[dataCPI.size() - 1 - i];
        }
        //prediction += ((double) random()) / RAND_MAX * 2.0 - 1.0;

        //std::cout << "Prediction: " << prediction << "\n";
        std::cout << prediction << ", ";
        dataCPI.push_back(prediction);
    }

    std::cout << "\n";
}

void doExperiment()
{
    std::vector<std::vector<double>> dataCPI(ROWS, std::vector<double>(COLS));
    if (noArgument) {
        parseData("data_0.csv", 1, 1, dataCPI.size(), dataCPI[0].size(), dataCPI);
    } else if (experiment1) {
        parseData("data_ex1.csv", 1, 1, dataCPI.size(), dataCPI[0].size(), dataCPI);
    } else {
        parseData("data_ex2.csv", 1, 1, dataCPI.size(), dataCPI[0].size(), dataCPI);
    }
    for (unsigned i = 0; i < dataCPI.size(); i++)
        for (unsigned j = 0; j < dataCPI[0].size(); j++)
            dataCPI[i][j] -= 100;

    std::vector<std::vector<double>> dataSpotrebnyKos(ROWS - 1, std::vector<double>(COLS));
    parseData("data_1.csv", 1, 1, dataSpotrebnyKos.size(), dataSpotrebnyKos[0].size(), dataSpotrebnyKos);
    for (unsigned i = 0; i < dataSpotrebnyKos.size(); i++)
        for (unsigned j = 0; j < dataSpotrebnyKos[0].size(); j++)
            dataSpotrebnyKos[i][j] /= 100;

#ifdef DEBUG
    std::cout << "\n\t<<< BEGIN DEBUG >>>\n";
    std::cout << "dataCPI:\n";
    for (int i = 0; i < ROWS; i++)
    {
        std::cout << "  ";
        for (int j = 0; j < COLS; j++)
        {
            std::cout << dataCPI[i][j] << " - ";
        }
        std::cout << "\n";
    }
    std::cout << "\ndataSpotrebnyKos:\n";
    for (int i = 0; i < ROWS - 1; i++)
    {
        std::cout << "  ";
        for (int j = 0; j < COLS; j++)
        {
            std::cout << dataSpotrebnyKos[i][j] << " - ";
        }
        std::cout << "\n";
    }
    std::cout << "\t<<< END DEBUG >>>\n\n";
#endif

    std::vector<double> yearAverages(COLS);
    // Iterate years
    for (int year = 0; year < COLS; year++)
    {
        // Iterate categories
        for (int j = 1; j < ROWS; j++)
        {
            yearAverages[year] += dataCPI[j][year];
        }
        // Average
        yearAverages[year] /= ROWS - 1;
    }

    std::vector<double> predikcie(COLS);
    for (int year = 0; year < COLS; year++)
    {
        for (int j = 1; j < ROWS; j++)
        {
            predikcie[year] += dataCPI[j][year] * dataSpotrebnyKos[j - 1][year];
        }

        predikcie[year] /= ROWS - 1;
    }

    if (experiment1 || experiment2)
    {
        std::cout << "Predikcia: " << predikcie[0] <<"\t| Uhrn: " << dataCPI[0][0] << " \t| Priemer: " << yearAverages[0] << "\n";
    }
    else
    {
        for (int year = 0; year < COLS; year++)
        {
            std::cout << "Predikcia: " << predikcie[year] <<"\t| Uhrn: " << dataCPI[0][year] << " \t| Priemer: " << yearAverages[year] << "\n";
        }
    }
}

void parseData(std::string filename,
               unsigned startRow, unsigned startCol,
               unsigned endRow, unsigned endCol,
               std::vector<std::vector<double>> &v)
{
    auto lines = readFile(filename);

    bool inString = false;
    unsigned row = 0;

    for (const auto &line : lines)
    {
        unsigned col = 0;
        std::string temp = "";

        for (unsigned long j = 0; j < line.length(); j++)
        {
            char c = line[j];

            if (c == '"')
            {
                if (inString)
                {
                    inString = false;
                }
                else
                {
                    inString = true;
                }
            }
            else if (c == ',')
            {
                if (!inString)
                {
                    if (row >= startRow && col >= startCol &&
                        row <= endRow   && col <= endCol)
                    {
                        writeData(row - 1, col - 1, std::stod(temp), v);
                    }
                    temp = "";
                    col++;
                }
            }
            else
            {
                if (!inString)
                {
                    temp += c;
                }
            }
        }

        if (row >= startRow && col >= startCol &&
            row <= endRow   && col <= endCol)
        {
            writeData(row - 1, col - 1, std::stod(temp), v);
        }

        row++;
    }
}

void writeData(unsigned row, unsigned col, double data, std::vector<std::vector<double>> &v)
{
    if (row >= v.size() || col > v[row].size())
    {
        std::cout   << "Invalid dimensions '" << row << "x" << col << "' for an 2d vector of size '"
                    << v.size() << "x" << v[0].size() << "'\n";

        return;
    }
    v[row][col] = data;
}

std::vector<std::string> readFile(std::string filename)
{
    std::vector<std::string> lines;
    std::string line;

    std::ifstream input_file(filename);
    if (!input_file.is_open())
    {
        std::cerr << "Could not open the file - '" << filename << "'" << std::endl;
        return lines;
    }

    while (getline(input_file, line)){
        lines.push_back(line);
    }

    input_file.close();

    return lines;
}
