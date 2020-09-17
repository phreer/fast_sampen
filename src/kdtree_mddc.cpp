#include <iostream>
#include <string>
#include <time.h>

#include "sample_entropy_calculator_direct.h"
#include "sample_entropy_calculator_kd.h"
#include "utils.h"

using namespace kdtree_mddc;
using std::cout;
using std::endl;
using std::string;

#define SAMPLE_ENTROPY(T,m,fast_direct,direct) SampleEntropy<(T),(m)>(\
    arg.filename,arg.input_type,arg.r,arg.data_length,(fast_direct),(direct),\
    arg.output_level)

char *usage =\
"Usage: kdtree_mddc -input <INPUT> -input-type {simple, multirecord}\\\n"\
"                   -r <THRESHOLD> -m <TEMPLATE_LENGTH>\\\n"\
"                   -n <N> [-output-level {1,2,3}]\n\n"\
"-input <INPUT>          The file name of the input file\n"\
"-input-type <TYPE>      The format of the input file. Should be either simple\n"\
"                        or multirecord. If set to simple, then each line of the\n"\
"                        input file contains exactly one column; if set to\n"\
"                        multirecord, then each line contains <NUM_RECORD> + 1\n"\
"                        columns, of which the first indicates the line number and\n"\
"                        and the remaining columns are instances of the records.\n"\
"                        The defualt value is simple.\n"\
"-r <R>                  The threshold argument in sample entropy.\n"\
"-m <M>                  The template length argument of sample entropy. Note\n"\
"                        that this program only supports 2 <= m <= 10.\n"
"-n <N>                  If the length of the signal specified by <FILENAME> is\n"
"                        greater than <N>, then it would be truncated to be of\n"
"                        length <N>. If <N> is 0, then the the original length\n"
"                        is employed. The default value is 0.\n";

template<typename T>
void PrintSampenSetting(unsigned line_offset, unsigned n, T r, unsigned K, 
                        std::string filename, double var);

void ParseArgument(int argc, char *argv[]);

struct Argument
{
    unsigned template_length;
    unsigned line_offset;
    string filename;
    string input_type;
    string type;
    unsigned data_length;
    double r;
    OutputLevel output_level;
} arg;

template<typename T, unsigned K>
void SampleEntropy(const string &filename, const string &input_type, double r, 
                   unsigned line_offset, unsigned n, bool fast_direct, 
                   bool direct, OutputLevel output_level);

int main(int argc, char *argv[])
{
    //std::string filename = R"(C:\Users\95774\OneDrive\Codes\workspace\sampen\data\power_two\ecg_114157_2p17.txt)";
    //std::string filename = R"(C:\Users\95774\Desktop\ex8\data\Health\Filt-time-1514.01.rr.txt)";
    //std::string filename = R"(C:\Users\95774\Desktop\ex8\data\AF\fa002.rr.txt)";
    //vector<double> data({ 6, 3, 7, 2, 1, 9, 4, 6, 3, 7, 2, 1, 5, 0 });
#ifdef DEBUG
    std::cout << "This is a debug version. " << std::endl;
#endif
    ParseArgument(argc, argv);
    if (arg.template_length < 2 || arg.template_length > 10)
    {
        cerr << "The template length " << arg.template_length;
        cerr << " should be integers between 2 and 10 .\n";
        exit(-1);
    }

    if (arg.output_level) std::cout << "Data type: " << arg.type << endl;
    bool fast_direct = true;
    bool direct = true;

    if (arg.type == "double")
    {
        using Type = double;
        switch (arg.template_length)
        {
        case 2:
            SampleEntropy<Type, 2>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 3:
            SampleEntropy<Type, 3>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 4:
            SampleEntropy<Type, 4>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 5:
            SampleEntropy<Type, 5>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 6:
            SampleEntropy<Type, 6>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 7:
            SampleEntropy<Type, 7>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 8:
            SampleEntropy<Type, 8>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 9:
            SampleEntropy<Type, 9>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 10: 
            SampleEntropy<Type, 10>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        default:
            cerr << "Invalid argument: -input-type " << arg.template_length << ". \n";
            cerr << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;
            exit(-1);
        }
    } 
    else if (arg.type == "int")
    {
        using Type = int;
        switch (arg.template_length)
        {
        case 2:
            SampleEntropy<Type, 2>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 3:
            SampleEntropy<Type, 3>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 4:
            SampleEntropy<Type, 4>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 5:
            SampleEntropy<Type, 5>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 6:
            SampleEntropy<Type, 6>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 7:
            SampleEntropy<Type, 7>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 8:
            SampleEntropy<Type, 8>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 9:
            SampleEntropy<Type, 9>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        case 10: 
            SampleEntropy<Type, 10>(arg.filename, arg.input_type, arg.r, 
                                   arg.line_offset, arg.data_length, 
                                   fast_direct, direct, arg.output_level); 
            break;
        default:
            cerr << "Invalid argument: -input-type " << arg.template_length << ". \n";
            cerr << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;
            exit(-1);
        }
    }
    else
    {
        cerr << "Invalid argument: -type " << arg.type << ".\n";
        cerr << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;
        exit(-1);
    }

    return 0;
}

template<typename T>
void PrintSampenSetting(unsigned line_offset, unsigned n, T r, unsigned m, 
                        std::string filename, double var)
{
    std::cout << "Sample Entropy Computation Setting: " << std::endl;
    std::cout << "\tfilename: " << filename << std::endl;
    std::cout << "\tline_offset: " << line_offset << std::endl;
    std::cout << "\tdata length: " << n << std::endl;
    std::cout << "\ttemplate length: " << m << std::endl;
    std::cout << "\tthreshold: " << r << std::endl;
    std::cout << "\tvariance: " << var << std::endl;
}

void ParseArgument(int argc, char *argv[])
{
    ArgumentParser parser(argc, argv);
    arg.filename = parser.getArg("-input");
    if (arg.filename.size() == 0)
    {
        cerr << "Specify a filename with -input <INPUT>." << endl;
        cerr << usage;
        exit(-1);
    }

    arg.input_type = parser.getArg("-input-type");
    if (arg.input_type.size() == 0) arg.input_type = "simple";
    if (arg.input_type != "simple" && arg.input_type != "multirecord") 
    {
        cerr << "Invalid argument: input-type, should be either simple or "; 
        cerr << "multirecord. \n";
        exit(-1);
    }

    arg.type = parser.getArg("-type");
    if (arg.type.size() == 0) arg.type = "double";

    string r = parser.getArg("-r");
    if (r.size() == 0)
    {
        cerr << "Specify a threshold with -r <R>. " << endl;
        cerr << usage;
        exit(-1);
    }
    arg.r = std::stod(r);

    string template_length = parser.getArg("-m");
    if (template_length.size() == 0)
    {
        cerr << "Specify template length with -m <M> (1 <= M <= 10." << endl;
        cerr << usage;
        exit(-1);
    } else
    {
        arg.template_length = static_cast<unsigned>(std::stoi(template_length));
    }

    string data_length = parser.getArg("-n");
    if (data_length.size() == 0) arg.data_length = 0;
    else arg.data_length = static_cast<unsigned>(std::stoi(data_length));

    string line_offset = parser.getArg("--line-offset");
    if (line_offset.size() == 0) arg.line_offset = 0;
    else arg.line_offset = static_cast<unsigned>(std::stoi(line_offset));
    
    string output_level = parser.getArg("--output-level");
    if (output_level.size() == 0) arg.output_level = Info;
    else if (output_level == "1") arg.output_level = Silent;
    else if (output_level == "2") arg.output_level = Info;
    else if (output_level == "3") arg.output_level = Debug;
    else
    {
        cerr << "Invalid argument -output-level " << output_level;
        cerr << ", must be one of the following: 1, 2, 3. \n";
        cerr << usage;
        exit(-1);
    }
}


template<typename T, unsigned K>
void SampleEntropy(const string &filename, const string &input_type, double r, 
                   unsigned line_offset, unsigned n, bool fast_direct, 
                   bool direct, OutputLevel output_level)
{
    vector<T> data = ReadData<T>(arg.filename, input_type, n, line_offset);
    n = data.size();
    if (n <= K)
    {
        std::cerr << "Data length n " << n << " is too short (K = " << K; 
        std::cerr << "). \n";
        exit(-1);
    }

    double var = ComputeVarience(data);
    T r_scaled = static_cast<T>(sqrt(var) * r);

    if (r_scaled < 0)
    {
        cerr << "Invalid r [scaled]: " << r_scaled << ".\n";
        exit(-1);
    }

    PrintSampenSetting(line_offset, n, r_scaled, K, arg.filename, var);

    // Compute sample entropy. 
    double sampen = 0;
    clock_t t;

    SampleEntropyCalculatorMao<T, K> sec(output_level);
    t = clock();
    sampen = sec.ComputeSampleEntropy(data.cbegin(), data.cend(), r_scaled);
    t = clock() - t;
    cout << "kd tree (Mao): Sampen(" << n << ", " << K;
    cout << ", " << r << "): " << sampen;
    cout << ", time: " << static_cast<double>(t) / CLOCKS_PER_SEC << " seconds. \n";

    SampleEntropyCalculatorLiu<T, K> secl(output_level);
    t = clock();
    sampen = secl.ComputeSampleEntropy(data.cbegin(), data.cend(), r_scaled);
    t = clock() - t;
    cout << "kd tree (Liu): Sampen(" << n << ", " << K;
    cout << ", " << r << "): " << sampen;
    cout << ", time: " << static_cast<double>(t) / CLOCKS_PER_SEC << " seconds. \n";

    if (fast_direct)
    {
        SampleEntropyCalculatorFastDirect<T, K> secfd(output_level);
        t = clock();
        sampen = secfd.ComputeSampleEntropy(data.cbegin(), data.cend(), r_scaled);
        t = clock() - t;
        cout << "kd tree (Fast Direct): Sampen(" << n << ", " << K;
        cout << ", " << r << "): " << sampen;
        cout << ", time: " << static_cast<double>(t) / CLOCKS_PER_SEC << " seconds. \n";
    }

    if (direct)
    {
        SampleEntropyCalculatorDirect<T, K> secd(output_level);
        t = clock();
        sampen = secd.ComputeSampleEntropy(data.cbegin(), data.cend(), r_scaled);
        t = clock() - t;
        cout << "kd tree (Direct): Sampen(" << n << ", " << K;
        cout << ", " << r << "): " << sampen;
        cout << ", time: " << static_cast<double>(t) / CLOCKS_PER_SEC << " seconds. \n";
    }
}
