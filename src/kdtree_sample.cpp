#include <iostream>
#include <string>
#include <time.h>

#include "random_sampler.h"
#include "sample_entropy_calculator_direct.h"
#include "sample_entropy_calculator_kd.h"
#include "utils.h"

using namespace kdtree_mddc;
using std::cout;
using std::endl;
using std::string;
using std::vector; 

char usage[] =\
"Usage: kdtree_mddc --input <INPUT> --input-type {simple, multirecord}\\\n"
"                   -r <THRESHOLD> -m <TEMPLATE_LENGTH>\\\n"
"                   -n <N> [-output-level {1,2,3}]\\\n"
"                   --sample-size <SAMPLE_SIZE> --sample-num <SAMPLE_NUM>\n\n"
"--input <INPUT>         The file name of the input file.\n"
"--input-format <FORMAT> The format of the input file. Should be either simple\n"
"                        or multirecord. If set to simple, then each line of the\n"
"                        input file contains exactly one column; if set to\n"
"                        multirecord, then each line contains <NUM_RECORD> + 1\n"
"                        columns, of which the first indicates the line number\n"
"                        and the remaining columns are instances of the records.\n"
"                        The default value is simple.\n"
"--input-type <TYPE>     The data type of the input data, either int or float.\n"
"                        Default: double.\n"
"-r <R>                  The threshold argument in sample entropy.\n"
"-m <M>                  The template length argument of sample entropy. Note\n"
"                        that this program only supports 2 <= m <= 10.\n"
"-n <N>                  If the length of the signal specified by <FILENAME> is\n"
"                        greater than <N>, then it would be truncated to be of\n"
"                        length <N>. If <N> is 0, then the the original length\n"
"                        is employed. The default value is 0.\n"
"--sample-size <N0>      The number of points to sample.\n"
"--sample-num <N1>       The number of computations where the average is taken.\n"
"--random                If this option is enabled, the random seed will be set\n"
"                        randomly.\n"
"-q                      If this option is enabled, the quasi-Monte Carlo based\n"
"                        method is conducted.\n"
"--quasi-type <TYPE>     The type of the quasi-random sequence for sampling,\n"
"                        can be one of the following: sobol, halton,\n"
"                        reversehalton, niederreiter_2 or grid. Default: sobol.\n"
"-u                      If this option is enabled, the Monte Carlo based\n"
"                        using uniform distribution is conducted.\n"
"--output-level <LEVEL>  The amount of information printed. Should be one of\n"
"                        {0,1,2}. Level 0 is most silent while level 2 is for\n"
"                        debugging.\n";

template<typename T>
void PrintSampenSetting(unsigned line_offset, unsigned n, T r, unsigned K, 
                        std::string filename, double var);

void ParseArgument(int argc, char *argv[]);

struct Argument
{
    unsigned template_length;
    unsigned line_offset;
    string filename;
    string input_format; 
    string input_type;
    unsigned data_length;
    unsigned sample_size; 
    unsigned sample_num; 
    double r;
    OutputLevel output_level;
    bool fast_direct = true; 
    bool direct = false; 
    bool random_; 
    bool q, u;
    RandomType rtype; 
    void PrintArguments() const;
} arg;

template<typename T, unsigned K>
void SampleEntropy();

int main(int argc, char *argv[])
{
#ifdef DEBUG
    cout << "This is a debug version. " << std::endl;
#endif
    ParseArgument(argc, argv);

    if (arg.input_type == "double")
    {
        using Type = double;
        switch (arg.template_length)
        {
        case 2:
            SampleEntropy<Type, 2>(); 
            break;
        case 3:
            SampleEntropy<Type, 3>(); 
            break;
        case 4:
            SampleEntropy<Type, 4>(); 
            break;
        case 5:
            SampleEntropy<Type, 5>(); 
            break;
        case 6:
            SampleEntropy<Type, 6>(); 
            break;
        case 7:
            SampleEntropy<Type, 7>(); 
            break;
        case 8:
            SampleEntropy<Type, 8>(); 
            break;
        case 9:
            SampleEntropy<Type, 9>(); 
            break;
        case 10: 
            SampleEntropy<Type, 10>(); 
            break;
        default:
            cerr << "Invalid argument: -input-type " << arg.template_length << ". \n";
            cerr << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;
            exit(-1);
        }
    } 
    else if (arg.input_type == "int")
    {
        using Type = int;
        switch (arg.template_length)
        {
        case 2:
            SampleEntropy<Type, 2>(); 
            break;
        case 3:
            SampleEntropy<Type, 3>(); 
            break;
        case 4:
            SampleEntropy<Type, 4>(); 
            break;
        case 5:
            SampleEntropy<Type, 5>(); 
            break;
        case 6:
            SampleEntropy<Type, 6>(); 
            break;
        case 7:
            SampleEntropy<Type, 7>(); 
            break;
        case 8:
            SampleEntropy<Type, 8>(); 
            break;
        case 9:
            SampleEntropy<Type, 9>(); 
            break;
        case 10: 
            SampleEntropy<Type, 10>(); 
            break;
        default:
            cerr << "Invalid argument: -input-type " << arg.template_length << ". \n";
            cerr << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;
            exit(-1);
        }
    }
    else
    {
        cerr << "Invalid argument: -type " << arg.input_type << ".\n";
        cerr << "File: " << __FILE__ << ", Line: " << __LINE__ << std::endl;
        exit(-1);
    }

    return 0;
}

void Argument::PrintArguments() const 
{
    std::cout << "Sample Entropy Computation Setting: " << std::endl;
    std::cout << "\tfilename: " << arg.filename << std::endl;
    std::cout << "\tinput type: " << arg.input_type << std::endl; 
    std::cout << "\tline offset: " << arg.line_offset << std::endl;
    std::cout << "\tdata length: " << arg.data_length << std::endl;
    std::cout << "\ttemplate length: " << arg.template_length << std::endl;
    std::cout << "\tthreshold: " << arg.r << std::endl;
    std::cout << "\trandom: " << arg.random_ << std::endl; 
    std::cout << "\tquasi type: " << random_type_names[arg.rtype] << std::endl;
    std::cout << "\tsample num: " << arg.sample_num << std::endl; 
    std::cout << "\tsample size: " << arg.sample_size << std::endl; 
}

void ParseArgument(int argc, char *argv[])
{
    ArgumentParser parser(argc, argv);
    long result_long; 
    arg.filename = parser.getArg("--input");
    if (arg.filename.size() == 0)
    {
        cerr << "Specify a filename with --input <INPUT>." << endl;
        cerr << usage;
        exit(-1);
    }

    arg.input_format = parser.getArg("--input-format");
    if (arg.input_format.size() == 0) arg.input_format = "simple";
    if (arg.input_format != "simple" && arg.input_format != "multirecord") 
    {
        cerr << "Invalid argument: --input-format " << arg.input_format; 
        cerr << ", should be either simple or multirecord. \n";
        exit(-1);
    }

    arg.input_type = parser.getArg("--input-type");
    if (arg.input_type.size() == 0) arg.input_type = "double";

    arg.r = parser.getArgDouble("-r", -1.); 
    if (arg.r < 0)
    {
        cerr << "Specify a positive threshold with -r <R>. " << endl;
        cerr << usage;
        exit(-1);
    }

    string template_length = parser.getArg("-m");
    if (template_length.size() == 0)
    {
        cerr << "Specify template length with -m <M> (1 <= M <= 10." << endl;
        cerr << usage;
        exit(-1);
    } else
    {
        arg.template_length = static_cast<unsigned>(std::stoi(template_length));
        if (arg.template_length < 1 || arg.template_length > 10) 
        {
            cerr << "Specify template length with -m <M> (1 <= M <= 10." << endl;
            exit(-1); 
        }
    }

    arg.data_length = static_cast<unsigned>(parser.getArgLong("-n", 0));
    arg.line_offset = static_cast<unsigned>(parser.getArgLong("--line-offset", 0));
    
    string output_level = parser.getArg("--output-level");
    if (output_level.size() == 0) arg.output_level = Info;
    else if (output_level == "0") arg.output_level = Silent;
    else if (output_level == "1") arg.output_level = Info;
    else if (output_level == "2") arg.output_level = Debug;
    else
    {
        cerr << "Invalid argument --output-level " << output_level;
        cerr << ", must be one of the following: {1, 2, 3}. \n";
        cerr << usage;
        exit(-1);
    }

    arg.q = parser.isOption("-q"); 
    arg.u = parser.isOption("-u"); 

    if (arg.q || arg.u) 
    {
        arg.random_ = parser.isOption("--random"); 
        if (arg.q) 
        {
            std::string rtype = parser.getArg("--quasi-type"); 
            if (rtype.size() == 0 || rtype == "sobol")
                arg.rtype = SOBOL; 
            else if (rtype == "halton") 
                arg.rtype = HALTON; 
            else if (rtype == "reverse_halton") 
                arg.rtype = REVERSE_HALTON; 
            else if (rtype == "niederreiter_2") 
                arg.rtype = NIEDERREITER_2; 
            else if (rtype == "grid") 
                arg.rtype = GRID; 
            else 
            {
                cerr << "Invalid argument --quasi-random " << rtype << ". "; 
                cerr << "Should be one of the following: sobol, halton, "
                    "reverse_halton, niederreiter_2 or grid. \n"; 
                exit(-1); 
            }
        }
        result_long = parser.getArgLong("--sample-size", 0);
        if (result_long <= 0)
        {
            cerr << "Specify a positive integer for sample size with ";
            cerr << "--sample-size <SAMPLE_SIZE>. \n";
            cerr << usage;
            exit(-1);
        }
        arg.sample_size = static_cast<unsigned>(result_long);

        result_long = parser.getArgLong("--sample-num", 0);
        if (result_long <= 0)
        {
            cerr << "Specify a positive integer for sample num with ";
            cerr << "--sample-num <SAMPLE_SIZE>. \n";
            cerr << usage;
            exit(-1);
        }
        arg.sample_num = static_cast<unsigned>(result_long);
    }
}


template<typename T, unsigned K>
void SampleEntropy()
{
    vector<T> data = ReadData<T>(
        arg.filename, arg.input_format, arg.data_length, arg.line_offset);
    unsigned n = data.size();
    if (n <= K)
    {
        std::cerr << "Data length n " << n << " is too short (K = " << K; 
        std::cerr << "). \n";
        exit(-1);
    }

    double var = ComputeVarience(data);
    T r_scaled = static_cast<T>(sqrt(var) * arg.r);

    if (r_scaled < 0)
    {
        cerr << "Invalid r [scaled]: " << r_scaled << ".\n";
        exit(-1);
    }

    cout.precision(4); 
    cout << std::scientific; 
    cout << "========================================"; 
    cout << "========================================\n"; 
    arg.PrintArguments(); 
    std::cout << "\tvariance: " << var << std::endl; 
    std::cout << "\tr (scaled): " << r_scaled << std::endl; 

    // Compute sample entropy. 
    SampleEntropyCalculatorMao<T, K> sec(data.cbegin(), data.cend(), 
                                         r_scaled, arg.output_level);
    sec.ComputeSampleEntropy(); 
    cout << sec.get_result_str(); 

    if (arg.u) 
    {
        SampleEntropyCalculatorSamplingKDTree<T, K> secs(
            data.cbegin(), data.cend(), r_scaled, 
            arg.sample_size, arg.sample_num, 
            sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), UNIFORM,
            arg.random_, arg.output_level); 
        secs.ComputeSampleEntropy(); 
        cout << secs.get_result_str(); 

        SampleEntropyCalculatorSamplingDirect<T, K> secds(
            data.cbegin(), data.cend(), r_scaled, 
            arg.sample_size, arg.sample_num, 
            sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), UNIFORM, 
            arg.random_, arg.output_level); 
        secds.ComputeSampleEntropy(); 
        cout << secds.get_result_str(); 
    }

    if (arg.q) 
    {
        SampleEntropyCalculatorSamplingKDTree<T, K> secs(
            data.cbegin(), data.cend(), r_scaled, 
            arg.sample_size, arg.sample_num, 
            sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), arg.rtype, 
            arg.random_, arg.output_level); 
        secs.ComputeSampleEntropy(); 
        cout << secs.get_result_str(); 

        SampleEntropyCalculatorSamplingDirect<T, K> secds(
            data.cbegin(), data.cend(), r_scaled, 
            arg.sample_size, arg.sample_num, 
            sec.get_entropy(), sec.get_a_norm(), sec.get_b_norm(), arg.rtype, 
            arg.random_, arg.output_level);
        secds.ComputeSampleEntropy(); 
        cout << secds.get_result_str(); 
    }

    if (arg.fast_direct)
    {
        SampleEntropyCalculatorFastDirect<T, K> secfd(
            data.cbegin(), data.cend(), r_scaled, arg.output_level);
        secfd.ComputeSampleEntropy();
        cout << secfd.get_result_str(); 
    }

    if (arg.direct)
    {
        SampleEntropyCalculatorDirect<T, K> secd(
            data.cbegin(), data.cend(), r_scaled, arg.output_level);
        secd.ComputeSampleEntropy();
        cout << secd.get_result_str(); 
    }
    cout << "========================================"; 
    cout << "========================================\n"; 
}
