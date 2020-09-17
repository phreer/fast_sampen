#ifndef __SAMPLE_ENTROPY_CALCULATOR__
#define __SAMPLE_ENTROPY_CALCULATOR__

#include <vector>
#include <type_traits>

#include "utils.h"

namespace kdtree_mddc
{
using std::vector;

template<typename T, unsigned K, 
         typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
class SampleEntropyCalculator
{
public:
    SampleEntropyCalculator(OutputLevel output_level) : _output_level(output_level) {}
    virtual double ComputeSampleEntropy(typename vector<T>::const_iterator first,
                                        typename vector<T>::const_iterator last,
                                        T r) = 0;
protected:
    OutputLevel _output_level;
};

} // namespace kdtree_mddc

#endif // !__SAMPLE_ENTROPY_CALCULATOR__