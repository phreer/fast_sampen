#ifndef __KDPOINT__
#define __KDPOINT__

#include <vector>
#include <iostream>

namespace kdtree_mddc
{
using std::vector;

template<typename T, unsigned K>
class KDPoint
{
public:
    KDPoint() = default;
    KDPoint(typename vector<T>::const_iterator first,
            typename vector<T>::const_iterator last,
            int count = 0);
    int count() const { return _count; }
    int value() const { return _value; }
    void set_value(int value) { _value = value; }
    void set_count(int count) { _count = count; }
    void increase_count(int count) { _count += count; }
    unsigned dim() const { return K; }
    bool Within(const KDPoint<T, K> &p, T r, unsigned m = K)
    {
        for (unsigned i = 0; i < m; i++)
        {
            T diff = _data[i] - p[i];
            if (diff > r || -diff > r) return false;
        }
        return true;
    }
    // Overloads operators. 
    bool operator < (const KDPoint &p) const
    {
        for (unsigned i = 0; i < K; i++)
        {
            if (_data[i] < p[i]) return true;
            else if (_data[i] > p[i]) return false;
        }
        return false;
    }
    bool operator == (const KDPoint &p) const
    {
        for (unsigned i = 0; i < K; i++)
        {
            if (_data[i] != p[i]) return false;
        }
        return true;
    }
    T operator [] (unsigned n) const
    {
        return _data[n];
    }
    T& operator [] (unsigned n)
    {
        return _data[n];
    }
    void Print() const
    {
        std::cout << "(";
        for (unsigned i = 0; i < K; i++)
        {
            std::cout << _data[i] << ", ";
        }
        std::cout << ")\n";
    }
private:
    int _count;
    T _data[K];
    int _value;
};

//////////////////////////////////////////////////////////////////////////
// Implementations 
//////////////////////////////////////////////////////////////////////////

template<typename T, unsigned K>
KDPoint<T, K>::KDPoint(typename vector<T>::const_iterator first,
                       typename vector<T>::const_iterator last,
                       int count)
    : _count(count)
{
    memset(_data, 0, K * sizeof(unsigned));
    unsigned n = K < (last - first) ? K : (last - first);
    std::copy(first, first + n, _data);
}
} //namespace kdtree_mddc

#endif // !__KDPOINT__