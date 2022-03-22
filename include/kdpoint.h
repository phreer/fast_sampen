#ifndef __KDPOINT__
#define __KDPOINT__

#include <iostream>
#include <string.h>
#include <vector>

namespace sampen {
using std::vector;

template <typename T>
class LastAxisTreeNode;

template <typename T>
class RangeKDTree2KNode;

template <typename T> class KDPoint {
public:
  KDPoint(unsigned m = 0): _count(0), _data(m), _value(0) {}
  KDPoint(KDPoint &&point) = default;
  KDPoint(const KDPoint &point) = default;
  KDPoint(typename vector<T>::const_iterator first,
          typename vector<T>::const_iterator last,
          unsigned m, int count = 0)
      : _count(count), _value(0)  {
    const unsigned n = last - first;
    assert(m >= n);
    _data = std::vector<T>(m, static_cast<T>(0));
    for (unsigned i = 0; i < n; ++i) {
      _data[i] = *(first + i);
    }
  }
  KDPoint& operator=(const KDPoint<T> &point) = default;
  int count() const { return _count; }
  int value() const { return _value; }
  void set_value(int value) { _value = value; }
  void set_count(int count) { _count = count; }
  void increase_count(int count) { _count += count; }
  unsigned dim() const { return _data.size(); }
  bool Within(const KDPoint<T> &p, T r, unsigned m) const {
    for (unsigned i = 0; i < m; i++) {
      T diff = _data[i] - p[i];
      if (diff > r || -diff > r)
        return false;
    }
    return true;
  }
  // Overloads operators.
  bool operator<(const KDPoint &p) const {
    const unsigned K = dim();
    for (unsigned i = 0; i < K; i++) {
      if (_data[i] < p[i])
        return true;
      else if (_data[i] > p[i])
        return false;
    }
    return false;
  }
  bool operator==(const KDPoint &p) const {
    const unsigned K = dim();
    for (unsigned i = 0; i < K; i++) {
      if (_data[i] != p[i])
        return false;
    }
    return true;
  }
  const T& operator[](unsigned n) const { return _data[n]; }
  T& operator[](unsigned n) { return _data[n]; }
  void Print() const {
    const unsigned K = dim();
    std::cout << "(";
    for (unsigned i = 0; i < K; i++) {
      std::cout << _data[i] << ", ";
    }
    std::cout << ")\n";
  }

private:
  int _count;
  std::vector<T> _data;
  int _value;
};

template <typename T>
class KDPointRKD: public KDPoint<T> {
public:
  KDPointRKD() = default;
  KDPointRKD(typename vector<T>::const_iterator first,
             typename vector<T>::const_iterator last,
             unsigned m, int count = 0)
      : KDPoint<T>(first, last, m, count), _rank_last_axis(-1),
      _subtree_nodes(0) {}
  KDPointRKD(KDPoint<T> &&point)
      : KDPoint<T>(std::move(point)), _rank_last_axis(-1),
      _subtree_nodes(0) {}
  KDPointRKD(const KDPoint<T> &point)
      : KDPoint<T>(point), _rank_last_axis(-1), _subtree_nodes(0) {}

  void AddSubtreeNode(LastAxisTreeNode<T> *node) {
    _subtree_nodes.push_back(node);
  }
  vector<LastAxisTreeNode<T> *>& subtree_nodes() {
    return _subtree_nodes;
  }
  void set_rank_last_axis(int rank) { _rank_last_axis = rank; }
  int rank_last_axis() const { return _rank_last_axis; }

private:
  // For range kd tree (rkd) method.
  int _rank_last_axis;
  std::vector<LastAxisTreeNode<T> *> _subtree_nodes;
};


} // namespace sampen

#endif // !__KDPOINT__
