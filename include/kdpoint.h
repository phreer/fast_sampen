#ifndef __KDPOINT__
#define __KDPOINT__

#include <iostream>
#include <string.h>
#include <vector>

namespace sampen {
using std::vector;

template <typename T>
class LastAxisTreeNode;

template <typename T, unsigned K>
class RangeKDTree2KNode;

template <typename T, unsigned K> class KDPoint {
public:
  KDPoint() = default;
  KDPoint(KDPoint &&point) = default;
  KDPoint(const KDPoint &point) = default;
  KDPoint(typename vector<T>::const_iterator first,
          typename vector<T>::const_iterator last, int count = 0);
  KDPoint& operator=(const KDPoint<T, K> &point) = default;
  int count() const { return _count; }
  int value() const { return _value; }
  void set_value(int value) { _value = value; }
  void set_count(int count) { _count = count; }
  void increase_count(int count) { _count += count; }
  unsigned dim() const { return K; }
  bool Within(const KDPoint<T, K> &p, T r, unsigned m = K) const {
    for (unsigned i = 0; i < m; i++) {
      T diff = _data[i] - p[i];
      if (diff > r || -diff > r)
        return false;
    }
    return true;
  }
  // Overloads operators.
  bool operator<(const KDPoint &p) const {
    for (unsigned i = 0; i < K; i++) {
      if (_data[i] < p[i])
        return true;
      else if (_data[i] > p[i])
        return false;
    }
    return false;
  }
  bool operator==(const KDPoint &p) const {
    for (unsigned i = 0; i < K; i++) {
      if (_data[i] != p[i])
        return false;
    }
    return true;
  }
  T operator[](unsigned n) const { return _data[n]; }
  T &operator[](unsigned n) { return _data[n]; }
  void Print() const {
    std::cout << "(";
    for (unsigned i = 0; i < K; i++) {
      std::cout << _data[i] << ", ";
    }
    std::cout << ")\n";
  }

private:
  int _count;
  T _data[K];
  int _value;
};

template <typename T, unsigned K>
class KDPointRKD: public KDPoint<T, K> {
public:
  KDPointRKD() = default;
  KDPointRKD(typename vector<T>::const_iterator first,
             typename vector<T>::const_iterator last, int count = 0)
      : KDPoint<T, K>(first, last, count), _rank_last_axis(-1),
      _subtree_nodes(0), _rkd_node(nullptr) {}
  KDPointRKD(KDPoint<T, K> &&point)
      : KDPoint<T, K>(std::move(point)), _rank_last_axis(-1),
      _subtree_nodes(0), _rkd_node(nullptr) {}
  KDPointRKD(const KDPoint<T, K> &point)
      : KDPoint<T, K>(point), _rank_last_axis(-1),
      _subtree_nodes(0), _rkd_node(nullptr) {}

  void AddSubtreeNode(LastAxisTreeNode<T> *node) {
    _subtree_nodes.push_back(node);
  }
  vector<LastAxisTreeNode<T> *>& subtree_nodes() {
    return _subtree_nodes;
  }
  void set_rkd_node(RangeKDTree2KNode<T, K - 1> *node) { _rkd_node = node; }
  RangeKDTree2KNode<T, K - 1>* rkd_node() { return _rkd_node; }
  void set_rank_last_axis(int rank) { _rank_last_axis = rank; }
  int rank_last_axis() const { return _rank_last_axis; }

private:
  // For range kd tree (rkd) method.
  int _rank_last_axis;
  std::vector<LastAxisTreeNode<T> *> _subtree_nodes;
  RangeKDTree2KNode<T, K - 1> *_rkd_node;
};

//////////////////////////////////////////////////////////////////////////
// Implementations
//////////////////////////////////////////////////////////////////////////

template <typename T, unsigned K>
KDPoint<T, K>::KDPoint(typename vector<T>::const_iterator first,
                       typename vector<T>::const_iterator last,
                       int count)
    : _count(count) {
  memset(_data, 0, K * sizeof(unsigned));
  unsigned n = K < (last - first) ? K : (last - first);
  std::copy(first, first + n, _data);
}


} // namespace sampen

#endif // !__KDPOINT__
