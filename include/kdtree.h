/**
 * @file kdtree.h
 * @author Liu Weifeng
 * @date 2020/05/08
 *
 * @brief Implements the kd tree structure.
 *
 * @details
 */

#ifndef __FAST_SAMPEN_KDTREE__
#define __FAST_SAMPEN_KDTREE__

#include <iostream>
#include <stdlib.h>

#include "utils.h"

namespace sampen {

template <typename T, unsigned K, unsigned D>
Range<T, K> GetRange(typename vector<KDPoint<T, D>>::const_iterator first,
                     typename vector<KDPoint<T, D>>::const_iterator last);

// template<typename T, unsigned K, unsigned D>
// Range<T, K> GetRange(typename vector<KDPointRKD<T, D>>::const_iterator first,
//                      typename vector<KDPointRKD<T, D>>::const_iterator last);

template <typename T, unsigned K>
class KDCountingTreeNode {
public:
  /**
   * @brief Construct kd tree nodes recursively.
   *
   * @param depth: The depth of the current node.
   * @param father: The father node of the current node.
   * @param[out] leaves: The vector of all leaves of this node.
   */
  KDCountingTreeNode(unsigned depth, KDCountingTreeNode *father,
                     vector<KDCountingTreeNode *> &leaves,
                     typename vector<KDPoint<T, K>>::iterator first,
                     typename vector<KDPoint<T, K>>::iterator last);

  ~KDCountingTreeNode() {
    if (_left_child)
      delete _left_child;
    if (_right_child)
      delete _right_child;
  }

  long long CountRange(const Range<T, K> &range, long long &num_nodes) const;

  void UpdateCount(int d) {
    KDCountingTreeNode<T, K> *node = this;
    while (node) {
      node->_weighted_count += d;
      node = node->_father;
    }
  }

  int weighted_count() const { return _weighted_count; }

  unsigned num_nodes() const {
    unsigned result = 1;
    if (_left_child)
      result += _left_child->num_nodes();
    if (_right_child)
      result += _right_child->num_nodes();
    return result;
  }

  void Print();

private:
  unsigned _depth;
  // The number of points.
  unsigned _count;
  int _weighted_count;
  Range<T, K> _range;
  // The weighted count.
  KDCountingTreeNode *_father;
  KDCountingTreeNode *_left_child;
  KDCountingTreeNode *_right_child;
};

template <typename T, unsigned K>
class KDCountingTree {
public:
  KDCountingTree(const vector<KDPoint<T, K>> &points, OutputLevel output_level)
      : _leaves(0), _root(nullptr), _points(points), _index2leaf(points.size()),
        _output_level(output_level) {
    const size_t n = points.size();
    if (n == 0)
      return;

    for (unsigned i = 0; i < n; ++i) {
      _points[i].set_value(i);
    }
    _root = new KDCountingTreeNode<T, K>(0, nullptr, _leaves, _points.begin(),
                                         _points.end());
    for (unsigned i = 0; i < n; ++i) {
      _index2leaf[_points[i].value()] = i;
    }
  }
  ~KDCountingTree() {
    if (_root)
      delete _root;
  }
  long long CountRange(const Range<T, K> &range, long long &num_nodes) const;
  /*
   * @brief Update the counting of a given node and its ancestors.
   *
   * @param position: The index of the leaf node n to be updated.
   * @param d: The value to be added from the leaf node n and the ancestral node
   * of n.
   */
  void UpdateCount(unsigned position, int d);
  void Close(unsigned position);
  unsigned count() const { return _leaves.size(); }
  unsigned num_nodes() const {
    if (_root)
      return _root->num_nodes();
    return 0;
  }
  /// @brief Print the tree recursively.
  void Print();

private:
  vector<KDCountingTreeNode<T, K> *> _leaves;
  KDCountingTreeNode<T, K> *_root;
  vector<KDPoint<T, K>> _points;
  vector<unsigned> _index2leaf;
  OutputLevel _output_level;
};

template <typename T, unsigned K>
class KDCountingTree2KNode {
public:
  /*
   * @brief Construct kd tree nodes recursively.
   *
   * @param depth: The depth of the current node.
   * @param father: The father node of the current node.
   * @param[out] leaves: The vector of all leaves of this node.
   */
  KDCountingTree2KNode(unsigned depth, KDCountingTree2KNode *father,
                       vector<KDCountingTree2KNode *> &leaves,
                       typename vector<KDPoint<T, K>>::iterator first,
                       typename vector<KDPoint<T, K>>::iterator last);

  ~KDCountingTree2KNode() {
    for (unsigned i = 0; i < _num_child; i++) {
      if (_children[i])
        delete _children[i];
    }
  }

  long long CountRange(const Range<T, K> &range, long long &num_nodes,
                       vector<const KDCountingTree2KNode *> &q1,
                       vector<const KDCountingTree2KNode *> &q2) const;

  double CountRangeEstimate(const Range<T, K> &range, long long &num_nodes,
                            vector<const KDCountingTree2KNode *> &q1,
                            vector<const KDCountingTree2KNode *> &q2,
                            unsigned max_depth) const;
  void UpdateCount(int d) {
    KDCountingTree2KNode *node = this;
    while (node) {
      node->_weighted_count += d;
      node = node->_father;
    }
  }

  unsigned num_child() const { return _num_child; }
  int weighted_count() const { return _weighted_count; }

  unsigned num_nodes() const {
    unsigned result = 1;
    for (unsigned i = 0; i < _num_child; i++) {
      result += _children[i]->num_nodes();
    }
    return result;
  }

private:
  unsigned _depth;
  Range<T, K> _range;
  // The number of points.
  unsigned _count;
  // The weighted count.
  int _weighted_count;
  KDCountingTree2KNode *_father;
  vector<KDCountingTree2KNode *> _children;
  unsigned _num_child;
};

template <typename T, unsigned K>
class KDCountingTree2K {
public:
  KDCountingTree2K(const vector<KDPoint<T, K>> &points,
                   OutputLevel output_level)
      : _leaves(0), _root(nullptr), _points(points), _index2leaf(points.size()),
        _q1(points.size()), _q2(points.size()), _output_level(output_level) {
    clock_t t = clock();
    const size_t n = points.size();
    if (n == 0)
      return;

    for (unsigned i = 0; i < n; ++i) {
      _points[i].set_value(i);
    }
    _root = new KDCountingTree2KNode<T, K>(0, nullptr, _leaves, _points.begin(),
                                           _points.end());
    for (unsigned i = 0; i < n; ++i) {
      _index2leaf[_points[i].value()] = i;
    }

    t = clock() - t;
    if (_output_level == Debug) {
      std::cout << "[DEBUG] The time consumed to build a KDCountingTree (K = "
                << K << "): ";
      std::cout << static_cast<double>(t) / CLOCKS_PER_SEC << " seconds. \n";
    }
  }

  ~KDCountingTree2K() {
    if (_root)
      delete _root;
  }

  long long CountRange(const Range<T, K> &range, long long &num_nodes) {
    if (_root)
      return _root->CountRange(range, num_nodes, _q1, _q2);
    return 0;
  }

  /**
   * @brief Update the counting of a given node and its ancestors.
   *
   * @param position: The index of the leaf node n to be updated.
   * @param d: The value to be added from the leaf node n and the ancestral node
   * of n.
   */
  void UpdateCount(unsigned position, int d) {
    assert(position < count() && "position >= count()");
    position = _index2leaf[position];
    if (d)
      _leaves[position]->UpdateCount(d);
  }

  void Close(unsigned position) {
    assert(position < count() && "position >= count()");
    position = _index2leaf[position];
    int w = _leaves[position]->weighted_count();
    if (w != 0)
      _leaves[position]->UpdateCount(-w);
  }

  unsigned count() const { return _leaves.size(); }

  unsigned num_nodes() const {
    if (_root)
      return _root->num_nodes();
    return 0;
  }

private:
  vector<KDCountingTree2KNode<T, K> *> _leaves;
  KDCountingTree2KNode<T, K> *_root;
  vector<KDPoint<T, K>> _points;
  vector<unsigned> _index2leaf;
  vector<const KDCountingTree2KNode<T, K> *> _q1;
  vector<const KDCountingTree2KNode<T, K> *> _q2;
  OutputLevel _output_level;
};

template <typename T, unsigned K>
class KDTree2KNode {
public:
  KDTree2KNode(unsigned depth, KDTree2KNode *father,
               vector<KDTree2KNode *> &leaves,
               typename vector<KDPoint<T, K + 1>>::iterator first,
               typename vector<KDPoint<T, K + 1>>::iterator last,
               unsigned leaf_left);
  ~KDTree2KNode() {
    for (unsigned i = 0; i < _num_child; i++)
      delete _children[i];
  }
  vector<long long> CountRange(const Range<T, K + 1> &range,
                               long long &num_nodes,
                               const vector<KDTree2KNode *> &leaves,
                               vector<const KDTree2KNode *> &q1,
                               vector<const KDTree2KNode *> &q2) const;
  void UpdateCount(int d) {
    KDTree2KNode *node = this;
    while (node) {
      node->_weighted_count += d;
      node = node->_father;
    }
  }
  unsigned count() const { return _count; }
  int weighted_count() const { return _weighted_count; }
  // TODD: This can be optimized.
  unsigned num_child() const { return _num_child; }
  unsigned num_nodes() const {
    unsigned result = 1;
    for (unsigned i = 0; i < num_child(); ++i)
      result += _children[i]->num_nodes();
    return result;
  }

private:
  KDTree2KNode *_father;
  vector<KDTree2KNode *> _children;
  unsigned _depth;
  // The number of points (whether account or not) in this node.
  unsigned _count;
  // Or more precisely, effective_count. Some points in the tree may be closed
  // so as to skipped counting these unnecessary points, without deleting them
  // from the tree.
  int _weighted_count;
  // The index of the first leaf in the current node.
  unsigned _leaf_left;
  Range<T, K> _range;
  unsigned _num_child;
  // Meaningful only for leaf nodes.
  T _last_axis;
};


template <typename T, unsigned K>
class KDTree2K {
public:
  KDTree2K(const vector<KDPoint<T, K + 1>> &points, OutputLevel output_level)
      : _root(nullptr), _leaves(0), _points(points), _index2leaf(points.size()),
        _q1(points.size()), _q2(points.size()), _output_level(output_level) {
    clock_t t = clock();

    const size_t n = points.size();
    if (n == 0)
      return;

    for (unsigned i = 0; i < n; ++i) {
      _points[i].set_value(i);
    }
    _root = new KDTree2KNode<T, K>(0, nullptr, _leaves, _points.begin(),
                                   _points.end(), 0);
    for (unsigned i = 0; i < n; ++i) {
      _index2leaf[_points[i].value()] = i;
    }

    t = clock() - t;
    if (_output_level == Debug) {
      std::cout << "[DEBUG] The time consumed to build a KDCountingTree (K = "
                << K << "): ";
      std::cout << static_cast<double>(t) / CLOCKS_PER_SEC << " seconds. \n";
    }
  }
  ~KDTree2K() {
    if (_root)
      delete _root;
  }
  vector<long long> CountRange(const Range<T, K + 1> &range,
                               long long &num_nodes) {
    if (_root) {
      return _root->CountRange(range, num_nodes, _leaves, _q1, _q2);
    }
    return vector<long long>({0, 0});
  }
  void UpdateCount(unsigned position, int d) {
    assert(position < count() && "position >= count()");
    position = _index2leaf[position];
    if (d)
      _leaves[position]->UpdateCount(d);
  }

  void Close(unsigned position) {
    assert(position < count() && "position >= count()");
    position = _index2leaf[position];
    int w = _leaves[position]->weighted_count();
    if (w != 0)
      _leaves[position]->UpdateCount(-w);
  }

  unsigned count() const { return _leaves.size(); }

  unsigned num_nodes() const {
    if (_root)
      return _root->num_nodes();
    return 0;
  }

private:
  KDTree2KNode<T, K> *_root;
  vector<KDTree2KNode<T, K> *> _leaves;
  vector<KDPoint<T, K + 1>> _points;
  vector<unsigned> _index2leaf;
  vector<const KDTree2KNode<T, K> *> _q1;
  vector<const KDTree2KNode<T, K> *> _q2;
  OutputLevel _output_level;
};


template <typename T>
class LastAxisTreeNode {
public:
  LastAxisTreeNode(typename vector<LastAxisTreeNode<T> >::iterator first,
                   typename vector<LastAxisTreeNode<T> >::iterator last,
                   T lower, T upper, LastAxisTreeNode<T> *parent)
      : _is_leaf(false), _weighted_count(0), _count(last - first),
      _lower(lower), _upper(upper), _parent(parent) {
    const unsigned count = this->count();
    const unsigned median = count / 2;
    if (count > 3) {
      _left_child = new LastAxisTreeNode<T>(
          first, first + median, this->lower(), (first + median - 1)->upper(),
          this);
      _right_child = new LastAxisTreeNode<T>(
          first + median, last, (first + median)->lower(), this->upper(), this);
    } else if (count == 3) {
      first->_parent = this;
      _left_child = &(*first);
      _right_child = new LastAxisTreeNode<T>(
          first + median, last, (first + median)->lower(), this->upper(), this);
    } else if (count == 2) {
      first->_parent = this;
      _left_child = &(*first);
      (first + 1)->_parent = this;
      _right_child = &(*(first + 1));
    } else {
      _left_child = nullptr;
      first->_parent = this;
      _right_child = &(*first);
    }
  }
  LastAxisTreeNode(T value)
      : _is_leaf(true), _weighted_count(0), _count(1),
      _lower(value), _upper(value), _parent(nullptr),
      _left_child(nullptr), _right_child(nullptr) {} 
  ~LastAxisTreeNode() {
    if (_left_child && !_left_child->is_leaf()) {
      delete _left_child;
    }
    if (_right_child && !_right_child->is_leaf()) {
      delete _right_child;
    }
  }
  int CountRange(T lower, T upper) const {
    if (lower <= _lower && _upper <= upper) {
      return weighted_count();
    }
    if (upper < _lower || _upper < lower) {
      return 0;
    }
    return _left_child->CountRange(lower, upper) +\
        _right_child->CountRange(lower, upper);
  }

  T lower() const { return _lower; }
  T upper() const { return _upper; }
  bool is_leaf() const { return _is_leaf; }
  int weighted_count() const { return _weighted_count; }
  unsigned count() const { return _count; }
  void UpdateCount(int x) {
    LastAxisTreeNode<T> *curr = this;
    while (curr) {
      curr->_weighted_count += x;
      curr = curr->_parent;
    }
  }
private:
  bool _is_leaf;
  int _weighted_count;
  unsigned _count;
  T _lower;
  T _upper;
  LastAxisTreeNode<T> *_parent;
  LastAxisTreeNode<T> *_left_child;
  LastAxisTreeNode<T> *_right_child;
};


template <typename T>
class LastAxisTree {
public:
  // Note that nodes must be increasingly ordered.
  LastAxisTree(vector<LastAxisTreeNode<T> > &&nodes)
      : _leaf_nodes(std::move(nodes)), _root(nullptr) {
    if (_leaf_nodes.size() > 1) {
      _root = new LastAxisTreeNode<T>(_leaf_nodes.begin(),
                                      _leaf_nodes.end(),
                                      _leaf_nodes.front().lower(),
                                      _leaf_nodes.back().upper(),
                                      nullptr);
    } else if (_leaf_nodes.size() == 1) {
      _root = &_leaf_nodes[0];
    }
  }
  LastAxisTree(const vector<LastAxisTreeNode<T> > &nodes)
      : _leaf_nodes(nodes), _root(nullptr) {
    if (_leaf_nodes.size() > 1) {
      _root = new LastAxisTreeNode<T>(_leaf_nodes.begin(),
                                      _leaf_nodes.end(),
                                      _leaf_nodes.front().lower(),
                                      _leaf_nodes.back().upper(),
                                      nullptr);
    } else if (_leaf_nodes.size() == 1) {
      _root = &_leaf_nodes[0];
    }
  }
  ~LastAxisTree() {
    if (_root && !_root->is_leaf()) {
      delete _root;
    }
  }
  int CountRange(T lower, T upper) const {
    if (_root) {
      return _root->CountRange(lower, upper);
    }
    return 0;
  }
  const vector<LastAxisTreeNode<T> >& leaf_nodes() const { return _leaf_nodes; }
  vector<LastAxisTreeNode<T> >& leaf_nodes() { return _leaf_nodes; }
  int weighted_count() const { return _root->weighted_count(); }
private:
  vector<LastAxisTreeNode<T> > _leaf_nodes;
  LastAxisTreeNode<T> *_root;
};


template <typename T, unsigned K>
class RangeKDTree2KNode {
public:
  RangeKDTree2KNode(
      unsigned depth, RangeKDTree2KNode *father,
     vector<RangeKDTree2KNode *> &leaves,
     typename vector<KDPointRKD<T, K + 1>>::iterator first,
     typename vector<KDPointRKD<T, K + 1>>::iterator last,
     unsigned leaf_left);
  ~RangeKDTree2KNode() {
    for (unsigned i = 0; i < _num_child; i++)
      delete _children[i];
    if (_subtree) {
      delete _subtree;
    }
  }
  vector<long long> CountRange(const Range<T, K + 1> &range,
                               long long &num_nodes,
                               const vector<RangeKDTree2KNode *> &leaves,
                               vector<const RangeKDTree2KNode *> &q1,
                               vector<const RangeKDTree2KNode *> &q2) const;
  void UpdateCount(int d) {
    RangeKDTree2KNode *node = this;
    while (node) {
      node->_weighted_count += d;
      node = node->_father;
    }
  }
  unsigned count() const { return _count; }
  int weighted_count() const { return _weighted_count; }
  // TODD: This can be optimized.
  unsigned num_child() const { return _num_child; }
  unsigned num_nodes() const {
    unsigned result = 1;
    for (unsigned i = 0; i < num_child(); ++i)
      result += _children[i]->num_nodes();
    return result;
  }

private:
  RangeKDTree2KNode *_father;
  vector<RangeKDTree2KNode *> _children;
  unsigned _depth;
  // The number of points (whether account or not) in this node.
  unsigned _count;
  // Or more precisely, effective_count. Some points in the tree may be closed
  // so as to skipped counting these unnecessary points, without deleting them
  // from the tree.
  int _weighted_count;
  // The index of the first leaf in the current node.
  unsigned _leaf_left;
  Range<T, K> _range;
  unsigned _num_child;
  // For non-leaf nodes for fast searching.
  // vector<T> _last_axis_array;
  LastAxisTree<T> *_subtree;
};


template <typename T, unsigned K>
class RangeKDTree2K {
public:
  RangeKDTree2K(const vector<KDPointRKD<T, K + 1>> &points,
                OutputLevel output_level)
      : _root(nullptr),
        _leaves(0),
        _points(points),
        _index2leaf(points.size()),
        _q1(points.size()),
        _q2(points.size()),
        _output_level(output_level) {
    clock_t t = clock();

    const size_t n = points.size();
    if (n == 0)
      return;
    
    // Sort points according to the last axis and store the rank to the point.
    std::vector<int> order_last_axis(n);
    for (size_t i = 0; i < n; ++i) {
      order_last_axis[i] = i;
    }
    std::sort(order_last_axis.begin(), order_last_axis.end(), 
              [&points] (int i, int j) { return points[i][K] < points[j][K]; });
    for (unsigned i = 0; i < n; ++i) {
      _points[order_last_axis[i]].set_rank_last_axis(i);
    }

    // For mapping from index to leaf.
    for (unsigned i = 0; i < n; ++i) {
      _points[i].set_value(i);
    }
    _root = new RangeKDTree2KNode<T, K>(
        0, nullptr, _leaves, _points.begin(), _points.end(), 0);
    for (unsigned i = 0; i < n; ++i) {
      _index2leaf[_points[i].value()] = i;
    }

    t = clock() - t;
    if (_output_level == Debug) {
      std::cout << "[DEBUG] The time consumed to build a KDCountingTree (K = "
                << K << "): ";
      std::cout << static_cast<double>(t) / CLOCKS_PER_SEC << " seconds. \n";
    }
  }
  ~RangeKDTree2K() {
    if (_root)
      delete _root;
  }
  vector<long long> CountRange(const Range<T, K + 1> &range,
                               long long &num_nodes) {
    if (_root) {
      return _root->CountRange(range, num_nodes, _leaves, _q1, _q2);
    }
    return vector<long long>({0, 0});
  }
  void UpdateCount(unsigned position, int d) {
    assert(position < count() && "position >= count()");
    position = _index2leaf[position];
    if (d) {
      _leaves[position]->UpdateCount(d);
      auto subtree_nodes = _points[position].subtree_nodes();
      for (LastAxisTreeNode<T> *node: subtree_nodes) {
        node->UpdateCount(d);
      }
    }
  }

  void Close(unsigned position) {
    assert(position < count() && "position >= count()");
    position = _index2leaf[position];
    int w = _leaves[position]->weighted_count();
    if (w != 0) {
      _leaves[position]->UpdateCount(-w);
      auto subtree_nodes = _points[position].subtree_nodes();
      for (LastAxisTreeNode<T> *node: subtree_nodes) {
        node->UpdateCount(-w);
      }
    }
  }

  unsigned count() const { return _leaves.size(); }

  unsigned num_nodes() const {
    if (_root)
      return _root->num_nodes();
    return 0;
  }

private:
  RangeKDTree2KNode<T, K> *_root;
  vector<RangeKDTree2KNode<T, K> *> _leaves;
  vector<KDPointRKD<T, K + 1>> _points;
  vector<unsigned> _index2leaf;
  // Buffers for searching without recursion.
  vector<const RangeKDTree2KNode<T, K> *> _q1;
  vector<const RangeKDTree2KNode<T, K> *> _q2;
  OutputLevel _output_level;
};
} // namespace sampen

#endif // !__FAST_SAMPEN_KDTREE__
