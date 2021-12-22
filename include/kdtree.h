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
  unsigned _count;
  int _weighted_count;
  // The index of the first leaf in the current node.
  unsigned _leaf_left;
  Range<T, K> _range;
  unsigned _num_child;
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

} // namespace sampen

#endif // !__FAST_SAMPEN_KDTREE__
