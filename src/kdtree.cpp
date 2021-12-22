//
// Created by Phree on 2021/12/9.
//
#include "kdtree.h"
namespace sampen {

template<typename T, unsigned K>
void KDCountingTree<T, K>::Print() {
  if (_root) {
    std::cout << "KD Counting Tree: " << std::endl;
    _root->Print();
  } else {
    std::cout << "Empty KD Counting Tree. " << std::endl;
  }
}

template<typename T, unsigned K, unsigned D>
Range<T, K> GetRange(typename vector<KDPoint<T, D>>::const_iterator first,
                     typename vector<KDPoint<T, D>>::const_iterator last) {
  const size_t n = last - first;
  assert(n > 0);
  Range<T, K> range;

  T maximum[K];
  T minimum[K];
  for (unsigned i = 0; i < K; ++i) {
    minimum[i] = (*first)[i];
    maximum[i] = minimum[i];
  }
  for (size_t j = 0; j < (n - 1) / 2; ++j) {
    const KDPoint<T, D> &point1 = *(first + 2 * j + 1);
    const KDPoint<T, D> &point2 = *(first + 2 * (j + 1));
    for (unsigned i = 0; i < K; ++i) {
      T curr1 = point1[i];
      T curr2 = point2[i];
      if (curr1 < curr2) {
        if (maximum[i] < curr2)
          maximum[i] = curr2;
        if (minimum[i] > curr1)
          minimum[i] = curr1;
      } else {
        if (maximum[i] < curr1)
          maximum[i] = curr1;
        if (minimum[i] > curr2)
          minimum[i] = curr2;
      }
    }
  }

// last one
  if (n % 2 == 0) {
    const KDPoint<T, D> &point = *(first + n - 1);
    for (unsigned i = 0; i < K; ++i) {
      T curr = point[i];
      if (maximum[i] < curr)
        maximum[i] = curr;
      if (minimum[i] > curr)
        minimum[i] = curr;
    }
  }
  for (unsigned i = 0; i < K; ++i) {
    range.lower_ranges[i] = minimum[i];
    range.upper_ranges[i] = maximum[i];
  }
  return range;
}

template<typename T, unsigned K>
long long KDCountingTree<T, K>::CountRange(const Range<T, K> &range,
                                                  long long &num_nodes) const {
  if (!_root)
    return 0;

  return _root->CountRange(range, num_nodes);
}

template<typename T, unsigned K>
void KDCountingTree<T, K>::UpdateCount(unsigned position, int d) {
  assert(position < count() && "position >= count()");
  position = _index2leaf[position];
  if (d)
    _leaves[position]->UpdateCount(d);
}

template<typename T, unsigned K>
void KDCountingTree<T, K>::Close(unsigned position) {
  assert(position < count() && "position >= count()");
  position = _index2leaf[position];
  int w = _leaves[position]->weighted_count();
  if (w != 0)
    _leaves[position]->UpdateCount(-w);
}

template<typename T, unsigned K>
void KDCountingTreeNode<T, K>::Print() {
  std::cout << "Tree Node: " << std::endl;
  std::cout << "\tcount: " << _count << std::endl;
  std::cout << "\tweighted_count: " << _weighted_count << std::endl;
  std::cout << "\tdepth: " << _depth << std::endl;
  std::cout << "\tranges: [";
  for (unsigned i = 0; i < K; ++i) {
    std::cout << _range.lower_ranges[i] << ", " << _range.upper_ranges[i];
    std::cout << "; ";
  }
  std::cout << "]\n\n";

  if (_left_child)
    _left_child->Print();
  if (_right_child)
    _right_child->Print();
}

template<typename T, unsigned K>
KDCountingTreeNode<T, K>::KDCountingTreeNode(
    unsigned depth, KDCountingTreeNode *father,
    vector<KDCountingTreeNode *> &leaves,
    typename vector<KDPoint<T, K>>::iterator first,
    typename vector<KDPoint<T, K>>::iterator last)
    : _depth(depth), _count(last - first), _weighted_count(0), _father(father),
      _left_child(nullptr), _right_child(nullptr) {
  _range = GetRange<T, K, K>(first, last);
  if (_count == 1) {
    leaves.push_back(this);
    return;
  }

  std::nth_element(first, first + _count / 2, last,
                   [depth](const KDPoint<T, K> &p1, const KDPoint<T, K> &p2) {
                     const unsigned dim = depth % K;
                     return p1[dim] < p2[dim];
                   });

  _left_child = new KDCountingTreeNode<T, K>(_depth + 1, this, leaves, first,
                                             first + _count / 2);
  _right_child = new KDCountingTreeNode<T, K>(_depth + 1, this, leaves,
                                              first + _count / 2, last);
}

template<typename T, unsigned K>
long long KDCountingTreeNode<T, K>::CountRange(const Range<T, K> &range,
                                               long long &num_nodes) const {
  if (_weighted_count == 0)
    return 0;

  enum CASE { NOT_INTER, WITHIN, INTER };

  long long result = 0;

  // Queue of nodes to count.
  auto q1 =
      static_cast<const KDCountingTreeNode **>(malloc(_count * sizeof(void *)));
  auto q2 =
      static_cast<const KDCountingTreeNode **>(malloc(_count * sizeof(void *)));
  if (!q1 || !q2) {
    std::cerr << "Unable to allocate memory, size: ";
    std::cerr << _count * sizeof(void *) << std::endl;
    exit(-1);
  }
  q1[0] = this;
  unsigned n1 = 1, n2 = 0;

  T a, b, c, d;
  while (n1) {
    num_nodes += n1;
    for (unsigned j = 0; j < n1; j++) {
      const KDCountingTreeNode *curr = q1[j];
      enum CASE _case = WITHIN;
      for (unsigned i = 0; i < K; i++) {
        a = curr->_range.lower_ranges[i];
        b = curr->_range.upper_ranges[i];
        c = range.lower_ranges[i];
        d = range.upper_ranges[i];
        if (a > d || b < c) {
          _case = NOT_INTER;
          break;
        }
        if (a < c || b > d) {
          _case = INTER;
        }
      }

      switch (_case) {
        case WITHIN: {
          result += static_cast<long long>(curr->_weighted_count);
          break;
        }
        case INTER: {
          if (curr->_left_child) {
            if (curr->_left_child->_weighted_count > 0) {
              q2[n2] = curr->_left_child;
              ++n2;
            }
          };
          if (curr->_right_child) {
            if (curr->_right_child->_weighted_count > 0) {
              q2[n2] = curr->_right_child;
              ++n2;
            }
          }
          break;
        }
        default: break;
      }
    }
    std::swap(q1, q2);
    n1 = n2;
    n2 = 0;
  }

  free(q2);
  free(q1);

  return result;
}

template<typename T, unsigned K>
KDCountingTree2KNode<T, K>::KDCountingTree2KNode(
    unsigned depth, KDCountingTree2KNode *father,
    vector<KDCountingTree2KNode *> &leaves,
    typename vector<KDPoint<T, K>>::iterator first,
    typename vector<KDPoint<T, K>>::iterator last)
    : _depth(depth), _count(last - first), _weighted_count(0), _father(father) {
  _range = GetRange<T, K, K>(first, last);
  if (_count == 1) {
    _num_child = 0;
    leaves.push_back(this);
    return;
  }

  unsigned splitters[1u << (K + 1)];
  splitters[0] = 0;
  splitters[1u << K] = _count;

  unsigned median, splitter1, splitter2;
  for (unsigned i = 0; i < K; i++) {
    const unsigned spacing = 1u << (K - i);
    for (unsigned j = 0; j < (1u << i); j++) {
      splitter1 = splitters[j * spacing];
      splitter2 = splitters[(j + 1) * spacing];

      median = splitter1 + (splitter2 - splitter1) / 2;
      splitters[j * spacing + spacing / 2] = median;
      std::nth_element(first + splitter1, first + median, first + splitter2,
                       [&i](const KDPoint<T, K> &p1, const KDPoint<T, K> &p2) {
                         return p1[i] < p2[i];
                       });
    }
  }

  unsigned k = 0;
  for (unsigned i = 0; i < (1u << K); i++) {
    splitter1 = splitters[i];
    splitter2 = splitters[i + 1];
    if (splitter1 != splitter2) {
      KDCountingTree2KNode<T, K> *child = new KDCountingTree2KNode<T, K>(
          _depth + 1, this, leaves, first + splitter1, first + splitter2);
      _children.push_back(child);
      k++;
    }
  }
  _num_child = k;
}

template<typename T, unsigned K>
long long KDCountingTree2KNode<T, K>::CountRange(
    const Range<T, K> &range, long long &num_nodes,
    vector<const KDCountingTree2KNode *> &q1,
    vector<const KDCountingTree2KNode *> &q2) const {
  if (weighted_count() == 0)
    return 0;
  enum CASE { NOT_INTER, WITHIN, INTER };

  long long result = 0;
  // Nodes to count.
  q1[0] = this;
  unsigned n1 = 1, n2 = 0;

  T a, b, c, d;
  while (n1) {
    num_nodes += n1;
    for (unsigned j = 0; j < n1; j++) {
      const KDCountingTree2KNode *curr = q1[j];
      enum CASE _case = WITHIN;
      for (unsigned i = 0; i < K; ++i) {
        a = curr->_range.lower_ranges[i];
        b = curr->_range.upper_ranges[i];
        c = range.lower_ranges[i];
        d = range.upper_ranges[i];
        if (a > d || b < c) {
          _case = NOT_INTER;
          break;
        }
        if (a < c || b > d) {
          _case = INTER;
        }
      }

      switch (_case) {
        case WITHIN: {
          result += static_cast<long long>(curr->_weighted_count);
          break;
        }
        case INTER: {
          for (unsigned i = 0; i < curr->num_child(); ++i) {
            // This criterion is critical!
            if (curr->_children[i]->_weighted_count) {
              q2[n2] = curr->_children[i];
              ++n2;
            }
          }
          break;
        }
        case NOT_INTER:
        default:break;
      }
    }
    std::swap(q1, q2);
    n1 = n2;
    n2 = 0;
  }
  return result;
}

template<typename T, unsigned K>
KDTree2KNode<T, K>::KDTree2KNode(
    unsigned depth, KDTree2KNode *father, vector<KDTree2KNode *> &leaves,
    typename vector<KDPoint<T, K + 1>>::iterator first,
    typename vector<KDPoint<T, K + 1>>::iterator last, unsigned leaf_left)
    : _father(father), _depth(depth), _count(last - first), _weighted_count(0),
      _leaf_left(leaf_left) {
  assert(_count > 0);
  _range = GetRange<T, K, K + 1>(first, last);

  if (_count == 1) {
    _num_child = 0;
    _last_axis = (*first)[K];
    leaves.push_back(this);
    return;
  }

  unsigned splitters[1u << (K + 1)];
  splitters[0] = 0;
  splitters[1u << K] = _count;

  unsigned median, splitter1, splitter2;
  for (unsigned i = 0; i < K; i++) {
    const unsigned spacing = 1u << (K - i);
    for (unsigned j = 0; j < (1u << i); j++) {
      splitter1 = splitters[j * spacing];
      splitter2 = splitters[(j + 1) * spacing];

      median = splitter1 + (splitter2 - splitter1) / 2;
      splitters[j * spacing + spacing / 2] = median;
      std::nth_element(
          first + splitter1, first + median, first + splitter2,
          [&i](const KDPoint<T, K + 1> &p1, const KDPoint<T, K + 1> &p2) {
            return p1[i] < p2[i];
          });
    }
  }

  unsigned k = 0;
  for (unsigned i = 0; i < (1u << K); i++) {
    splitter1 = splitters[i];
    splitter2 = splitters[i + 1];
    if (splitter1 != splitter2) {
      KDTree2KNode<T, K> *child =
          new KDTree2KNode<T, K>(_depth + 1, this, leaves, first + splitter1,
                                 first + splitter2, leaf_left + splitter1);
      _children.push_back(child);
      k++;
    }
  }
  _num_child = k;
}

// Non-recursive version.
template<typename T, unsigned K>
vector<long long> KDTree2KNode<T, K>::CountRange(
    const Range<T, K + 1> &range, long long &num_nodes,
    const vector<KDTree2KNode *> &leaves, vector<const KDTree2KNode *> &q1,
    vector<const KDTree2KNode *> &q2) const {
  vector<long long> result({0, 0});
  if (weighted_count() == 0)
    return result;

  enum CASE { NOT_INTER, WITHIN, INTER };

  // Nodes to count.
  q1[0] = this;
  unsigned n1 = 1, n2 = 0;

  T a, b, c, d;
  while (n1) {
    num_nodes += n1;
    for (unsigned j = 0; j < n1; j++) {
      const KDTree2KNode *curr = q1[j];
      enum CASE _case = WITHIN;
      for (unsigned i = 0; i < K; ++i) {
        a = curr->_range.lower_ranges[i];
        b = curr->_range.upper_ranges[i];
        c = range.lower_ranges[i];
        d = range.upper_ranges[i];
        if (a > d || b < c) {
          _case = NOT_INTER;
          break;
        }
        if (a < c || b > d) {
          _case = INTER;
        }
      }

      switch (_case) {
        case WITHIN: {
          result[1] += static_cast<long long>(curr->_weighted_count);
          // Check last coordinate.
          for (unsigned i = 0; i < curr->_count; ++i) {
            if (leaves[curr->_leaf_left + i]->weighted_count() == 0)
              continue;
            T last_axis = leaves[curr->_leaf_left + i]->_last_axis;
            if (range.lower_ranges[K] <= last_axis &&
                last_axis <= range.upper_ranges[K]) {
              result[0] += 1;
            }
          }
          break;
        }
        case INTER: {
          for (unsigned i = 0; i < curr->num_child(); ++i) {
            // This criterion is critical!
            if (curr->_children[i]->_weighted_count) {
              q2[n2] = curr->_children[i];
              ++n2;
            }
          }
          break;
        }
        case NOT_INTER:
        default:break;
      }
    }
    std::swap(q1, q2);
    n1 = n2;
    n2 = 0;
  }

  return result;
}

#define INSTANTIATE_KDTREE(K) \
template class KDCountingTree2K<int, K>; \
template class KDCountingTree<int, K>; \
template class KDTree2K<int, K>; \
template class KDCountingTree2K<double, K>; \
template class KDCountingTree<double, K>; \
template class KDTree2K<double, K>; \
template class KDCountingTree2K<unsigned, K>; \
template class KDCountingTree<unsigned, K>; \
template class KDTree2K<unsigned, K>;


INSTANTIATE_KDTREE(1)
INSTANTIATE_KDTREE(2)
INSTANTIATE_KDTREE(3)
INSTANTIATE_KDTREE(4)
INSTANTIATE_KDTREE(5)
INSTANTIATE_KDTREE(6)
INSTANTIATE_KDTREE(7)
INSTANTIATE_KDTREE(8)
INSTANTIATE_KDTREE(9)
INSTANTIATE_KDTREE(10)
INSTANTIATE_KDTREE(11)
} // namespace sampen