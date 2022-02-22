#include "gtest/gtest.h"
#include <vector>

#include "utils.h"

TEST(TestBinarySearchIndex, Test1) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::BinarySearchIndexNoCheck<int>(data, 10);
  EXPECT_EQ(result, 10);
}

TEST(TestBinarySearchIndex, Test2) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::BinarySearchIndexNoCheck<int>(data, 14);
  EXPECT_EQ(result, 11);
}

TEST(TestBinarySearchIndex, Test3) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::BinarySearchIndexNoCheck<int>(data, 4);
  EXPECT_EQ(result, 3);
}

TEST(TestBinarySearchIndex, Test4) {
  std::vector<int> data{0, 1, 8};
  int result = sampen::BinarySearchIndexNoCheck<int>(data, 0);
  EXPECT_EQ(result, 0);
}

TEST(TestBinarySearchIndex, Test5) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::BinarySearchIndexNoCheck<int>(data, 0);
  EXPECT_EQ(result, 0);
}

// Test sampen::CountRangeLastAxis
TEST(TestBinaryCountRange, Large) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::CountRangeLastAxis(4, 10, data);
  EXPECT_EQ(result, 9); 
}

TEST(TestBinaryCountRange, Test2) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::CountRangeLastAxis(4, 4, data);
  EXPECT_EQ(result, 2); 
}

TEST(TestBinaryCountRange, OutOfRangeUpper) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::CountRangeLastAxis(18, 20, data);
  EXPECT_EQ(result, 0); 
}

TEST(TestBinaryCountRange, OutOfRangeLower) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::CountRangeLastAxis(-10, -1, data);
  EXPECT_EQ(result, 0); 
}

TEST(TestBinaryCountRange, All) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::CountRangeLastAxis(-10, 20, data);
  EXPECT_EQ(result, data.size()); 
}

TEST(TestBinaryCountRange, All2) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::CountRangeLastAxis(0, 15, data);
  EXPECT_EQ(result, data.size()); 
}

TEST(TestBinaryCountRange, NearlyAll) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::CountRangeLastAxis(1, 15, data);
  EXPECT_EQ(result, data.size() - 1); 
}

TEST(TestBinaryCountRange, None) {
  std::vector<int> data{0, 1, 4, 4, 5, 8, 9, 9, 10, 10, 10, 13, 15};
  int result = sampen::CountRangeLastAxis(11, 12, data);
  EXPECT_EQ(result, 0); 
}