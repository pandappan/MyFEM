#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include "SkylineMatrix.h"
#include "../build/_deps/googletest-src/googletest/include/gtest/gtest.h"

static std::unique_ptr<SkylineMatrix<double>> Build3x3FullSkyline() {
    auto K = std::unique_ptr<SkylineMatrix<double>>(new SkylineMatrix<double>(3));

    std::vector<unsigned int> lm = {1, 2, 3};
    K->CalculateColumnHeight(lm);
    K->CalculateMaximumHalfBandwidth();
    K->Diagonal();
    K->Allocate();
    return K;
}

// 测试对称性
TEST(SkyLineMatrix, SymmetricAccessSameLocation) {
    auto K = Build3x3FullSkyline();
    (*K)(1,3) = 42.0;
    EXPECT_EQ((*K)(3,1), 42.0);
    (*K)(2,3) = 1.0;
    EXPECT_EQ((*K)(3,2), 1.0);
}

// 测试列高计算
TEST(SkyLineMatrix, ColumnHeightsAreCorrect) {
    auto K = Build3x3FullSkyline();
    const std::vector<unsigned int>& ch = K->GetColumnHeights();
    EXPECT_EQ(ch[0], 0);
    EXPECT_EQ(ch[1], 1);
    EXPECT_EQ(ch[2], 2);
    EXPECT_EQ(K->GetMaximumHalfBandwidth(), 3u);
}

// 测试对角元地址
TEST(SkylineMatrix, DiagonalAddressLayout) {
    auto K = Build3x3FullSkyline();

    const std::vector<unsigned int>& da = K->GetDiagonalAddress();
    // M(0) = 1
    // M(1) = M(0) + H(0) + 1 = 1 + 0 + 1 = 2
    // M(2) = M(1) + H(1) + 1 = 2 + 1 + 1 = 4
    // M(3) = M(2) + H(2) + 1 = 4 + 2 + 1 = 7
    EXPECT_EQ(da[0], 1u);
    EXPECT_EQ(da[1], 2u);
    EXPECT_EQ(da[2], 4u);
    EXPECT_EQ(da[3], 7u);

    // NWK = M(NEQ) - M(0) = 7 - 1 = 6
    EXPECT_EQ(K->size(), 6u);
}

