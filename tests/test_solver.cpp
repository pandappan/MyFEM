//
// Created by Administrator on 2026/7/8.
//
#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include "SkylineMatrix.h"
#include <Solver.h>

// 手动装一个 3x3 SPD 矩阵 K，验证 LDLT 解 K·u = f
// K = [[ 4, -1,  0],
//      [-1,  4, -1],
//      [ 0, -1,  4]]
// f = [3, 2, 3]
// 解析解 u = [1, 1, 1]
TEST(LDLTSolver, SolvesTridiagonalSystem) {
    // 建 3x3 SkylineMatrix，让 col 1 相连到 col 3
    // 但 K(1,3)=0，所以列高只到 col 2
    // 实际上我们可以只让 col 1-2 相连、col 2-3 相连
    auto K = std::unique_ptr<SkylineMatrix<double>>(new SkylineMatrix<double>(3));

    // 两个"单元"分别连接 [1,2] 和 [2,3]
    std::vector<unsigned int> lm1 = {1, 2};
    std::vector<unsigned int> lm2 = {2, 3};
    K->CalculateColumnHeight(lm1);
    K->CalculateColumnHeight(lm2);
    K->CalculateMaximumHalfBandwidth();
    K->Diagonal();
    K->Allocate();

    // 填入 K（1-based）
    (*K)(1, 1) =  4.0;
    (*K)(2, 2) =  4.0;
    (*K)(3, 3) =  4.0;
    (*K)(1, 2) = -1.0;
    (*K)(2, 3) = -1.0;
    // K(1,3) 不存在（列高不够，也不需要）

    // 右端项 f = [3, 2, 3]
    std::vector<double> f = {3.0, 2.0, 3.0};

    // 求解
    CLDLTSolver solver(*K);
    solver.LDLT();
    solver.BackSubstitution(f);

    // 期望 u = [1, 1, 1]
    EXPECT_NEAR(f[0], 1.0, 1e-10);
    EXPECT_NEAR(f[1], 1.0, 1e-10);
    EXPECT_NEAR(f[2], 1.0, 1e-10);
}