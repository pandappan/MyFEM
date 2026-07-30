#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "DenseMatrix.h"
#include "Node.h"
#include "Element/Bar3D.h"
#include "Material/BarMaterial.h"

class Bar3DAxialFixture : public ::testing::Test {
protected:
    std::vector<Node> nodes_;
    std::unique_ptr<BarMaterial> mat_;
    std::unique_ptr<Bar3D> elem_;

    static constexpr double E_    = 1000.0;
    static constexpr double A_    = 2.0;
    static constexpr double L_    = 3.0;   // 沿 X 方向长度

    void SetUp() override {
        // 两个节点，沿 X 方向
        nodes_.emplace_back(0.0, 0.0, 0.0);
        nodes_.emplace_back(L_, 0.0, 0.0);
        nodes_[0].Index = 0;
        nodes_[1].Index = 1;

        // 材料
        mat_.reset(new BarMaterial());
        mat_->nset = 1;
        mat_->E    = E_;
        mat_->Area = A_;

        // 单元
        elem_.reset(new Bar3D());
        std::vector<Node*> nodePtrs = {&nodes_[0], &nodes_[1]};
        elem_->SetupForTesting(nodePtrs, mat_.get());
    }
};
constexpr double Bar3DAxialFixture::E_;
constexpr double Bar3DAxialFixture::A_;
constexpr double Bar3DAxialFixture::L_;
// Bar 刚度矩阵沿 X 方向：K11 = K44 = EA/L, K14 = K41 = -EA/L
// 3D Bar 是 6x6：dof 0,1,2 是 node1 的 XYZ，3,4,5 是 node2
TEST_F(Bar3DAxialFixture, StiffnessMatrixAxialAlignment) {
    DenseMatrix<double> Ke(6, 6);
    elem_->ElementStiffness(Ke);

    const double k = E_ * A_ / L_;   // 解析刚度系数

    // X 方向 (dof 0 和 dof 3)
    EXPECT_NEAR(Ke(0, 0),  k, 1e-8);
    EXPECT_NEAR(Ke(3, 3),  k, 1e-8);
    EXPECT_NEAR(Ke(0, 3), -k, 1e-8);
    EXPECT_NEAR(Ke(3, 0), -k, 1e-8);

    // Y 和 Z 方向应该是 0（因为沿 X 轴）
    EXPECT_NEAR(Ke(1, 1), 0.0, 1e-8);
    EXPECT_NEAR(Ke(2, 2), 0.0, 1e-8);
}
TEST_F(Bar3DAxialFixture, StiffnessMatrixSymmetric) {
    DenseMatrix<double> Ke(6, 6);
    elem_->ElementStiffness(Ke);

    for (unsigned int i = 0; i < 6; ++i)
        for (unsigned int j = i+1; j < 6; ++j)
            EXPECT_NEAR(Ke(i, j), Ke(j, i), 1e-10)
                << "at (" << i << "," << j << ")";
}
// 施加拉伸位移，检查应力
TEST_F(Bar3DAxialFixture, StressUnderPrescribedElongation) {
    // 手动设置节点位移：Node 2 沿 X 方向拉伸 Δ = 0.01
    const double delta = 0.01;
    nodes_[0].displacement[UX] = 0.0;
    nodes_[0].displacement[UY] = 0.0;
    nodes_[0].displacement[UZ] = 0.0;
    nodes_[1].displacement[UX] = delta;
    nodes_[1].displacement[UY] = 0.0;
    nodes_[1].displacement[UZ] = 0.0;
    
    // 解析应力：σ = E·ε = E·(Δ/L)
    const double expected = E_ * delta / L_;
    EXPECT_NEAR(elem_->ElementStress(), expected, 1e-8);
}