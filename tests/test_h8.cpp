#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include "Node.h"
#include "Element/H8.h"
#include "Material/Solid3DMaterial.h"
class H8UnitCubeFixture : public ::testing::Test {
protected:
    std::vector<CNode> nodes_;
    std::unique_ptr<CSolid3DMaterial> mat_;
    std::unique_ptr<CH8> elem_;
    
    void SetUp() override {
        // 单位立方体：8 个节点
        nodes_.emplace_back(0.0, 0.0, 0.0);   // N1
        nodes_.emplace_back(1.0, 0.0, 0.0);   // N2
        nodes_.emplace_back(1.0, 1.0, 0.0);   // N3
        nodes_.emplace_back(0.0, 1.0, 0.0);   // N4
        nodes_.emplace_back(0.0, 0.0, 1.0);   // N5
        nodes_.emplace_back(1.0, 0.0, 1.0);   // N6
        nodes_.emplace_back(1.0, 1.0, 1.0);   // N7
        nodes_.emplace_back(0.0, 1.0, 1.0);   // N8
        for (unsigned int i = 0; i < 8; ++i)
            nodes_[i].NodeNumber = i + 1;
        
        mat_.reset(new CSolid3DMaterial());
        mat_->nset = 1;
        mat_->E    = 1000.0;
        mat_->nu   = 0.3;
        mat_->rho  = 0.0;
        
        elem_.reset(new CH8());
        std::vector<CNode*> nps;
        for (auto& n : nodes_) nps.push_back(&n);
        elem_->SetupForTesting(nps, mat_.get());
        elem_->InitializeIntegrationPoints();
    }
};

TEST_F(H8UnitCubeFixture, JacobianForUnitCube) {
    const auto& ips = elem_->GetIntegrationPoints();
    ASSERT_EQ(ips.size(), 8u);
    // 单位立方体，Jacobian = 0.5·I, detJ = 0.125
    // detJ * w = 0.125 * 1 = 0.125
    for (const auto& ip : ips)
        EXPECT_NEAR(ip.detJ_times_weight, 0.125, 1e-12);

    // 累加 = 1.0 (体积)
    EXPECT_NEAR(elem_->GetVolume(), 1.0, 1e-12);
}

TEST_F(H8UnitCubeFixture, StiffnessMatrixIsSymmetric) {
    DenseMatrix<double> Ke(24, 24);
    elem_->ElementStiffness(Ke);
    for (unsigned int i = 0; i < 24; ++i)
        for (unsigned int j = i + 1; j < 24; ++j)
            EXPECT_NEAR(Ke(i, j), Ke(j, i), 1e-8)
                << "at (" << i << "," << j << ")";
}

TEST_F(H8UnitCubeFixture, RigidBodyTranslation) {
    DenseMatrix<double> Ke(24, 24);
    elem_->ElementStiffness(Ke);

    for (int dir = 0; dir < 3; ++dir) {
        std::vector<double> u(24, 0.0);
        for (int I = 0; I < 8; ++I) u[I * 3 + dir] = 1.0;
        auto f = Ke.DotVec(u);
        for (unsigned int i = 0; i < 24; ++i)
            EXPECT_NEAR(f[i], 0.0, 1e-6)
                << "dir=" << dir << " dof=" << i;
    }
}

// 外推矩阵行和为1
TEST_F(H8UnitCubeFixture, ExtrapolationMatrixRowSumIsOne) {
    const DenseMatrix<double>& E = elem_->GetExprapolationMatrix();
    for (unsigned int i = 0; i < 8; ++i) {
        double rowTotal = 0.0;
        for (unsigned int j = 0; j < 8; ++j) {
            rowTotal += E(i, j);
        }
        EXPECT_NEAR(rowTotal, 1.0, 1e-6);
    }
}

// -------- 均匀 GP 应力外推后节点应力也均匀 --------
TEST_F(H8UnitCubeFixture, UniformStressExtrapolatesToSameValue) {
    auto E = elem_->GetExprapolationMatrix();
    // 给 8 个 GP 都赋 σ = 100
    // 每个节点值 = Σ_gp E[node, gp] · 100 = 100 · (行和=1) = 100
    for (unsigned int i = 0; i < 8; ++i) {
        double val = 0.0;
        for (unsigned int j = 0; j < 8; ++j) val += E(i, j) * 100.0;
        EXPECT_NEAR(val, 100.0, 1e-10) << "node " << i;
    }
}
