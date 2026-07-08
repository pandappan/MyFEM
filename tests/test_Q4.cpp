//
// Created by Administrator on 2026/7/8.
//
#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include "DenseMatrix.h"
#include "Node.h"
#include "Element/Q4.h"
#include "Material/CPlaneStressMaterial.h"

class Q4UnitSqareFixture : public::testing::Test {
protected:
    std::vector<CNode> nodes_;
    std::unique_ptr<CPlaneStressMaterial> material_;
    std::unique_ptr<CQ4> elem_;

    void SetUp() override {
        nodes_.reserve(4);
        nodes_.emplace_back(0.0,0.0,0.0);
        nodes_.emplace_back(1.0,0.0,0.0);
        nodes_.emplace_back(1.0,1.0,0.0);
        nodes_.emplace_back(0.0,1.0,0.0);

        material_.reset(new CPlaneStressMaterial());
        material_->nset = 1;
        material_->E = 1000.0;
        material_->nu = 0.3;
        material_->thk = 1.0;
        material_->rho = 0.0;

        elem_.reset(new CQ4);
        std::vector<CNode*> nodesPtrs = {&nodes_[0], &nodes_[1], &nodes_[2], &nodes_[3]};
        elem_->SetUpForTesting(nodesPtrs, material_.get());
        elem_->InitializeIntegrationPoints();
    }
};

TEST_F(Q4UnitSqareFixture, FixtureComplies) {
    EXPECT_EQ(nodes_.size(), 4);
    EXPECT_NE(material_, nullptr);
    EXPECT_NE(elem_, nullptr);
}

TEST_F(Q4UnitSqareFixture, JacobianDeterminantIsQuarterOfArea) {
    const auto& ips = elem_->GetIntegrationPoints();
    ASSERT_EQ(ips.size(), 4u);

    for (const auto& ip: ips) {
        EXPECT_NEAR(ip.detJ_times_weight, 0.25, 1e-12);
    }

    double totalArea = 0.0;
    for (const auto& ip: ips) totalArea += ip.detJ_times_weight;
    EXPECT_NEAR(totalArea, 1.0, 1.0e-12);
}

TEST_F(Q4UnitSqareFixture, StiffnessMatrixIsSymmetric) {
    DenseMatrix<double> K(8, 8);
    elem_->ElementStiffness(K);
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            EXPECT_NEAR(K(i, j), K(j, i), 1.0e-12)
            << "Asymmetric at (" << i << ", " << j << " )" << std::endl;
        }
    }
}

TEST_F(Q4UnitSqareFixture, RigidBodyTranslationX) {
    DenseMatrix<double> K(8, 8);
    elem_->ElementStiffness(K);

    std::vector<double> u(8);
    for (int i = 0; i < 4; ++i) {
        u[2 * i] = 1.0;
        u[2 * i + 1] = 0.0;
    }

    std::vector<double> Fin = K.DotVec(u);

    for (int i = 0; i < 8; i++) {
        EXPECT_NEAR(Fin[i], 0.0, 1.0e-12);
    }
}

TEST_F(Q4UnitSqareFixture, RigidBodyTranslationY) {
    DenseMatrix<double> K(8, 8);
    elem_->ElementStiffness(K);

    std::vector<double> u(8);
    for (int i = 0; i < 4; ++i) {
        u[2 * i] = 0.0;
        u[2 * i + 1] = 1.0;
    }

    std::vector<double> Fin = K.DotVec(u);

    for (int i = 0; i < 8; i++) {
        EXPECT_NEAR(Fin[i], 0.0, 1.0e-12);
    }
}

TEST_F(Q4UnitSqareFixture, RigidBodyRoation) {
    DenseMatrix<double> K(8, 8);
    elem_->ElementStiffness(K);

    std::vector<double> u = {
        0.0, 0.0,   // node 1
        0.0, 1.0,   // node 2
        -1.0, 1.0,   // node 3
        -1.0, 0.0   // node 4
    };

    std::vector<double> Fin = K.DotVec(u);

    for (int i = 0; i < 8; i++) {
        EXPECT_NEAR(Fin[i], 0.0, 1.0e-12);
    }
}

