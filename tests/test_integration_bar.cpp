//
// Created by Administrator on 2026/7/8.
//
#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include <iostream>
#include "Model.h"
#include "Node.h"
#include "Element/Bar3D.h"
#include "Element/ElementGroup.h"
#include "Material/BarMaterial.h"
#include "Assembly.h"
#include "Solver.h"

TEST(IntegrationTest, BarAxialTension_UsingFullPipeline) {
    // -------- 参数 --------
    const double E = 210e9;   // Pa
    const double A = 1e-4;    // m^2
    const double L = 2.0;     // m
    const double F = 1000.0;  // N
    const double expected_ux = F * L / (E * A);
    const double expected_stress = F / A;
    const double extEnergy = F * expected_ux;

    // -------- 构建 Model --------
    Model model;
    model.dimension = 3;

    // 节点：注意 CNode 默认构造后 bcode 全 0（自由）
    model.nodes.emplace_back(0.0, 0.0, 0.0);
    model.nodes.emplace_back(L,   0.0, 0.0);
    model.nodes[0].NodeNumber = 1;
    model.nodes[1].NodeNumber = 2;

    // Node 1：全固定
    model.nodes[0].bcode[UX] = 1;
    model.nodes[0].bcode[UY] = 1;
    model.nodes[0].bcode[UZ] = 1;

    // Node 2：X 自由，YZ 固定
    model.nodes[1].bcode[UX] = 0;
    model.nodes[1].bcode[UY] = 1;
    model.nodes[1].bcode[UZ] = 1;

    // 载荷：X 方向 F
    model.nodes[1].NodeForce[UX] = F;

    // -------- 建 Element Group --------
    CElementGroup group;
    group.SetTypeForTesting(ElementTypes::Bar3D);

    // 材料
    auto mat = std::unique_ptr<CBarMaterial>(new CBarMaterial());
    mat->nset = 1;
    mat->E    = E;
    mat->Area = A;
    group.AddMaterialForTesting(std::move(mat));

    // 单元
    auto elem = std::unique_ptr<CBar3D>(new CBar3D());
    std::vector<CNode*> nodes = {&model.nodes[0], &model.nodes[1]};
    elem->SetupForTesting(
        nodes,
        &group.GetMaterial(0)
    );
    group.AddElementForTesting(std::move(elem));

    model.groups.push_back(std::move(group));

    // -------- 执行完整分析流程 --------
    Assembler::CalculateEquationNumber(model);
    ASSERT_EQ(model.neq, 1u) << "Only 1 free DOF expected";

    Assembler::CalculateLocationMatrix(model);
    Assembler::AllocateLinearSystem(model);
    Assembler::AssembleForce(model);
    Assembler::AssembleStiffnessAndConstraintCorrection(model);

    // 验证：Force 向量应该有 F
    ASSERT_EQ(model.force.size(), 1u);
    EXPECT_NEAR(model.force[0], F, 1e-6);

    // 求解
    CLDLTSolver solver(*model.K);
    solver.LDLT();
    solver.BackSubstitution(model.force);

    Assembler::WriteDisplacementToNodes(model);

    // -------- 断言结果 --------
    // Node 2 X 位移 = F·L/(E·A)
    EXPECT_NEAR(model.nodes[1].Displacement[UX], expected_ux, 1e-10);

    // Node 1 位移 = 0（约束）
    EXPECT_DOUBLE_EQ(model.nodes[0].Displacement[UX], 0.0);
    EXPECT_DOUBLE_EQ(model.nodes[0].Displacement[UY], 0.0);
    EXPECT_DOUBLE_EQ(model.nodes[0].Displacement[UZ], 0.0);

    // 单元应力 = F/A
    auto* bar = dynamic_cast<CBar3D*>(&model.groups[0].GetElement(0));
    ASSERT_NE(bar, nullptr);
    EXPECT_NEAR(bar->ElementStress(), expected_stress, 1e-4);

    // 计算单元应变能
    double energy = model.groups[0].GetElement(0).CalculateElementEnergy();
    EXPECT_NEAR(2.0 * energy, extEnergy, 1e-8);
}