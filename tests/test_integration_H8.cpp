#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include "Model.h"
#include "Element/ElementGroup.h"
#include "Element/H8.h"
#include "Material/Solid3DMaterial.h"
#include "Assembly.h"
#include "Solver.h"

TEST(IntegrationTest, H8_UniaxialTension) {
    // 单位立方体，E=100, ν=0
    // Z 方向拉伸，各节点施加 F/4=1
    // 期望：εzz = 1·L/(E·A) = 1/(100·1·1) = 0.01（每节点 1 力，4 节点合力 4，面积 1）
    //       σzz = F_total/A = 4/1 = 4， εzz = 0.04
    //       Uz(top) = 0.04

    const double E = 100.0;
    const double F_per_node = 1.0;

    Model model;
    model.dimension = 3;

    // 8 节点，同 fixture
    model.nodes.emplace_back(0.0, 0.0, 0.0);
    model.nodes.emplace_back(1.0, 0.0, 0.0);
    model.nodes.emplace_back(1.0, 1.0, 0.0);
    model.nodes.emplace_back(0.0, 1.0, 0.0);
    model.nodes.emplace_back(0.0, 0.0, 1.0);
    model.nodes.emplace_back(1.0, 0.0, 1.0);
    model.nodes.emplace_back(1.0, 1.0, 1.0);
    model.nodes.emplace_back(0.0, 1.0, 1.0);
    for (unsigned int i = 0; i < 8; ++i) model.nodes[i].NodeNumber = i + 1;

    // 底面 z=0 (N1-N4): Uz=0
    // N1: 全约束（防止刚体平移+旋转）
    model.nodes[0].bcode[UX] = 1;
    model.nodes[0].bcode[UY] = 1;
    model.nodes[0].bcode[UZ] = 1;
    // N2: X 方向也约束（防绕 Z 旋转），Uz=0
    model.nodes[1].bcode[UY] = 1;
    model.nodes[1].bcode[UZ] = 1;
    // N3, N4: Uz=0
    model.nodes[2].bcode[UZ] = 1;
    model.nodes[3].bcode[UX] = 1;
    model.nodes[3].bcode[UZ] = 1;

    // 顶面 (N5-N8) 各施加 Fz = 1
    for (int i = 4; i < 8; ++i)
        model.nodes[i].NodeForce[UZ] = F_per_node;

    CElementGroup group;
    group.SetTypeForTesting(ElementTypes::H8);

    auto mat = std::unique_ptr<CSolid3DMaterial>(new CSolid3DMaterial());
    mat->nset = 1; mat->E = E; mat->nu = 0.0;
    group.AddMaterialForTesting(std::move(mat));

    auto elem = std::unique_ptr<CH8>(new CH8());
    std::vector<CNode*> nps;
    for (unsigned int i = 0; i < 8; ++i) nps.push_back(&model.nodes[i]);
    elem->SetupForTesting(nps, &group.GetMaterial(0));
    elem->InitializeIntegrationPoints();
    group.AddElementForTesting(std::move(elem));

    model.groups.push_back(std::move(group));

    // Pipeline
    Assembler::CalculateEquationNumber(model);
    Assembler::CalculateLocationMatrix(model);
    Assembler::AllocateLinearSystem(model);
    Assembler::AssembleForce(model);
    Assembler::AssembleStiffnessAndConstraintCorrection(model);

    CLDLTSolver solver(*model.K);
    solver.LDLT();
    solver.BackSubstitution(model.force);
    Assembler::WriteDisplacementToNodes(model);

    // σzz = 4/1 = 4, εzz = 4/100 = 0.04, Uz(top) = 0.04
    for (int i = 4; i < 8; ++i)
        EXPECT_NEAR(model.nodes[i].Displacement[UZ], 0.04, 1e-8)
            << "Node " << i + 1;

    // ν=0 → x, y 方向无位移（除约束外）
    // Node 6 (1, 0, 1): Ux 自由 → 应约等于 0
    EXPECT_NEAR(model.nodes[5].Displacement[UX], 0.0, 1e-8);

    // 计算单元应变能
    double energy = model.groups[0].GetElement(0).CalculateElementEnergy();
    EXPECT_NEAR(energy, 0.08, 1e-8);
    // 验证能量守恒，Fext*U
    std::vector<double> Fext(24);
    unsigned int index = 0;
    for (auto& node: model.nodes) {
        for (unsigned int i = 0; i < model.dimension; i++) {
            unsigned int dofIndex = index * model.dimension + i;
            Fext[dofIndex] = node.NodeForce[i];
        }
        index++;
    }
    std::vector<double> U(24);
    std::vector<int> nodesBcode(24);
    model.groups[0].GetElement(0).GetElementNodesDisp(U, nodesBcode);
    double extEnergy = 0.0;
    for (unsigned int i = 0; i < 24; i++) {
        extEnergy += Fext[i] * U[i];
    }
    EXPECT_NEAR(extEnergy, 0.16, 1e-8);
    EXPECT_NEAR(2.0 * energy, extEnergy, 1e-8);
}
