#include <gtest/gtest.h>
#include "Assembly.h"
#include "Model.h"
#include "Node.h"
#include "Element/Bar3D.h"

class ExpansionTest : public::testing::Test {
protected:
    Model model;
    CBar3D elem;
    void SetUp() override {
        model.nodes.resize(3);
        // Node A 全部自由
        CNode& nodeA = model.nodes[0];
        nodeA.Index = 0;
        nodeA.bcode[UX] = 0;
        nodeA.bcode[UY] = 0;
        nodeA.bcode[UZ] = 0;
        nodeA.eqn[UX] = 1;
        nodeA.eqn[UY] = 2;
        nodeA.eqn[UZ] = 3;

        CNode& nodeB = model.nodes[1];
        nodeB.Index = 1;
        nodeB.bcode[UX] = 3;
        nodeB.bcode[UY] = 0;
        nodeB.bcode[UZ] = 0;
        nodeB.eqn[UX] = 0;
        nodeB.eqn[UY] = 4;
        nodeB.eqn[UZ] = 5;

        CNode& nodeC = model.nodes[2];
        nodeC.Index = 2;
        nodeC.bcode[UX] = 1; // 固定
        nodeC.bcode[UY] = 2; // 指定位移
        nodeC.bcode[UZ] = 0; // 自由
        nodeC.eqn[UX] = 0;
        nodeC.eqn[UY] = 0;
        nodeC.eqn[UZ] = 6;
        nodeC.Displacement[UY] = 0.05;

        // MPC: u_{B,x} = 2.0 * u_{A,x}(自由) + 3.0 * u_{C.x}(固定) + 0.5 * u_{C,y}(指定位移) + 0.1
        // 期望terms = {eqn 1, 2.0} const = 0.1 + 0.5 * 0.05 = 0.125
        MPC mpc;
        mpc.slaveNode_0 = 1;
        mpc.slaveDof_0 = 0;
        mpc.masters.push_back({0,0,2.0});
        mpc.masters.push_back({2,0,3.0});
        mpc.masters.push_back({2,1,0.5});
        mpc.beta = 0.1;
        model.mpcs.push_back(mpc);
        model.slaveDofToMpc[1 * CNode::NDF + UX] = 0;
        elem.SetupForTesting({&model.nodes[1], &model.nodes[2]},nullptr);
    }
};

// 自由
TEST_F(ExpansionTest, FreeDofMapsToItself) {
    auto e = Assembler::GetLocalDofExpansion(elem, 1, model); // B.y
    ASSERT_EQ(e.terms.size(), 1u);
    EXPECT_EQ(e.terms[0].globalEqn, 4u);
    EXPECT_DOUBLE_EQ(e.terms[0].coeff, 1.0);
    EXPECT_DOUBLE_EQ(e.constant, 0.0);
}

// 固定
TEST_F(ExpansionTest, FixedDofIsEmpty) {
    auto e = Assembler::GetLocalDofExpansion(elem, 3, model); // C.x
    EXPECT_TRUE(e.terms.empty());
    EXPECT_DOUBLE_EQ(e.constant, 0.0);
}

// 指定位移
TEST_F(ExpansionTest, PrescribedDfBecomesConstant) {
    auto e = Assembler::GetLocalDofExpansion(elem, 4, model);
    EXPECT_TRUE(e.terms.empty());
    EXPECT_DOUBLE_EQ(e.constant, 0.05);
}

// slave自由度
TEST_F(ExpansionTest, SlaveExpandsToMaster) {
    auto e = Assembler::GetLocalDofExpansion(elem, 0, model); // B.x
    ASSERT_EQ(e.terms.size(), 1u);
    EXPECT_EQ(e.terms[0].globalEqn, 1u);
    EXPECT_DOUBLE_EQ(e.terms[0].coeff, 2.0);
    EXPECT_DOUBLE_EQ(e.constant, 0.125);
}

TEST_F(ExpansionTest, SlaveWithMissingMpcThrows) {
    // 手工制造一个没有 MPC 记录的从自由度：C.z 强设为 slave
    model.nodes[2].bcode[UZ] = 3;
    EXPECT_THROW(Assembler::GetLocalDofExpansion(elem, 5, model),
                 std::runtime_error);
}