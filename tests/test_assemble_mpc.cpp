#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include <cmath>
#include "Assembly.h"
#include "Solver.h"
#include "Model.h"
#include "Node.h"
#include "JsonReader.h"
#include "../Core/SkylineMatrix.h"

namespace {
std::string WriteTempJson(const std::string& name, const std::string& content) {
    std::ofstream ofs(name);
    ofs << content;
    ofs.close();
    return name;
}

// 跑完整求解主流程
void Solve(Model& model) {
    Assembler::CalculateEquationNumber(model);
    Assembler::CalculateLocationMatrix(model);
    Assembler::AllocateLinearSystem(model);
    Assembler::InitializeElementMap(model);
    Assembler::AssembleForce(model);
    Assembler::AssembleStiffnessAndConstraintCorrection(model);
    CLDLTSolver solver(*model.K);
    solver.LDLT();
    solver.BackSubstitution(model.force);
    Assembler::WriteDisplacementToNodes(model);
    Assembler::RecoverSlaveDisplacement(model);
}
}

// 等位移 MPC：u_{N2.x} = u_{N3.x}，把 2-3 段变刚性
TEST(AssembleMpcTest, EqualDisplacement) {
    // E*A = 100, 每段长 1 -> k = 100
    // N1 固定，N3.x 受力 F=10
    // 约束 u2 = u3 后，等效单弹簧：u2 = u3 = F/k = 0.1
    const char* jsonText = R"({
  "title": "equal displacement mpc",
  "dimension": 3,
  "nodes": [
    {"id": 1, "x": 0, "y": 0, "z": 0},
    {"id": 2, "x": 1, "y": 0, "z": 0},
    {"id": 3, "x": 2, "y": 0, "z": 0}
  ],
  "materials": [
    {"type": "bar", "E": 100, "area": 1.0}
  ],
  "element_groups": [
    {"type": "Bar3D", "material": 1,
     "elements": [
       {"id": 1, "connectivity": [1, 2]},
       {"id": 2, "connectivity": [2, 3]}
     ]}
  ],
  "boundary_conditions": {
    "fixed": [
      {"nodes": [1], "dof": ["x", "y", "z"]},
      {"nodes": [2, 3], "dof": ["y", "z"]}
    ]
  },
  "loads": {
    "concentrated": [
      {"node": 3, "dof": "x", "value": 10.0}
    ]
  },
  "constraints": {
    "mpc": [
      {
        "slave": {"node": 2, "dof": "x"},
        "masters": [{"node": 3, "dof": "x", "coeff": 1.0}],
        "beta": 0.0
      }
    ]
  }
})";

    Model model;
    JsonReader reader;
    std::string path = WriteTempJson("mpc_equal_disp.json", jsonText);
    ASSERT_TRUE(reader.Read(path, model));

    // N2.x 是 slave
    EXPECT_EQ(model.nodes[1].bcode[UX], 3u);

    Solve(model);

    // master 是 N3.x，应等于 F/k = 0.1
    EXPECT_NEAR(model.nodes[2].displacement[UX], 0.1, 1e-9);
}

// 等位移：slave 位移应等于 master
TEST(RecoverTest, SlaveEqualsMaster) {
    const char* jsonText = R"({
  "title": "recover slave equal master",
  "dimension": 3,
  "nodes": [
    {"id": 1, "x": 0, "y": 0, "z": 0},
    {"id": 2, "x": 1, "y": 0, "z": 0},
    {"id": 3, "x": 2, "y": 0, "z": 0}
  ],
  "materials": [
    {"type": "bar", "E": 100, "area": 1.0}
  ],
  "element_groups": [
    {"type": "Bar3D", "material": 1,
     "elements": [
       {"id": 1, "connectivity": [1, 2]},
       {"id": 2, "connectivity": [2, 3]}
     ]}
  ],
  "boundary_conditions": {
    "fixed": [
      {"nodes": [1], "dof": ["x", "y", "z"]},
      {"nodes": [2, 3], "dof": ["y", "z"]}
    ]
  },
  "loads": {
    "concentrated": [
      {"node": 3, "dof": "x", "value": 10.0}
    ]
  },
  "constraints": {
    "mpc": [
      {
        "slave": {"node": 2, "dof": "x"},
        "masters": [{"node": 3, "dof": "x", "coeff": 1.0}],
        "beta": 0.0
      }
    ]
  }
})";
    Model model;
    JsonReader reader;
    ASSERT_TRUE(reader.Read(WriteTempJson("recover_equal.json", jsonText), model));
    Solve(model);
    double um = model.nodes[2].displacement[UX];  // master N3.x
    double us = model.nodes[1].displacement[UX];  // slave  N2.x
    EXPECT_NEAR(um, 0.1, 1e-9);
    EXPECT_NEAR(us, um, 1e-12);   // slave == master
}
// 带系数：u_slave = 0.5 * u_master
// N1-N2-N3 串联，N1 固定，N3 受力。约束 u2 = 0.5*u3。
// 展开后 2-3 段刚度按 0.5 系数进入 N3 方程。
TEST(RecoverTest, SlaveWithCoefficient) {
    const char* jsonText = R"({
  "title": "recover slave with coeff",
  "dimension": 3,
  "nodes": [
    {"id": 1, "x": 0, "y": 0, "z": 0},
    {"id": 2, "x": 1, "y": 0, "z": 0},
    {"id": 3, "x": 2, "y": 0, "z": 0}
  ],
  "materials": [
    {"type": "bar", "E": 100, "area": 1.0}
  ],
  "element_groups": [
    {"type": "Bar3D", "material": 1,
     "elements": [
       {"id": 1, "connectivity": [1, 2]},
       {"id": 2, "connectivity": [2, 3]}
     ]}
  ],
  "boundary_conditions": {
    "fixed": [
      {"nodes": [1], "dof": ["x", "y", "z"]},
      {"nodes": [2, 3], "dof": ["y", "z"]}
    ]
  },
  "loads": {
    "concentrated": [
      {"node": 3, "dof": "x", "value": 10.0}
    ]
  },
  "constraints": {
    "mpc": [
      {
        "slave": {"node": 2, "dof": "x"},
        "masters": [{"node": 3, "dof": "x", "coeff": 0.5}],
        "beta": 0.0
      }
    ]
  }
})";
    Model model;
    JsonReader reader;
    ASSERT_TRUE(reader.Read(WriteTempJson("recover_coeff.json", jsonText), model));
    Solve(model);
    double um = model.nodes[2].displacement[UX];  // N3.x
    double us = model.nodes[1].displacement[UX];  // N2.x
    EXPECT_NEAR(us, 0.5 * um, 1e-12);   // 关系恒成立
}