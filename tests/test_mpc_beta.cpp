#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include "Assembly.h"
#include "Solver.h"
#include "Model.h"
#include "Node.h"
#include "JsonReader.h"

namespace {
std::string WriteTempJson(const std::string& name, const std::string& content) {
    std::ofstream ofs(name);
    ofs << content;
    ofs.close();
    return name;
}

void Solve(Model& model) {
    Assembler::CalculateEquationNumber(model);
    Assembler::CalculateLocationMatrix(model);
    Assembler::AllocateLinearSystem(model);
    Assembler::InitializeElementMap(model);
    Assembler::AssembleForce(model);
    Assembler::AssembleStiffnessAndConstraintCorrection(model);
    // neq 可能为 0 的极端情形下跳过求解
    if (model.neq > 0) {
        CLDLTSolver solver(*model.K);
        solver.LDLT();
        solver.BackSubstitution(model.force);
    }
    Assembler::WriteDisplacementToNodes(model);
    Assembler::RecoverSlaveDisplacement(model);
}
}

// beta 常数项：u_{N2.x} = u_{N3.x} + beta，N3 固定 -> u2 = beta
TEST(BetaTest, SlaveDrivenByBeta) {
    const char* jsonText = R"({
  "title": "beta constant term",
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
      {"nodes": [1, 3], "dof": ["x", "y", "z"]},
      {"nodes": [2], "dof": ["y", "z"]}
    ]
  },
  "constraints": {
    "mpc": [
      {
        "slave": {"node": 2, "dof": "x"},
        "masters": [{"node": 3, "dof": "x", "coeff": 1.0}],
        "beta": 0.05
      }
    ]
  }
})";

    Model model;
    JsonReader reader;
    ASSERT_TRUE(reader.Read(WriteTempJson("beta_driven.json", jsonText), model));
    Solve(model);

    // N3 固定 -> u3=0，故 u2 = 1.0*0 + beta = 0.05
    EXPECT_NEAR(model.nodes[1].Displacement[UX], 0.05, 1e-12);
}


// 从节点外力分摊：力加在 slave N2 上，应经 Tᵀ 传到 master N3
TEST(BetaTest, LoadOnSlaveNode) {
  const char* jsonText = R"({
  "title": "load on slave node",
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
      {"node": 2, "dof": "x", "value": 10.0}
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
  ASSERT_TRUE(reader.Read(WriteTempJson("load_on_slave.json", jsonText), model));
  Solve(model);

  double um = model.nodes[2].Displacement[UX];  // master N3.x
  double us = model.nodes[1].Displacement[UX];  // slave  N2.x
  EXPECT_NEAR(um, 0.1, 1e-9);   // 力经 Tᵀ 传到 master，等效 F/k
  EXPECT_NEAR(us, 0.1, 1e-9);   // slave 回代 = master
}


// 从节点外力分摊：力加在 slave N2 上，应经 Tᵀ 传到 master N3
TEST(BetaTest, LoadOnSlaveNodeForce) {
  const char* jsonText = R"({
  "title": "load on slave node",
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
      {"node": 2, "dof": "x", "value": 10.0}
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
  ASSERT_TRUE(reader.Read(WriteTempJson("load_on_slave.json", jsonText), model));
  Solve(model);

  double um = model.nodes[2].Displacement[UX];  // master N3.x
  double us = model.nodes[1].Displacement[UX];  // slave  N2.x
  EXPECT_NEAR(um, 0.1, 1e-9);   // 力经 Tᵀ 传到 master，等效 F/k
  EXPECT_NEAR(us, 0.1, 1e-9);   // slave 回代 = master
}