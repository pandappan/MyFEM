#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include "Assembly.h"
#include "Model.h"
#include "Node.h"
#include "json.hpp"
#include "JsonReader.h"
#include "SkylineMatrix.h"

using namespace nlohmann::literals;

namespace {
// 把 JSON 文本写到临时文件，返回完整路径（Read 需要含扩展名的完整文件名）
std::string WriteTempJson(const std::string& name, const std::string& content) {
    std::string path = name;
    std::ofstream ofs(path);
    ofs << content;
    ofs.close();
    return path;
}
}

// -------- 无 MPC：验证 profile 与改动前一致（零回归基线） --------
TEST(ProfileTest, NoMpcBaseline) {

    const char* jsonText =  R"({ "title": "profile baseline without mpc"
  , "dimension": 3,
  "nodes": [
    {"id": 1, "x": 0, "y": 0, "z": 0},
    {"id": 2, "x": 1, "y": 0, "z": 0},
    {"id": 3, "x": 0, "y": 1, "z": 0},
    {"id": 4, "x": 0, "y": 0, "z": 1}
  ],
  "materials": [
    {"type": "bar", "E": 200, "area": 0.05}
  ],
  "element_groups": [
    {"type": "Bar3D", "material": 1,
     "elements": [
       {"id": 1, "connectivity": [1, 2]},
       {"id": 2, "connectivity": [2, 3]},
       {"id": 3, "connectivity": [3, 4]}
     ]}
  ],
  "boundary_conditions": {
    "fixed": [
      {"nodes": [1, 4], "dof": ["x", "y", "z"]}
    ],
    "prescribed": [
      {"node": 2, "dof": "z", "value": 0.1}
    ]
  }
})";

    Model model;
    JsonReader reader;
    std::string path = WriteTempJson("profile_baseline.json", jsonText);
    ASSERT_TRUE(reader.Read(path, model));

    Assembler::CalculateEquationNumber(model);
    Assembler::CalculateLocationMatrix(model);
    Assembler::AllocateLinearSystem(model);

    // 方程编号推演：
    //  N1 固定 -> 无方程
    //  N2 x,y 自由 -> eqn 1,2 ; z 指定位移 -> 无方程
    //  N3 x,y,z 自由 -> eqn 3,4,5
    //  N4 固定 -> 无方程
    EXPECT_EQ(model.neq, 5u);

    // 列高（1..5，0-based 存储）: [0,1,2,3,4]
    auto& H = model.K->GetColumnHeights();
    EXPECT_EQ(H[0], 0u);
    EXPECT_EQ(H[1], 1u);
    EXPECT_EQ(H[2], 2u);
    EXPECT_EQ(H[3], 3u);
    EXPECT_EQ(H[4], 4u);

    // NWK = sum(H) + neq = 10 + 5 ; MK = max(H)+1 = 5
    EXPECT_EQ(model.K->size(), 15u);
    EXPECT_EQ(model.K->GetMaximumHalfBandwidth(), 5u);
}

// -------- 含 MPC：验证远端 master 耦合确实抬高了 profile --------
TEST(ProfileTest, MpcRaisesColumnHeight) {
    // 链式网格 1-2-3-4，仅固定 N1。
    // N4.x 设为从自由度，master = N2.x（编号很靠前）。
    // 这样连接 N3、N4 的 elem3 会把 N4.x 展开成 N2.x(eqn1)，
    // 与 N4 的其它自由度 (eqn7,8) 产生远端耦合。
    const char* jsonText = R"({
  "title": "profile with mpc",
  "dimension": 3,
  "nodes": [
    {"id": 1, "x": 0, "y": 0, "z": 0},
    {"id": 2, "x": 1, "y": 0, "z": 0},
    {"id": 3, "x": 2, "y": 0, "z": 0},
    {"id": 4, "x": 3, "y": 0, "z": 0}
  ],
  "materials": [
    {"type": "bar", "E": 200, "area": 0.05}
  ],
  "element_groups": [
    {"type": "Bar3D", "material": 1,
     "elements": [
       {"id": 1, "connectivity": [1, 2]},
       {"id": 2, "connectivity": [2, 3]},
       {"id": 3, "connectivity": [3, 4]}
     ]}
  ],
  "boundary_conditions": {
    "fixed": [
      {"nodes": [1], "dof": ["x", "y", "z"]}
    ]
  },
  "constraints": {
    "mpc": [
      {
        "slave": {"node": 4, "dof": "x"},
        "masters": [{"node": 2, "dof": "x", "coeff": 1.0}],
        "beta": 0.0
      }
    ]
  }
})";

    Model model;
    JsonReader reader;
    std::string path = WriteTempJson("profile_mpc.json", jsonText);
    ASSERT_TRUE(reader.Read(path, model));

    Assembler::CalculateEquationNumber(model);
    Assembler::CalculateLocationMatrix(model);
    Assembler::AllocateLinearSystem(model);

    // 方程编号推演：
    //  N1 固定 -> 无方程
    //  N2 free  -> eqn 1,2,3
    //  N3 free  -> eqn 4,5,6
    //  N4 x=slave(无方程), y,z free -> eqn 7,8
    EXPECT_EQ(model.neq, 8u);

    // 从自由度 N4.x 应被标为 bcode=3
    EXPECT_EQ(model.nodes[3].bcode[UX], 3u);

    // elem3 展开后的等效方程集合 = {4,5,6, 1(来自master), 7, 8}
    // nfirstrow=1，于是 N4 的 y/z 列 (eqn7,8) 高度被抬到 6、7
    auto& H = model.K->GetColumnHeights();
    EXPECT_EQ(H[6], 6u);   // 列 7：7-1
    EXPECT_EQ(H[7], 7u);   // 列 8：8-1

    // 完整列高 [0,1,2,3,4,5,6,7]
    EXPECT_EQ(H[0], 0u);
    EXPECT_EQ(H[1], 1u);
    EXPECT_EQ(H[2], 2u);
    EXPECT_EQ(H[3], 3u);
    EXPECT_EQ(H[4], 4u);
    EXPECT_EQ(H[5], 5u);

    // NWK = 28 + 8 ; MK = 7+1
    EXPECT_EQ(model.K->size(), 36u);
    EXPECT_EQ(model.K->GetMaximumHalfBandwidth(), 8u);
}