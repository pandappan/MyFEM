#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include "Reader.h"
#include "Model.h"
#include "Node.h"

TEST(Reader, PrescibedDisplacementSection) {
    // 设置输入
    const std::string filename = "test_reader_disp_input.dat";
    {
        std::ofstream f(filename);
        f << "Reader test\n";
        f << "3 2 1 1\n";                    // dim=3, numnp=2, numeg=1, modex=1
        f << "1 0 2 1  0.0 0.0 0.0\n";       // Node 1
        f << "2 1 2 2  1.0 0.0 0.0\n";       // Node 2
        f << "2 1 1\n";                      // group
        f << "1 0.0 1000.0 1.0\n";           // material 1
        f << "1 1 2 1\n";                    // element 1
        f << "0\n";                          // NL=0（无集中力）
        f << "3\n";                          // NDISP=3
        f << "1 2  0.01\n";                  // Node 1, Uy, 0.01
        f << "2 2  0.01\n";                  // Node 2, UY, 0.01
        f << "2 3 -0.02\n";                  // Node 2, UZ, -0.02
    }
    // 读取数据
    Model model;
    Reader reader;
    // 验证
    ASSERT_TRUE(reader.Read(filename, model));
    ASSERT_EQ(model.predisplacements.size(), 3u);
    EXPECT_EQ(model.nodes[0].bcode[1], 2);
    EXPECT_EQ(model.nodes[0].Displacement[1], 0.01);
    EXPECT_EQ(model.nodes[1].bcode[1], 2);
    EXPECT_EQ(model.nodes[1].Displacement[1], 0.01);
    EXPECT_EQ(model.nodes[1].bcode[2], 2);
    EXPECT_EQ(model.nodes[1].Displacement[2], -0.02);
    std::remove(filename.c_str());
}