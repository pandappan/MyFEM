//
// Created by Administrator on 2026/7/6.
//

#include "Model/Model.h"
#include "IO/Reader.h"
#include "IO/Writer.h"
#include "IO/VtuExporter.h"
#include "Analysis/Assembly.h"
#include "Analysis/Solver.h"
#include <iostream>

int main(int argc, char* argv[]) {
    // 输入和输出文件名，不加尾缀
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " InputFile\n";
        return 1;
    }
    std::string inFile  = argv[1];
    std::string outFile = inFile + ".out";
    std::string vtkFile = inFile + ".vtu";
    inFile = inFile + ".dat";
    // 模型信息
    Model model;
    Reader reader;
    if (!reader.Read(inFile, model)) return 1;
    // 输出模型信息
    Writer writer(outFile);
    writer.OutputHeading(model);
    writer.OutputNodeInfo(model);
    // 初始化刚度矩阵与残差
    Assembler::CalculateEquationNumber(model);
    writer.OutputEquationNumber(model);
    writer.OutputElementInfo(model);
    Assembler::CalculateLocationMatrix(model);
    writer.OutputNodeForce(model);
    Assembler::AllocateLinearSystem(model);
    writer.OutputTotalSystemData(model);
    if (model.modex == 0) {
        std::cout << "Data check model. Exit! \n";
        return 0;
    }
    // 构建单元全局映射表，转化面载荷至节点载荷
    Assembler::BuildGlobalElementIndex(model);
    Assembler::ConvertSLoadsToCLoads(model);
    // 先装配外载荷，再装配刚度矩阵和右端修正项
    Assembler::AssembleForce(model);
    Assembler::AssembleStiffnessAndConstraintCorrection(model);
    CLDLTSolver solver(*model.K);
    solver.LDLT();
    solver.BackSubstitution(model.force);
    Assembler::WriteDisplacementToNodes(model);
    Assembler::CalculateNodalBCForce(model);
    // 输出结果
    writer.OutputNodalDisplacement(model);
    writer.OutputElementStress(model);
    writer.OutputNodalBCForce(model);
    // 可视化导出结果
    VtuExporter exporter;
    exporter.ExportVtk(vtkFile,model);
    return 0;
}