//
// Created by Administrator on 2026/7/6.
//

#include <iostream>
#include "Analysis/Assembly.h"
#include "Analysis/Solver.h"
#include "IO/JsonReader.h"
#include "IO/VtuExporter.h"
#include "IO/Writer.h"
#include "Model/Model.h"
#include "Element/CContinuumElement.h"

int main(int argc, char* argv[]) {
    // 输入和输出文件名，不加尾缀
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " InputFile\n";
        return 1;
    }
    std::string inFile  = argv[1];
    std::string outFile = inFile + ".out";
    std::string vtkFile = inFile + ".vtu";
    std::string vtkGsFile = inFile + ".gs.vtu";
    inFile = inFile + ".json";
    // 模型信息
    Model model;
    JsonReader jsonReader;
    if (!jsonReader.Read(inFile,model)) return 1;
    // 输出模型信息
    Writer writer(outFile);
    writer.OutputHeading(model);
    writer.OutputNodeInfo(model);
    // 初始化刚度矩阵与残差
    Assembler::CalculateEquationNumber(model);
    writer.OutputEquationNumber(model);
    writer.OutputElementInfo(model);
    Assembler::CalculateLocationMatrix(model);
    Assembler::AllocateLinearSystem(model);
    writer.OutputTotalSystemData(model);
    if (model.modex == 0) {
        std::cout << "Data check model. Exit! \n";
        return 0;
    }
    Assembler::InitializeElementMap(model);
    Assembler::ConvertSLoadsToCLoads(model);
    Assembler::ConvertBLoadsToCLoads(model);
    // 先装配外载荷，再装配刚度矩阵和右端修正项
    Assembler::AssembleForce(model);
    Assembler::AssembleStiffnessAndConstraintCorrection(model);
    CLDLTSolver solver(*model.K);
    solver.LDLT();
    solver.BackSubstitution(model.force);
    // 求解结果处理
    Assembler::WriteDisplacementToNodes(model);
    Assembler::CalculateNodalBCForce(model);
    Assembler::CalculateNodalStress(model);
    // 可视化导出结果
    VtuExporter exporter;
    exporter.ExportMesh(vtkFile,model);
    exporter.ExportGaussPoints(vtkGsFile,model);
    return 0;
}