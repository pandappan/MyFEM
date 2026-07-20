//
// Created by Administrator on 2026/7/6.
//

#pragma once
#include <vector>
class CElement;
class Model;
struct LocalDofExpansion {
    unsigned int globalEqn; // 1基全局方程号
    double coeff;
};
// 局部自由度展开为全局方程项目+常数
// 对应变换 u_i = Σ_p coeff_p · a_p + constant
struct DofExpansion {
    std::vector<LocalDofExpansion> terms;
    double constant = 0.0;
};
class Assembler {
public:
    static void CalculateEquationNumber(Model& model);
    static void CalculateLocationMatrix(Model& model);
    static void AllocateLinearSystem(Model& model);
    static void InitializeElementMap(Model& model);
    static void ConvertSLoadsToCLoads(Model& model);
    static void ConvertBLoadsToCLoads(Model& model);
    static void AssembleForce(Model& model);
    static void AssembleStiffnessAndConstraintCorrection(Model& model);
    static void WriteDisplacementToNodes(Model& model);
    static void CalculateNodalBCForce(Model& model);
    static void CalculateNodalStress(Model& model);
    // 从约束方程中计算从自由度的约束位移
    static void RecoverSlaveDisplacement(Model& model);
    // 单元展开后实际涉及的所有全局方程号（用于 skyline 列高）
    static std::vector<unsigned int> GetEffectiveEquations(const CElement& element,
                                                           const Model& model);
    static DofExpansion GetLocalDofExpansion(const CElement& element,
        unsigned int localDof, const Model& model);
};

