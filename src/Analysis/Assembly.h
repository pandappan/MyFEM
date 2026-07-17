//
// Created by Administrator on 2026/7/6.
//

#pragma once
class Model;
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
};

