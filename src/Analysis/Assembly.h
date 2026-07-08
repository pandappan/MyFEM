//
// Created by Administrator on 2026/7/6.
//

#pragma once
class Model;
class Assembler {
public:
    static void CalculateEquationNumber(Model& model);
    static void CalculateLocationMatrix(Model& model);
    static void AllocateStiffnessMatrix(Model& model);
    static void AssembleStiffnessMatrix(Model& model);
    static void AssembleForce(Model& model);
    static void WriteDisplacementToNodes(Model& model);
};

