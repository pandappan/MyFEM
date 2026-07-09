//
// Created by Administrator on 2026/7/6.
//

#include "Assembly.h"
#include "../Model/Element/Element.h"
#include "../Model/Model.h"
#include "../Core/DenseMatrix.h"

void Assembler::CalculateEquationNumber(Model& model) {
    model.neq = 0;
    for (auto& node : model.nodes)
        node.GenerateNodeEquation(model.neq);
}

void Assembler::CalculateLocationMatrix(Model &model) {
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            element.GenerateLocationMatrix();
        }
    }
}

void Assembler::AllocateLinearSystem(Model &model) {
    // Allocate for global force/displacement vector
    model.force.assign(model.neq, 0.0);

    // Create the banded stiffness matrix
    model.K.reset(new CSkylineMatrix<double>(model.neq));
    // 计算列高
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            model.K->CalculateColumnHeight(element.GetLocationMatrix());
        }
    }
    model.K->CalculateMaximumHalfBandwidth();
    model.K->CalculateDiagnoalAddress();
    model.K->Allocate();
}

// 循环装配单元刚度矩阵和指定位移约束的造成的右端修正项
void Assembler::AssembleStiffnessAndConstraintCorrection(Model &model) {
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            unsigned int nd = element.GetND();
            const std::vector<unsigned int>& lm = element.GetLocationMatrix();
            // 装配单元刚度矩阵
            DenseMatrix<double> ke(nd, nd);
            element.ElementStiffness(ke);
            model.K->Assembly(ke,lm);
            // 装配右端修正项
            std::vector<double> right(nd);
            element.ElementRight(ke, right);
            for (unsigned int i = 0; i < nd; i++) {
                if (lm[i] != 0) {
                    model.force[lm[i]-1] -= right[i];
                }
            }
        }
    }
}

void Assembler::AssembleForce(Model &model) {
    std::fill(model.force.begin(), model.force.end(), 0.0);
    for (auto& node : model.nodes) {
        for (unsigned int d = 0; d < CNode::NDF; d++) {
            unsigned int eq = node.eqn[d];
            if (eq) {
                model.force[eq - 1] += node.GetForce(d);
            }
        }
    }
}

void Assembler::WriteDisplacementToNodes(Model &model) {
    for (auto& node : model.nodes) {
        node.UpdataNodeDisplacement(model.force);
    }
}

void Assembler::CalculateNodalBCForce(Model &model) {
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            element.CalculateBCForce();
        }
    }
}
