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

void Assembler::AllocateStiffnessMatrix(Model &model) {
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

void Assembler::AssembleStiffnessMatrix(Model &model) {
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            unsigned int nd = element.GetND();
            DenseMatrix<double> ke(nd, nd);
            element.ElementStiffness(ke);
            model.K->Assembly(ke,element.GetLocationMatrix());
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