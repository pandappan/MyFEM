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

// 构建单元的全局映射表，用于面载荷寻找所属单元
void Assembler::BuildGlobalElementIndex(Model &model) {
    unsigned int totalElements = 0;
    for (auto& group : model.groups) {
        totalElements += group.GetNUME();
    }
    model.globalElementList.reserve(totalElements);
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            model.globalElementList.push_back(&group.GetElement(e));
        }
    }
}

// 面载荷转化为等效节点力，存入节点中
void Assembler::ConvertSLoadsToCLoads(Model &model) {
    // 循环所有面载荷，逐步转化为等效节点载荷，并存入节点中
    for (auto& sload: model.sloads) {
        unsigned int elemID_0based = sload.elemID - 1;
        unsigned int faceID_0based = sload.faceID - 1; // 转化为0基
        unsigned int dof_0based = sload.dof -1; // 转化为0基
        CElement* elem = model.globalElementList[elemID_0based];
        elem->CalculateSurfaceLoad(faceID_0based, dof_0based, sload.value);
    }
}

// 所有单元体载荷转化为等效节点力，存入节点中
void Assembler::ConvertBLoadsToCLoads(Model &model) {
    // 体力加速度为0，直接退出，避免后续浪费计算量
    const double* b = model.bodyForce;
    double totalForce = 0.0;
    totalForce = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    if (totalForce < 1.0e-12) return;
    // 正常转换体力
    for (auto& group : model.groups) {
        unsigned int nume = group.GetNUME();
        for (unsigned int e = 0; e < nume; e++) {
            CElement& element = group.GetElement(e);
            element.CalculateBodyForce(b);
        }
    }
}

// 装配节点上的所有力，包括点载荷，面载的等效点载，体载的等效点载
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