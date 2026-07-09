//
// Created by Administrator on 2026/7/5.
//

#include "Reader.h"
#include "../Model/Model.h"
#include "../Model/Node.h"
#include "../Model/Element/ElementGroup.h"
#include <iostream>

bool Reader::Read(const std::string& filename, Model& model) {
    std::ifstream input(filename.c_str());
    if (!input) {
        std::cerr << "Error while reading file " << filename << std::endl;
        return false;
    }
    if (!ReadHeader(input, model)) {return false;}
    if (!ReadNodes(input, model)) {return false;}
    if (!ReadGroups(input, model)) {return false;}
    if (!ReadLoads(input, model)) {return false;}
    if (!ReadPreDisp(input, model)) {return false;}
    return true;
}

bool Reader::ReadHeader(std::ifstream& input, Model& model) {
    char title[256];
    input.getline(title, 256);
    model.title = title;
    unsigned int numnp, numgp;
    input >> model.dimension >> numnp >> numgp >> model.modex;
    if (numnp == 0) {
        std::cerr << "Error input node numbers" << std::endl;
    }
    if (numgp == 0) {
        std::cerr << "Error input group numbers" << std::endl;
    }
    model.nodes.resize(numnp);
    model.groups.resize(numgp);
    return true;
}

bool Reader::ReadNodes(std::ifstream& input, Model& model) {
    for (unsigned int np = 0; np < model.nodes.size(); np++) {
        CNode& node = model.nodes[np];
        if (!node.Read(input, model.dimension)) return false;
        if (node.NodeNumber != np + 1) {
            std::cerr << "Error Nodes must be inputted in order\n"
            << " Expected: " << np + 1
            << ", Got: " << node.NodeNumber << std::endl;
            return false;
        }
    }
    return true;
}

bool Reader::ReadGroups(std::ifstream &input, Model &model) {
    for (auto& group : model.groups) {
        if (!group.Read(input, model.nodes)) return false;
    }
    return true;
}

bool Reader::ReadLoads(std::ifstream &input, Model &model) {
    unsigned int NL;
    input >> NL;
    model.cloads.reserve(NL);
    for (unsigned int i = 0; i < NL; i++) {
        ConcentratedLoad load{};
        input >> load.node >> load.dof >> load.value;
        model.cloads.push_back(load);
        // 点载荷直接加入节点中
        if (load.node >= 1 && load.node <= model.nodes.size()) {
            model.nodes[load.node - 1].AddForce((load.dof-1), load.value);
        } else {
            return false;
        }
    }
    return true;
}

bool Reader::ReadPreDisp(std::ifstream &input, Model &model) {
    unsigned int NP;
    // 向后兼容
    if (!(input >> NP)) return true;
    model.predisplacements.reserve(NP);
    for (unsigned int i = 0; i < NP; i++) {
        PreDisplacement pd{};
        input >> pd.node >> pd.dof >> pd.value;
        model.predisplacements.push_back(pd);
        if (pd.node >= 1 && pd.node <= model.nodes.size()) {
            if (!model.nodes[pd.node - 1].SetPreDisp(pd.dof-1, pd.value)) {
                std::cerr << "Error while setting PreDisplacement for node " << pd.node
                << ", dof" << pd.dof-1
                << ", origin node constrain input is "<< model.nodes[pd.node - 1].bcode[pd.dof-1]
                << std::endl;
                return false;
            }
        } else {
            return false;
        }
    }
    return true;
}

// 读取，存储是1基，后续传入函数中是0基
bool Reader::ReadSLoads(std::ifstream &input, Model &model) {
    unsigned int NS;
    if (!(input >> NS)) return true;
    model.sloads.reserve(NS);
    for (unsigned int i = 0; i < NS; i++) {
        SurfaceLoad sload{};
        input >> sload.elemID >> sload.dof >> sload.value;
        model.sloads.push_back(sload);
    }
    return true;
}
