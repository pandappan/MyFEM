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
    return true;
}

bool Reader::ReadHeader(std::ifstream& input, Model& model) {
    char title[256];
    input.getline(title, 256);
    model.title = title;
    unsigned int numnp, numgp;
    input >> model.dimension >> numnp >> numgp >> model.modex;
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
    model.loads.reserve(NL);
    for (unsigned int i = 0; i < NL; i++) {
        ConcentratedLoad load{};
        input >> load.node >> load.dof >> load.value;
        model.loads.push_back(load);
        if (load.node >= 1 && load.node <= model.nodes.size()) {
            model.nodes[load.node - 1].AddForce((load.dof-1), load.value);
        } else {
            return false;
        }
    }
    return true;
}
