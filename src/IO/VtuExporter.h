//
// Created by Administrator on 2026/7/6.
//
#pragma once
#include <string>
class Model;
class VtuExporter {
public:
    bool ExportVtk(const std::string& fileName, const Model& model);
};
