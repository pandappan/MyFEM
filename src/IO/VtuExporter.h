//
// Created by Administrator on 2026/7/6.
//
#pragma once
#include <string>
class Model;
class VtuExporter {
public:
    // 网格文件，节点数据
    bool ExportMesh(const std::string& fileName, const Model& model);
    // 积分点数据
    bool ExportGaussPoints(const std::string& fileName, const Model& model);
};
