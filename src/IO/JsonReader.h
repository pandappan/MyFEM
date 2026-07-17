#pragma once
#include <string>
#include "json.hpp"

class Model;
class JsonReader {
public:
    bool Read(const std::string& filename, Model& model);
private:
    using json = nlohmann::json;
    bool ParseHeader(const json& j, Model& model);
    bool ParseNodes(const json& j, Model& model);
    bool ParseMaterials(const json& j, Model& model);
    bool ParseElementGroups(const json& j, Model& model);
    bool ParseBoundaryConditions(const json& j, Model& model);
    bool ParseLoads(const json& j, Model& model);
    static std::vector<unsigned int> DofStringToInt(const std::vector<std::string>& dofsString);
    static unsigned int DofStringToInt(const std::string& dofsString);
    static unsigned int Base1ToBase0(unsigned int idx_1);
    static std::vector<unsigned int> Base1ToBase0(const std::vector<unsigned int> &idxs_1);
};
