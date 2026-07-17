#include <fstream>
#include <iostream>
#include <ostream>
#include "json.hpp"

using json = nlohmann::json;

int main(int argc, char* argv[]) {
    // 输入检测
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <json-file>" << std::endl;
        return 1;
    }
    // 文件检测
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "Could not open " << argv[1] << std::endl;
        return 1;
    }
    // json输入检测
    json j;
    try {
        in >> j;
    } catch (const json::parse_error& e) {
        std::cerr << "Json parse error: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Tiele: " << j.at("title").get<std::string>() << '\n';
    std::cout << "Dimension: " << j.at("dimension").get<int>() << '\n';
    std::cout << "\n Nodes \n";
    for (const auto& node: j.at("nodes")) {
        int id = node.at("id").get<int>();
        double x = node.at("x").get<double>();
        double y = node.at("y").get<double>();
        double z = node.value("z",0.0);
        std::cout << id << " " << x << " " << y << " " << z << '\n';
    }
    std::cout << "\n Materials \n";
    for (const auto& mat : j.at("materials")) {
        std::cout << "  id=" << mat.at("id").get<int>()
                  << "  type=" << mat.at("type").get<std::string>()
                  << "  E=" << mat.at("E").get<double>()
                  << "  nu=" << mat.at("nu").get<double>()
                  << "  thickness=" << mat.value("thickness", 1.0)
                  << "\n";
    }
    std::cout << "\nElement groups:\n";
    for (const auto& grp : j.at("element_groups")) {
        std::cout << "  type=" << grp.at("type").get<std::string>()
                  << "  material=" << grp.at("material").get<int>()
                  << "  #elements=" << grp.at("elements").size() << "\n";
        for (const auto& elem : grp.at("elements")) {
            std::cout << "    elem id=" << elem.at("id").get<int>()
                      << "  nodes=[";
            for (const auto& n : elem.at("nodes"))
                std::cout << n.get<int>() << " ";
            std::cout << "]\n";
        }
    }
    std::cout << "\nBoundary conditions (fixed):\n";
    for (const auto& bc : j.at("boundary_conditions").at("fixed")) {
        std::cout << "  nodes=[";
        for (const auto& n : bc.at("nodes")) std::cout << n.get<int>() << " ";
        std::cout << "]  dof=[";
        for (const auto& d : bc.at("dof")) std::cout << d.get<std::string>() << " ";
        std::cout << "]\n";
    }
    // 集中载荷
    std::cout << "\nConcentrated loads:\n";
    for (const auto& load : j.at("loads").at("concentrated")) {
        std::cout << "  node=" << load.at("node").get<int>()
                  << "  dof=" << load.at("dof").get<std::string>()
                  << "  value=" << load.at("value").get<double>() << "\n";
    }
    std::cout << "\nSpike OK — nlohmann/json works and schema parses.\n";
    return 0;
}