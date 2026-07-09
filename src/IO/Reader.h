//
// Created by Administrator on 2026/7/5.
//

#ifndef MYFEM_READER_H
#define MYFEM_READER_H

#include <fstream>
#include <string>
class Model;
class Reader {
public:
    bool Read(const std::string& filename, Model& model);
    static bool TryConsumeKeyword(std::ifstream& input, const std::string& keyword);
private:
    bool ReadHeader(std::ifstream& input, Model& model);
    bool ReadNodes(std::ifstream& input, Model& model);
    bool ReadGroups(std::ifstream& input, Model& model);
    bool ReadLoads(std::ifstream& input, Model& model);
    bool ReadPreDisp(std::ifstream& input, Model& model);
    bool ReadSLoads(std::ifstream& input, Model& model);
    bool ReadBodyForce(std::ifstream& input, Model& model);
};


#endif //MYFEM_READER_H
