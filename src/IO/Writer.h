//
// Created by Administrator on 2026/7/6.
//
#pragma once
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// 作用：将处理后的输入信息打印，检查前处理是否存在问题
class Model;
class Writer {
private:
    std::ofstream outputFile_;
    void PrintTime();
    template<class T>
    void Tee(T fn) {
        fn(outputFile_);
        fn(std::cout);
    }
public:
    explicit Writer(const std::string& fileName);
    ~Writer() = default;
    void OutputHeading(const Model& model);
    void OutputNodeInfo(const Model& model);
    void OutputEquationNumber(const Model& model);
    void OutputElementInfo(const Model& model);
    void OutputTotalSystemData(const Model& model);
    // 操作符重载
    template<class T>
    Writer& operator<<(const T& value) {
        std::cout << value;
        outputFile_ << value;
        return *this;
    }
    // 支持endl等符号
    using CharOstream = std::basic_ostream<char>;
    Writer& operator<<(CharOstream& (*op)(CharOstream&)) {
        op(std::cout);
        op(outputFile_);
        return *this;
    }
};

