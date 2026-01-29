#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 
#include <map>
#include <string>
#include <iostream>

// Inclua o seu código fonte principal
// IMPORTANTE: Comente a função 'main' dentro do imm.cpp antes de compilar!
#include "imm.cpp" 

namespace py = pybind11;

// Wrapper simples
void run_algorithm(std::string graph_path, int k, std::string model, double eps, int numMC) {
    std::map<std::string, std::string> args;
    args["graph"] = graph_path;
    args["k"] = std::to_string(k);
    args["model"] = model;
    args["eps"] = std::to_string(eps);
    args["ell"] = "1.0"; 
    args["numMC"] = std::to_string(numMC);

    run(args); 
}

PYBIND11_MODULE(imm_module, m) {
    m.doc() = "Modulo C++ Otimizado (Windows/Linux)";
    
    m.def("run_cplusplus", &run_algorithm, 
          "Executa o IMM + Monte Carlo",
          py::arg("graph_path"), 
          py::arg("k"), 
          py::arg("model"), 
          py::arg("eps"), 
          py::arg("numMC"));
}