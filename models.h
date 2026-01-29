#ifndef MODELS_H
#define MODELS_H

#include <vector>
#include <random> // OBRIGATÓRIO para mt19937

using namespace std;

typedef struct {
	int u, v;
	double c; // Mantendo 'c' conforme seu código anterior
} edge;

// Apenas os protótipos do Monte Carlo (já que gen_RR está no imm.cpp)
int MonteCarlo_IC(int V, vector<vector<edge> >& es, vector<int>& S, std::mt19937& gen);
int MonteCarlo_LT(int V, vector<vector<edge> >& es, vector<int>& S, std::mt19937& gen);

#endif