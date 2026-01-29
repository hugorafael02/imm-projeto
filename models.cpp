#include "./models.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <bitset>
#include <algorithm>
#include <map>
#include <set>
#include <ctime>
#include <cmath>
#include <random>

#include <assert.h>
#include <unordered_set> // <--- Adicione isso no topo
#include <unordered_map> // <--- Adicione no topo se não tiver

#include "./mt19937ar.h"
#include "./tools.h"

using namespace std;

// ============================================================================
// MONTE CARLO IC - VERSÃO ESTÁVEL (SAFE)
// ============================================================================
int MonteCarlo_IC(int V, vector<vector<edge> >& es, vector<int>& S, mt19937& gen) {

	queue<int> Q;
	int count = 0;
	uniform_real_distribution<> dis_prob(0.0, 1.0);

	// ALOCAÇÃO SEGURA (LOCAL)
	// Para V=75k, isso é minúsculo (9KB). Não vai gargalar seu PC.
	// O sistema operacional gerencia isso automaticamente, sem risco de "sujeira" de execuções anteriores.
	vector<bool> active(V, false);

	for (int s : S) {
		// Proteção extra: ignora sementes inválidas
		if (s < V && !active[s]) {
			active[s] = true;
			Q.push(s);
			count++;
		}
	}

	while (!Q.empty()) {
		int u = Q.front();
		Q.pop();

		for (auto& e : es[u]) {
			int v = e.v;

			// Proteção de limites (Bound Check implícito pela lógica, mas garantido aqui)
			if (v < V && !active[v]) {
				if (dis_prob(gen) <= e.c) {
					active[v] = true;
					Q.push(v);
					count++;
				}
			}
		}
	}
	return count;
}

// ============================================================================
// MONTE CARLO LT - VERSÃO ESTÁVEL (SAFE)
// ============================================================================
int MonteCarlo_LT(int V, vector<vector<edge> >& es, vector<int>& S, mt19937& gen) {

	queue<int> Q;
	int count = 0;
	uniform_real_distribution<> dis_prob(0.0, 1.0);

	// ALOCAÇÃO SEGURA
	vector<bool> active(V, false);
	vector<double> incoming_weights(V, 0.0);
	vector<double> thresholds(V, -1.0); // -1 indica não gerado

	for (int s : S) {
		if (s < V && !active[s]) {
			active[s] = true;
			Q.push(s);
			count++;
		}
	}

	while (!Q.empty()) {
		int u = Q.front();
		Q.pop();

		for (auto& e : es[u]) {
			int v = e.v;

			if (v >= V) continue; // Proteção contra nós fora do limite
			if (active[v]) continue;

			// Lazy Threshold
			if (thresholds[v] < 0) {
				thresholds[v] = dis_prob(gen);
			}

			incoming_weights[v] += e.c;

			if (incoming_weights[v] >= thresholds[v]) {
				active[v] = true;
				Q.push(v);
				count++;
			}
		}
	}
	return count;
}