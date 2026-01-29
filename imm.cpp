#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <limits>
#include <iomanip>
#include <unordered_set> // <--- Adicione isso no topo
#include <utility> // Necessário para std::move

#include <omp.h> // <--- Adicione esta linha nos includes

#include "./models.h"
// #include "./mt19937ar.h"
#include "./tools.h"

using namespace std;

vector<int> imm(int V, vector<vector<edge> > &rs, 
		string model, int k, double eps, double ell);
int greedy(int V, vector<vector<int> > &h2v, vector<vector<int> > &v2h, int k,
		vector<int> &S);

double log_fact(int n) {
	double val = 1;
	for (int i = 1; i <= n; i++) {
		val += log(i);
	}
	return val;
}

// {n \choose k} = n!/k!(n-k)!
double log_nCk(int n, int k) {
	return log_fact(n) - log_fact(k) - log_fact(n - k);
}

vector<int> gen_RR_IC(int V, vector<vector<edge> >& rs, mt19937& gen) {
	uniform_int_distribution<> dis_node(0, V - 1);
	uniform_real_distribution<> dis_prob(0.0, 1.0);
	
	vector<int> RR;
	queue<int> Q;

	// USAMOS UM SET PARA NÃO ALOCAR VETOR GIGANTE DE 4MB A CADA ITERAÇÃO
	unordered_set<int> visited;

	int z = dis_node(gen);

	visited.insert(z); // Marca visitado
	RR.push_back(z);
	Q.push(z);

	while (!Q.empty()) {
		int u = Q.front();
		Q.pop();

		for (auto& e : rs[u]) {
			int v = e.u;

			// Se NÃO está no set (count == 0), então não foi visitado
			if (visited.count(v) == 0) {
				if (dis_prob(gen) <= e.c) {
					visited.insert(v);
					RR.push_back(v);
					Q.push(v);
				}
			}
		}
	}
	return RR;
}

// GERAÇÃO DE RR-SETS PARA LINEAR THRESHOLD (LT)
// Lógica: Random Walk Reverso (Caminhada Aleatória)
vector<int> gen_RR_LT(int V, vector<vector<edge> >& rs, mt19937& gen) {
	uniform_int_distribution<> dis_node(0, V - 1);
	uniform_real_distribution<> dis_prob(0.0, 1.0);

	vector<int> RR;

	// Otimização: Set para não alocar memória O(V)
	unordered_set<int> visited;

	// 1. Seleciona a Raiz
	int u = dis_node(gen);
	RR.push_back(u);
	visited.insert(u);

	// 2. Loop de Caminhada (Random Walk)
	// Diferente do IC, não precisamos de fila (Queue), pois o caminho é linear.
	// Só seguimos UM antecessor por vez.
	while (true) {

		// Se não há vizinhos de entrada, a caminhada morre
		if (rs[u].empty()) {
			break;
		}

		// 3. Roleta Russa (Weighted Selection)
		// Sorteamos um número entre 0 e 1
		double val = dis_prob(gen);
		double cumulative = 0.0;
		int next_node = -1;

		// Percorre os vizinhos acumulando os pesos (e.p)
		for (auto& e : rs[u]) {
			cumulative += e.c;
			if (val <= cumulative) {
				next_node = e.u; // Escolheu este vizinho
				break;
			}
		}

		// 4. Transição
		// Se selecionamos alguém E ele ainda não foi visitado (evita ciclos)
		if (next_node != -1 && visited.count(next_node) == 0) {
			visited.insert(next_node);
			RR.push_back(next_node);
			u = next_node; // Avança para o próximo nó
		}
		else {
			// Se a roleta caiu numa faixa "vazia" (val > soma dos pesos)
			// OU se encontramos um ciclo -> Paramos.
			break;
		}
	}

	return RR;
}

// --- GREEDY OTIMIZADO (CELF / LAZY) ---
// Substitui a versão ingênua que recalculava tudo
int greedy(int V, vector<vector<int> >& h2v, vector<vector<int> >& v2h, int k, vector<int>& S) {
	int H = (int)h2v.size();
	vector<bool> dead(H, false);
	vector<int> deg(V);

	// Fila de Prioridade para CELF
	priority_queue<pair<int, int> > Q;

	// 1. Inicialização
	for (int v = 0; v < V; v++) {
		deg[v] = (int)v2h[v].size();
		if (deg[v] > 0) {
			Q.push(make_pair(deg[v], v));
		}
	}

	int total_covered = 0;
	vector<bool> selected(V, false);

	// 2. Loop Lazy
	while (S.size() < k && !Q.empty()) {
		pair<int, int> top = Q.top();
		Q.pop();

		int v = top.second;
		int stored_deg = top.first;

		if (selected[v]) continue;

		// Se o grau armazenado é igual ao atual, é o vencedor real
		if (stored_deg == deg[v]) {
			S.push_back(v);
			selected[v] = true;
			total_covered += deg[v];

			// Marca os RR-Sets como cobertos e penaliza vizinhos
			for (int h_idx : v2h[v]) {
				if (!dead[h_idx]) {
					dead[h_idx] = true;
					for (int u : h2v[h_idx]) {
						deg[u]--;
					}
				}
			}
		}
		else {
			// Valor estava desatualizado, re-insere com valor novo
			Q.push(make_pair(deg[v], v));
		}
	}
	return total_covered;
}

// ... (outros includes permanecem os mesmos)

// =================================================================================
// FUNÇÃO AUXILIAR: GERAÇÃO E INDEXAÇÃO DE RR-SETS
// Faz o trabalho pesado: Gera 'needed' amostras, salva em h2v e atualiza v2h
// =================================================================================
void generate_samples(int V, vector<vector<edge> >& rs, string model, long long needed,
	vector<vector<int> >& h2v, vector<vector<int> >& v2h,
	long long& totW, int seed_salt) {

	if (needed <= 0) return;

	int initial_H_size = h2v.size();

	// 1. GERAÇÃO PARALELA (Com buffer local e std::move)
#pragma omp parallel
	{
		vector<vector<int> > local_h2v;

		unsigned long seed = (unsigned long)(time(NULL) ^ (omp_get_thread_num() * 12345 + seed_salt));
		mt19937 gen(seed);

#pragma omp for nowait
		for (long long j = 0; j < needed; j++) {
			vector<int> RR;
			if (model == "tvlt") {
				RR = gen_RR_LT(V, rs, gen);
			}
			else {
				// Fallback
				RR = gen_RR_IC(V, rs, gen);
			}
			local_h2v.push_back(RR);
		}

		// Merge Global Otimizado
#pragma omp critical
		{
			for (auto& rr : local_h2v) {
				h2v.push_back(std::move(rr));
			}
		}
	}

	// 2. INDEXAÇÃO SEQUENCIAL (v2h)
	// Processa apenas o que foi adicionado agora (do initial_H_size em diante)
	int final_H_size = h2v.size();
	for (int idx = initial_H_size; idx < final_H_size; idx++) {
		for (int v : h2v[idx]) {
			v2h[v].push_back(idx);
			totW++;
		}
	}
}

// =================================================================================
// FUNÇÃO PRINCIPAL IMM
// =================================================================================
vector<int> imm(int V, vector<vector<edge> >& rs, string model, int k, double eps, double ell) {
	const double e = exp(1);
	double log_VCk = log_nCk(V, k);

	ell = ell * (1 + log(2) / log(V));
	double eps_p = sqrt(2) * eps;

	printf("ell  = %f\n", ell);
	printf("eps' = %f\n", eps_p);
	printf("log{V c k} = %f\n", log_VCk);

	double OPT_lb = 1;

	int H = 0;
	vector<vector<int> > h2v;
	vector<vector<int> > v2h(V);
	long long int totW = 0;

	// --- FASE 1: Estimação ---
	for (int i = 1; i <= log2(V) - 1; i++) {
		double x = V / pow(2, i);
		double lambda_prime = (2 + 2.0 / 3.0 * eps_p)
			* (log_VCk + ell * log(V) + log(log2(V))) * V / (eps_p * eps_p);
		double theta_i = lambda_prime / x;

		printf("i = %d\n", i);
		printf("x  = %.0f\n", x);
		printf("theta_i = %.0f\n", theta_i);

		long long iterations_needed = (long long)(theta_i - H);

		// CHAMADA DA FUNÇÃO AUXILIAR
		generate_samples(V, rs, model, iterations_needed, h2v, v2h, totW, i);

		H = h2v.size(); // Atualiza H real

		printf("H  = %d\n", H);
		printf("totW = %lld\n", totW);

		vector<int> S;
		int degS = greedy(V, h2v, v2h, k, S);
		printf("deg(S) = %d\n", degS);
		printf("Inf(S) = %f\n", 1.0 * V * degS / H);
		printf("\n");

		if (1.0 * V * degS / theta_i >= (1 + eps_p) * x) {
			OPT_lb = (1.0 * V * degS) / ((1 + eps_p) * theta_i);
			break;
		}
	}

	// --- CÁLCULO DE PARÂMETROS ---
	double lambda_star;
	{
		double alpha = sqrt(ell * log(V) + log(2));
		double beta = sqrt((1 - 1 / e) * (log_VCk + ell * log(V) + log(2)));
		double c = (1 - 1 / e) * alpha + beta;
		lambda_star = 2 * V * c * c / (eps * eps);
	}
	double theta = lambda_star / OPT_lb;
	printf("OPT_ = %.0f\n", OPT_lb);
	printf("lambda* = %.0f\n", lambda_star);
	printf("theta = %.0f\n", theta);

	// --- FASE 2: Refinamento ---
	long long iterations_needed_2 = (long long)(theta - H);

	// CHAMADA DA FUNÇÃO AUXILIAR
	// Passamos um salt fixo grande (ex: 67890) para diferenciar da Fase 1
	generate_samples(V, rs, model, iterations_needed_2, h2v, v2h, totW, 67890);

	H = h2v.size();
	printf("H  = %d\n", H);

	vector<int> S;
	int degS = greedy(V, h2v, v2h, k, S);
	printf("deg(S) = %d\n", degS);
	printf("Inf(S) = %f\n", 1.0 * V * degS / H);

	return S;
}

void run(map<string, string> args) {
	string input = get_or_die(args, "graph");
	int k = atoi(get_or_die(args, "k").c_str());
	double eps = atof(get_or_die(args, "eps").c_str());
	double ell = atof(get_or_die(args, "ell").c_str());
	string model = get_or_die(args, "model");

	int numMC = atoi(get_or_die(args, "numMC").c_str());

	cout << "[1/4] Lendo arquivo..." << endl;
	ifstream is(input.c_str());
	if (!is.is_open()) {
		cerr << "Erro fatal: Nao foi possivel abrir " << input << endl;
		exit(1);
	}

	vector<edge> ps;
	int V = 0;
	int u, v;
	double p_val;

	while (is >> u >> v >> p_val) {
		if (u == v) continue;
		edge e = { u, v, p_val };
		V = max(V, max(u, v) + 1);
		ps.push_back(e);
	}
	is.close();

	// --- EXIBIÇÃO DAS ESTATÍSTICAS ---
	cout << "========================================" << endl;
	cout << "       IMM - CONFIGURACAO ATUAL" << endl;
	cout << "========================================" << endl;
	cout << "Dataset      : " << input << endl;
	cout << "Nos (V)      : " << V << endl;
	cout << "Arestas (E)  : " << ps.size() << endl;
	cout << "----------------------------------------" << endl;
	cout << "Modelo       : " << (model == "tvlt" ? "Linear Threshold (LT)" : "Independent Cascade (IC)") << endl;
	cout << "Sementes (k) : " << k << endl;
	cout << "Monte Carlo  : " << numMC << " simulacoes" << endl;
	cout << "========================================" << endl << endl;

	cout << "[2/4] Construindo grafo (V=" << V << ", E=" << ps.size() << ")..." << endl;
	vector<vector<edge> > rs(V), es(V);
	for (auto e : ps) {
		rs[e.v].push_back(e);
		es[e.u].push_back(e);
	}

	{
		vector<edge> empty;
		ps.swap(empty);
	}

	cout << "[3/4] Executando IMM (Selecao de Sementes)..." << endl;
	clock_t start_imm = clock(); // Cronômetro IMM
	vector<int> S = imm(V, rs, model, k, eps, ell);
	clock_t end_imm = clock();   // Fim IMM

	cout << "\n>>> TEMPO IMM: " << (double)(end_imm - start_imm) / CLOCKS_PER_SEC << "s <<<" << endl;
	cout << "Sementes Escolhidas: { ";
	for (size_t i = 0; i < S.size(); i++) cout << S[i] << (i < S.size() - 1 ? ", " : "");
	cout << " }" << endl;

	// --- MONTE CARLO COM ESTIMATIVA EM TEMPO REAL ---
	cout << "\n[4/4] Validacao Monte Carlo (" << numMC << " simulacoes)..." << endl;

	clock_t start_mc = clock();

	double total_inf = 0.0;     // Soma final segura
	double global_running_spread = 0.0; // Soma para visualização (Display)
	int progress_counter = 0;

#pragma omp parallel
	{
		unsigned long seed = (unsigned long)(time(NULL) ^ (omp_get_thread_num() * 99999));
		mt19937 gen(seed);

		double local_inf = 0;

#pragma omp for
		for (int sim = 0; sim < numMC; sim++) {
			int result = 0;
			if (model == "tvic" || model == "ic") {
				result = MonteCarlo_IC(V, es, S, gen);
			}
			else if (model == "tvlt" || model == "lt") {
				result = MonteCarlo_LT(V, es, S, gen);
			}

			local_inf += result;

			// Atualiza a soma global para exibição (Atômico é rápido o suficiente aqui)
#pragma omp atomic
			global_running_spread += result;

			// Barra de progresso
#pragma omp atomic
			progress_counter++;

			// Atualiza a cada 10 iterações (já que está levando ~1min para cada 1000)
			if (progress_counter % 100 == 0 || progress_counter == numMC) {
#pragma omp critical
				{
					double current_avg = global_running_spread / progress_counter;

					cout << "\rProgresso: " << progress_counter << "/" << numMC
						<< " (" << (int)(100.0 * progress_counter / numMC) << "%) "
						<< "| Spread Est.: " << fixed << setprecision(2) << current_avg << "   " << flush;
				}
			}
		}

#pragma omp atomic
		total_inf += local_inf;
	}

	clock_t end_mc = clock();

	cout << "\n\nCalculo Finalizado." << endl;
	double final_spread = total_inf / numMC;

	cout << "========================================" << endl;
	cout << "SPREAD MEDIO FINAL: " << final_spread << endl;
	cout << "TEMPO MONTE CARLO : " << (double)(end_mc - start_mc) / CLOCKS_PER_SEC << "s" << endl;
	cout << "========================================" << endl;

	cout << "Pressione ENTER para sair..." << endl;
	cin.get();
}

 //... (todo o código do imm.cpp e função run acima) ...
//int main() {
//	// 1. CAMINHO DO ARQUIVO REAL
//	// Usamos o Raw String Literal R"(...)" para o caminho que você definiu
//	string graph_file = R"(C:\Users\hugor\OneDrive\Documentos\Projetos Visual Studio\bases\imm\com-orkut.ungraph-imm.txt)";
//
//	// 2. PARÂMETROS DA EXECUÇÃO
//	int k_seeds = 3;  // Vamos buscar os 50 nós mais influentes (tamanho padrão para papers)
//	int sim_mc = 10000; // Número de simulações de Monte Carlo para validação (1000 é um bom balanço)
//
//	// 3. PREPARAR ARGUMENTOS
//	map<string, string> args;
//	args["graph"] = graph_file;
//	args["k"] = to_string(k_seeds);
//
//	args["model"] = "tvic";          // Placeholder (nossa lógica ignora o tempo)
//	args["eps"] = "0.3";             // Epsilon 0.1 é o padrão de mercado para precisão
//	args["ell"] = "1.0";
//	
//	args["numMC"] = to_string(sim_mc);
//
//	// Executa o algoritmo
//	// A função run() vai chamar nossa leitura de 3 colunas e processar tudo.
//	run(args);
//
//	cout << "\n========================================" << endl;
//	cout << "         FIM DA EXECUCAO" << endl;
//	cout << "========================================" << endl;
//
//	system("pause");
//	return 0;
//}