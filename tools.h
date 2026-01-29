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

// --- INICIO DA CORRECAO WINDOWS ---
#ifdef _WIN32
    // Se for Windows, usa bibliotecas nativas
#include <ctime>
#include <process.h> 
#define getpid _getpid // Mapeia a função getpid do Linux para o Windows
#else
    // Se for Linux (seu futuro SageMaker), mantém as originais
#include <sys/time.h>
#include <unistd.h>
#endif
// --- FIM DA CORRECAO ---#include <memory> // auto_ptr

#include "./mt19937ar.h"

using namespace std;

void init_args(int argc, char *argv[], map<string, string> &args);
string get_or_die(map<string, string> &argv, string key);
int genrand_int(int n);
