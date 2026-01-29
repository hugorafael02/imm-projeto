import imm_module
import time
import os

# --- CONFIGURAÇÃO ---
# Use r"" antes da string para o Windows não se confundir com as barras
arquivo_grafo = r"C:\Users\hugor\OneDrive\Documentos\Projetos Visual Studio\bases\imm\com-orkut.ungraph-imm.txt"

# Verificação de segurança
if not os.path.exists(arquivo_grafo):
    print(f"ERRO: O arquivo não foi encontrado em: {arquivo_grafo}")
    exit()

print("==================================================")
print(" INICIANDO PONTE PYTHON -> C++ ")
print("==================================================")
print(f"Lendo grafo: {arquivo_grafo}")
print("Configuração: k=5, Modelo=IC, Epsilon=0.5, MC=10000")

inicio = time.time()

# A MÁGICA ACONTECE AQUI
# O Python passa a bola para o C++, que vai usar todos os núcleos do seu processador
imm_module.run_cplusplus(
    arquivo_grafo, 
    5,     # k (sementes)
    "ic",   # modelo (ic ou lt)
    0.5,    # epsilon (precisão)
    10000    # numero de simulações Monte Carlo
)

fim = time.time()

print("==================================================")
print(f" SUCESSO! Tempo Total: {fim - inicio:.4f} segundos")
print("==================================================")