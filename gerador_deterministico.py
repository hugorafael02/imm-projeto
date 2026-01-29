import numpy as np
import random
import sys
import time

def generate_deterministic_scale_free(
    num_nodes, 
    edges_per_node, 
    seed, 
    output_file, 
    alpha=1, 
    beta_param=5
):
    print(f"--- Configuração Blindada (Cross-Platform) ---")
    print(f"Nós: {num_nodes:,}")
    print(f"Arestas/Nó (m): {edges_per_node}")
    print(f"Semente: {seed}")
    
    # 1. TRAVA AS SEMENTES
    random.seed(seed)
    np.random.seed(seed)

    start_time = time.time()

    # 2. CORREÇÃO CRÍTICA 1: newline='\n'
    # Isso força o Windows a usar quebra de linha do Linux (LF) em vez de (CRLF)
    # Garante que os bytes do arquivo sejam idênticos em qualquer OS.
    with open(output_file, 'w', newline='\n') as f:
        
        # Urna de preferential attachment
        targets = list(range(edges_per_node)) * edges_per_node
        
        buffer = []
        BUFFER_SIZE = 50000 
        
        # Gera arestas iniciais
        for i in range(edges_per_node):
            for j in range(i + 1, edges_per_node):
                w = np.random.beta(alpha, beta_param)
                buffer.append(f"{i} {j} {w:.6f}\n")

        # Loop principal
        for source in range(edges_per_node, num_nodes):
            
            neighbors = set()
            while len(neighbors) < edges_per_node:
                idx = random.randint(0, len(targets) - 1)
                neighbor = targets[idx]
                neighbors.add(neighbor)
            
            # 3. CORREÇÃO CRÍTICA 2: sorted(neighbors)
            # Sets não têm ordem garantida. O Windows pode iterar {1,5} como (1,5)
            # e o Linux como (5,1). Isso muda o hash. O sorted() resolve isso.
            for neighbor in sorted(neighbors):
                weight = np.random.beta(alpha, beta_param)
                buffer.append(f"{source} {neighbor} {weight:.6f}\n")
                
                targets.append(source)
                targets.append(neighbor)
            
            if len(buffer) >= BUFFER_SIZE:
                f.writelines(buffer)
                buffer = []
                
                if source % 100000 == 0:
                    elapsed = time.time() - start_time
                    perc = (source / num_nodes) * 100
                    print(f"\rProgresso: {perc:.1f}% | Nós: {source:,}", end="")

        if buffer:
            f.writelines(buffer)

    print(f"\n\nSucesso! Arquivo gerado: {output_file}")
    print(f"Tempo Total: {time.time() - start_time:.2f}s")

# --- EXECUÇÃO ---
generate_deterministic_scale_free(
    num_nodes=1_000_000,     
    edges_per_node=10,       
    seed=424242,             
    output_file="grafo_cross_platform.txt"
)
