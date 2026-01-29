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
    """
    Gera um grafo Scale-Free de forma determinística e eficiente em memória.
    Usa o modelo Barabási-Albert (Preferential Attachment).
    """
    
    print(f"--- Configuração ---")
    print(f"Nós: {num_nodes:,}")
    print(f"Arestas/Nó (m): {edges_per_node}")
    print(f"Semente (Seed): {seed}")
    print(f"Distribuição Beta: alpha={alpha}, beta={beta_param}")
    print(f"--------------------")

    # 1. TRAVA AS SEMENTES (O Segredo do Determinismo)
    # Isso garante que a topologia seja idêntica em qualquer máquina
    random.seed(seed)
    # Isso garante que os pesos sejam idênticos em qualquer máquina
    np.random.seed(seed)

    start_time = time.time()

    with open(output_file, 'w') as f:
        # Opcional: Escrever cabeçalho
        # f.write(f"{num_nodes} {num_nodes * edges_per_node}\n")

        # Lista que representa a urna de "preferential attachment"
        # Nós com mais conexões aparecem mais vezes aqui.
        # Iniciamos com um grafo completo pequeno de m nós
        targets = list(range(edges_per_node)) * edges_per_node
        
        # Buffer para escrita rápida
        buffer = []
        BUFFER_SIZE = 50000 
        
        # Gera arestas iniciais (clique inicial)
        for i in range(edges_per_node):
            for j in range(i + 1, edges_per_node):
                w = np.random.beta(alpha, beta_param)
                buffer.append(f"{i} {j} {w:.6f}\n")

        # Loop principal: Adiciona nós um a um
        # Começa do nó 'm' até 'num_nodes'
        for source in range(edges_per_node, num_nodes):
            
            # Escolhe 'm' vizinhos baseados na probabilidade (grau)
            # A lista 'targets' contém nós repetidos proporcionalmente ao grau
            # Como random.seed está travado, essa escolha é determinística
            neighbors = set()
            while len(neighbors) < edges_per_node:
                # Pegamos um indice aleatorio da lista targets
                # Usamos random.randint para ser rapido e deterministico
                idx = random.randint(0, len(targets) - 1)
                neighbor = targets[idx]
                neighbors.add(neighbor)
            
            # Processa as arestas criadas
            for neighbor in neighbors:
                # Gera peso determinístico
                weight = np.random.beta(alpha, beta_param)
                
                # Adiciona ao buffer
                buffer.append(f"{source} {neighbor} {weight:.6f}\n")
                
                # Atualiza a lista de preferential attachment
                # Adicionamos o novo nó e o vizinho escolhido
                targets.append(source)
                targets.append(neighbor)
            
            # Descarga o buffer periodicamente para não encher a RAM
            if len(buffer) >= BUFFER_SIZE:
                f.writelines(buffer)
                buffer = []
                
                # Log de progresso
                if source % 100000 == 0:
                    elapsed = time.time() - start_time
                    perc = (source / num_nodes) * 100
                    print(f"\rProgresso: {perc:.1f}% | Nós: {source:,} | Tempo: {elapsed:.0f}s", end="")

        # Grava o restante
        if buffer:
            f.writelines(buffer)

    print(f"\n\nSucesso! Arquivo gerado: {output_file}")
    print(f"Tempo Total: {time.time() - start_time:.2f}s")
    
    # Dica de Verificação
    print("DICA: Para verificar se os arquivos são iguais nas duas máquinas,")
    print("rode o comando de hash (md5sum no Linux ou CertUtil no Windows).")

# --- EXECUÇÃO ---
# Gere um grafo pequeno primeiro para validar
# Depois mude para 10_000_000 (10 Milhões)
generate_deterministic_scale_free(
    num_nodes=1_000_000,     # Mude para 10M quando quiser testar a escala
    edges_per_node=10,       # Gera ~10M arestas (ou 100M se nós=10M)
    seed=424242,             # MANTENHA IGUAL NAS DUAS MÁQUINAS
    output_file="grafo_teste_deterministico.txt",
    alpha=0.6, 
    beta_param=5

)