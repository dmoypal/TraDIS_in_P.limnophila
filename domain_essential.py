from collections import defaultdict

def load_gene_data(file_path):
    gene_positions = {}
    gene_chromosomes = {}
    
    with open(file_path, 'r') as f:
        for line in f:
            name, chromosome, _, _, start, end = line.strip().split('\t')
            gene_positions[name] = [int(start), int(end)]
            gene_chromosomes[name] = chromosome
    
    return gene_positions, gene_chromosomes

def load_insertion_data(file_path, main_chromosome, plasmid_chromosome):
    main_insertions = []
    plasmid_insertions = []
    
    with open(file_path, 'r') as f:
        for line in f:
            data = line.strip().split()
            chromosome = data[0]
            insertion_position = int(data[2])
    
            if chromosome == main_chromosome:
                main_insertions.append(insertion_position)
            elif chromosome == plasmid_chromosome:
                plasmid_insertions.append(insertion_position)
    
    return main_insertions, plasmid_insertions

def generate_insertion_info(gene_positions, gene_chromosomes, main_insertions, plasmid_insertions, window_size=300, increase=150):
    insertion_windows = defaultdict(list)
    
    with open('insertion_in_windows.tsv', 'w') as f:
        for gene, (start, end) in gene_positions.items():
            i = start
            if gene_chromosomes[gene] == 'CP001744':
                search_in = main_insertions
            elif gene_chromosomes[gene] == 'CP001745':
                search_in = plasmid_insertions
    
            while i <= end - window_size:
                insertions = [x for x in search_in if i <= x <= i + window_size]
                index_ratio = round(len(insertions) / window_size, 4)
                f.write(f'{gene}\t{start}-{end}\t{i}-{i + window_size}\t{len(insertions)}\t{index_ratio}\n') # genes with windows with len(insertions) = 0 are reported as domain_essential
                i += increase
    
            f.write(f'{gene}\t{start}-{end}\t-\t-\t-\t-\n')

if __name__ == "__main__":
    gene_positions, gene_chromosomes = load_gene_data('posicion_genes_plm.tsv')
    main_insertions, plasmid_insertions = load_insertion_data('insertions.bed', 'CP001744.1', 'CP001745.1')
    
    generate_insertion_info(gene_positions, gene_chromosomes, main_insertions, plasmid_insertions)
