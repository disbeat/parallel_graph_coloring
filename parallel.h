/*
 * parallel.h
 *
 * Created by msimoes and naml
 *
 */

#define TAG 42

void plassman_p_processors(int pid, int p, MPI_Comm * comm_workers);
void block_partition(int pid, int p, MPI_Comm * comm_workers);
void greedy_parallel(int pid);
void plassman_v_processors(int pid);

int find_my_color(PblSet * psetcolors, int nodes);
void print_adjacency_matrix(int pid, int * madjacency, int my_nodes, int nodes) ;