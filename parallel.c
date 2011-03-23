/*
 * parallel.c
 *
 * Created by msimoes and naml
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pbl.h>
#include <mpi.h>
#include <math.h>
#include "common.h"
#include "parallel.h"

#define MIN(x,y) (x<y?x:y)

int get_nodes_per_process(int nodes, int max_processes) {
	return floor( nodes * 1.0 / max_processes );
}

int get_pid_of_node(int node, int nodes, int max_processes) {
	return MIN(floor( node * 1.0 / get_nodes_per_process(nodes, max_processes ) ), max_processes-1) + 1;
}



int find_my_color(PblSet * psetcolors, int nodes) {
	int color;
	for (color = 1; pblSetContains(psetcolors, (void *) color) && color < nodes; color++);
	return color;
}

/*
 * Generate a random value between a and b inclusive. It suposes the srand function was already called.
 */
int rand_between(int low, int high) {

	return rand() % (high - low + 1) + low;
}


int color_node(int node, int pid, int nodes_per_block, int * madjacency, int* colors, PblSet * psetcolors, int nodes) {
	int othernode;
	pblSetClear(psetcolors);
	for (othernode = 0; othernode < nodes; othernode++) {
		if (othernode != (pid-1)*nodes_per_block + node) {
			if (madjacency[POS(othernode, node, nodes)] == 1) {
				if (colors[othernode] != 0)
				/* the neighbor already has a color */
					pblSetAdd(psetcolors, (void *) colors[othernode]);
			}
		}
	}

	return find_my_color(psetcolors, nodes);
}


int is_local(int node, int i, int * madjacency, int nodes, int p) {
	int othernode;
	int pid = get_pid_of_node(node, nodes, p);
	for (othernode = 0; othernode < nodes; othernode++) {
		if (othernode != node && madjacency[POS(othernode, i, nodes)] == 1) {
			if (get_pid_of_node(othernode, nodes, p) != pid)
				return 0;
		}
	}
	return 1;
}


void print_adjacency_matrix(int pid, int * madjacency, int my_nodes, int nodes) {
	int i, j;
	printf("[rank %d] adjacency %d:\n", pid, my_nodes*nodes);
	for ( i = 0; i < my_nodes; i++ ) {
		for ( j = 0; j < nodes; j++ )	
			printf( "%d ", madjacency[i*nodes + j]);
		printf( "[rank %d]\n",pid);
	}
	
	fflush(stdout);
}

void block_partition(int pid, int p, MPI_Comm * comm_workers) {
	
	int nodes, color;
	int my_nodes;
	int max_nodes_per_block, nodes_per_block;
	int * madjacency, *neighbours_colors, *my_colors, *conflicts;
	int n_conflicts = 0;

	int i, j, b, last_color = 0, neighbour_color = -1;

	MPI_Request req;
	MPI_Status status;
	char buffer[10000];
	PblSet * psetcolors = (PblSet *) pblSetNewHashSet();

	MPI_Recv( &nodes, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status );
	
	//printf("[rank %d] nodes: %d\n", pid, nodes);

	nodes_per_block = get_nodes_per_process(nodes, p);
	max_nodes_per_block = nodes - (pid - 1) * nodes_per_block;
	
	if (pid < p)
		my_nodes = nodes_per_block;
	else
		my_nodes = max_nodes_per_block;
		

	//printf("[rank %d] my_nodes: %d | max_nodes: %d | nodes = %d\n", pid, my_nodes, max_nodes_per_block, nodes);
	// alloc space for the neighbours information

	madjacency = (int*) malloc( nodes * my_nodes * sizeof(int) );

	neighbours_colors = (int*) malloc( nodes * sizeof(int) );

	my_colors = (int*) malloc( my_nodes * sizeof(int) );

	conflicts = (int*) malloc( nodes_per_block * sizeof(int) );


	for ( i = 0; i < nodes; i++ )
		neighbours_colors[i] = 0;
	
	MPI_Recv( madjacency, nodes*my_nodes, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status );

	
//	print_adjacency_matrix(pid, madjacency, my_nodes, nodes);

	last_color = 0;
	for ( i = 0; i < nodes_per_block; i++ ) {

		pblSetClear( psetcolors );

		for ( j = 0; j < nodes; j++ ) {
			if (j == (pid-1)*nodes_per_block + i )
				continue;


			if (madjacency[i*nodes + j] > 0 && neighbours_colors[j] > 0) {
//				printf("[rank %d] neighbour %d has color %d\n", pid, j, neighbours_colors[j]);
				pblSetAdd(psetcolors, (void *) 	neighbours_colors[j]);
			}
		}
		
		last_color = find_my_color(psetcolors, nodes);
//		printf("[rank %d] color of %d is now %d\n", pid, i, last_color);
		
		my_colors[i] = last_color;

		for ( b = 0; b < p; b++ ) {

			//printf("[rank %d] before - last_color = %d | neighbour_color = %d \n", pid, last_color, neighbour_color);

			if ( pid == b+1 )
				neighbour_color = last_color;

			MPI_Bcast( &neighbour_color, 1, MPI_INT, b, *comm_workers );
			
			neighbours_colors[b * nodes_per_block + i] = neighbour_color;

			//printf("[rank %d] received - neighbour %d has color %d\n", pid, b*nodes_per_block+i, neighbour_color);

		}
	}


	if (my_nodes > nodes_per_block){
		for (i = nodes_per_block; i < max_nodes_per_block; i++){
			pblSetClear( psetcolors );

			
			//printf("[rank %d] neighbours =", pid);
			 /*for ( j = 0; j < nodes; j++ )
				printf(" %d", neighbours_colors[j]);
			printf("\n");*/


			for ( j = 0; j < nodes; j++ ) {
				if (j == (pid-1)*nodes_per_block + i )
					continue;


				if (madjacency[i*nodes + j] > 0 && neighbours_colors[j] > 0) {
					//printf("[rank %d] neighbour %d has color %d\n", pid, j, neighbours_colors[j]);
					pblSetAdd(psetcolors, (void *) 	neighbours_colors[j]);
				}
			}
			
			my_colors[i] = find_my_color(psetcolors, nodes);
			neighbours_colors[(pid-1)*nodes_per_block + i] = my_colors[i];
		}
	}



	// STEP 2
	for ( i = 0; i < nodes_per_block; i++ ) {
		
		for ( j = 0; j < nodes; j++ ) {
			if (j == (pid-1)*nodes_per_block + i ) {
				continue;
			}
			//printf("[rank %d] %d %d > adjacency %d > j mod nodes_per_block == i %d > neighbours_colors[j] = %d > my_colors[i] %d\n", 
			//	pid, (pid-1)*nodes_per_block + i, j, madjacency[i*nodes + j], (j % nodes_per_block == i), neighbours_colors[j], my_colors[i]);
				
			if (madjacency[i * nodes + j] > 0 && (j % nodes_per_block == i) &&  neighbours_colors[j] == my_colors[i]) 
			{
				//printf("[rank %d] have conflict %d vs %d\n", pid, (pid-1)*nodes_per_block + i, j);
				if ( (pid-1) * nodes_per_block + i < j ) 
				{
					//printf("[rank %d] conflict++\n", pid);
					conflicts[n_conflicts] = (pid-1)*nodes_per_block + i;
					n_conflicts++;
					break;
				}
			}
		}
	}


	/*printf("[rank %d] my_colors = ", pid);
	for (i = 0; i < my_nodes; i++)
		printf(" %d", my_colors[i]);
	printf("\n");*/

	MPI_Send( my_colors, my_nodes, MPI_INT, 0, TAG, MPI_COMM_WORLD);
	MPI_Send( conflicts, n_conflicts, MPI_INT, 0, TAG, MPI_COMM_WORLD);
}

void plassman_v_processors(int pid) {
	int nodes;
	int * madjacency;
	int weight, n_wait = 0, n_send = 0, i, neigbour_weight, tmp_color, color;
	int * send_to;
	int * receive_from;
	int * colors;

	MPI_Request req;
	MPI_Status status;
	

	MPI_Recv( &nodes, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status );

	send_to = (int*) malloc( (nodes-1) * sizeof(int) );
	receive_from = (int*) malloc( (nodes-1) * sizeof(int) );
	colors = (int*) malloc( (nodes) * sizeof(int) );
	
	madjacency = (int*)malloc(nodes * sizeof(int));
	MPI_Recv( madjacency, nodes, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status );

	//defenir peso random
	srand( time( NULL ) * pid );
	weight = rand();

	
	// enviar peso para e receber peso dos vizinhos, definindo de quais vai esperar e de quais vai enviar
	for (i = 1; i <= nodes; i++) {

		if (i == pid || madjacency[i-1] == 0)
			continue;
		MPI_Isend( &weight, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &req );
		
		MPI_Recv( &neigbour_weight, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &status );

		if ( weight < neigbour_weight ) {
			receive_from[n_wait] = i;
			n_wait++;
		} 
		else {
			send_to[n_send] = i;
			n_send++;
		}

		MPI_Wait(&req, &status);
	}
	
	PblSet * psetcolors = (PblSet *) pblSetNewHashSet();

	// esperar pesos maiores
	for ( i = 0; i < n_wait; i++ ) {
		MPI_Recv( &tmp_color, 1, MPI_INT, receive_from[i], TAG, MPI_COMM_WORLD, &status );
		pblSetAdd(psetcolors, (void *) tmp_color);
	}
	
	// definir cor do no
	color = find_my_color(psetcolors, nodes);
	
	// enviar aos pesos menores
	for ( i = 0; i < n_send; i++) {
	    MPI_Send( &color, 1, MPI_INT, send_to[i], TAG, MPI_COMM_WORLD);
    }
	
	//enviar cor final para o master
	//printf("[rank %d] Send color %d to master\n", pid, color);
	MPI_Send( &color, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
}

void plassman_p_processors(int pid, int p, MPI_Comm * comm_workers) {
	int nodes;
	int * madjacency;
	int * all_weights, i, j;
	int my_nodes, nodes_per_block, max_nodes_per_block, n_locals;
	int * send_to, * n_wait, * n_send, * locals;
	int * receive_from, ended = 0;
	int * colors, my_colors[2], tmp = 0;

	MPI_Request req;
	MPI_Status status;
	PblSet * psetcolors = (PblSet *) pblSetNewHashSet();
	int good = 1;	

	MPI_Recv( &nodes, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status );

	nodes_per_block = get_nodes_per_process(nodes, p);
	max_nodes_per_block = nodes - (p - 1) * nodes_per_block;
	
	//printf("[rank %d] nodes: %d | nodes_per_block: %d | max_nodes_per_block: %d\n", pid, nodes, nodes_per_block, max_nodes_per_block);

	if (pid < p)
		my_nodes = nodes_per_block;
	else
		my_nodes = max_nodes_per_block;
		
	n_locals = 0;
	//printf("[rank %d] nodes: %d\n", pid, nodes);

	n_wait = (int*) malloc( my_nodes * sizeof(int) );
	n_send = (int*) malloc( my_nodes * sizeof(int) );
	send_to = (int*) malloc( my_nodes*(nodes-1) * sizeof(int) );
	receive_from = (int*) malloc( my_nodes*(nodes-1) * sizeof(int) );
	
	colors = (int*) malloc( (nodes) * sizeof(int) );
	all_weights = (int*) malloc( (nodes) * sizeof(int) );
	madjacency = (int*) malloc( nodes * my_nodes * sizeof(int) );
	locals = (int*) malloc( my_nodes * sizeof(int) );
	MPI_Recv( madjacency, nodes*my_nodes, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status );

	//defenir peso random
	srand( time( NULL ) / pid );
	for (i = 0; i < my_nodes; i++){
		all_weights[(pid-1)*nodes_per_block +i] = rand();
		n_wait[i] = 0;
		n_send[i] = 0;
		//print-f("[rank %d] weight[%d]: %d\n", pid, (pid-1)*nodes_per_block + i, all_weights[(pid-1)*nodes_per_block +i]);
	}
		

	//printf("[rank %d] will exchange weights\n", pid);
	//exchange all weights
	for ( i = 0; i < p-1; i++)
		MPI_Bcast( &all_weights[i*nodes_per_block], nodes_per_block, MPI_INT, i, *comm_workers );
	MPI_Bcast( &all_weights[(p-1)*nodes_per_block], max_nodes_per_block, MPI_INT, p-1, *comm_workers );
	
	//printf("[rank %d] weights exchanged\n", pid);

	fflush(stdout);

	for ( i = 0; i < my_nodes; i++ ) {
		for (j = 0; j < nodes; j++) {

			if ((pid-1)*nodes_per_block + i == j || madjacency[i*nodes + j] <= 0 || pid == get_pid_of_node(j, nodes, p))
				continue;

			if ( all_weights[(pid-1)*nodes_per_block + i] < all_weights[j] ) {
				receive_from[i*(nodes-1) + n_wait[i]] = j;
				n_wait[i]++;
			} 
			else {
				send_to[i*(nodes-1) + n_send[i]] = j;
				n_send[i]++;
			}
		}

		locals[i] = is_local( (pid-1)*nodes_per_block+i, i, madjacency, nodes, p);
		if (locals[i]) {
			n_locals++;
			continue;
		}

		//printf("[rank %d] vertice[%d]: wait %d | send %d\n", pid, (pid-1)*nodes_per_block+i, n_wait[i], n_send[i]);
		if(n_wait[i] == 0 ){
			colors[(pid-1)*nodes_per_block+i] = color_node(i, pid,nodes_per_block, madjacency, colors, psetcolors, nodes);
			ended++;
			
			my_colors[0] = (pid-1)*nodes_per_block+i;
			my_colors[1] = colors[my_colors[0]];
			for (j = 0; j < n_send[i]; j++) {
				//printf("[rank %d] send color %d do vertice %d para o %d do rank %d\n", pid, my_colors[1], (pid-1)*nodes_per_block+i,
				//		send_to[i*(nodes-1) +j], get_pid_of_node(send_to[i*(nodes-1) +j],nodes,p));

				MPI_Isend(my_colors, 2, MPI_INT, get_pid_of_node(send_to[i*(nodes-1) + j], nodes, p), TAG, MPI_COMM_WORLD, &req);
			}
		}
	}
	
	//print_adjacency_matrix(pid, madjacency, my_nodes, nodes);
	//printf("[rank %d] locals");
	//for (i = 0; i < my_nodes; i++)
		//printf(" %d", locals[i]);
	//printf("\n");

	//printf("[rank %d] ended: %d | locals = %d\n", pid, ended, n_locals);

	while(ended != my_nodes-n_locals) {
		//printf("[rank %d] mais uma voltinha...\n", pid);

		MPI_Recv(my_colors, 2, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, &status);

		colors[my_colors[0]] = my_colors[1];
		//printf("[rank %d] Reiceived color %d from vertice %d\n", pid, my_colors[1], my_colors[0]);
		//printf("[rank %d]", pid);
		//for( i = 0; i < nodes; i++)
			//printf(" %d", colors[i]);
		//printf("\n");


		for ( i = 0; i < my_nodes; i++ ) {
			if (n_wait[i] == 0 || locals[i])
				continue;
			//printf("[rank %d] %d not done...\n", pid, i);
			for (j = 0; j < n_wait[i]; j++)
				if (receive_from[i*(nodes-1)+j] == my_colors[0] ) {
					receive_from[i*(nodes-1)+j] = receive_from[i*(nodes-1)+n_wait[i]-1];
					n_wait[i]--;
					break;
				}
				
			if (n_wait[i] > 0)
				continue;
			//printf("[rank %d] %d done...\n", pid, i);
			colors[(pid-1)*nodes_per_block+i] = color_node(i, pid,nodes_per_block, madjacency, colors, psetcolors, nodes);
			ended++;
			
			my_colors[0] = (pid-1)*nodes_per_block+i;
			my_colors[1] = colors[my_colors[0]];
			for (j = 0; j < n_send[i]; j++) {
				//printf("[rank %d] send color do vertice %d para o %d do rank %d\n", pid, (pid-1)*nodes_per_block+i,
				//		send_to[i*(nodes-1) +j], get_pid_of_node(send_to[i*(nodes-1) +j],nodes,p));

				MPI_Isend(my_colors, 2, MPI_INT, get_pid_of_node(send_to[i*(nodes-1) + j], nodes, p), TAG, MPI_COMM_WORLD, &req);
			}
		}

	}

	for(i = 0; i < my_nodes; i++)
		if (locals[i])
			colors[(pid-1)*nodes_per_block + i] = color_node(i, pid,nodes_per_block, madjacency, colors, psetcolors, nodes);
	//enviar cor final para o master
//	for(i = 0; i < my_nodes; i++)
		//printf("[rank %d] color[%d] = %d\n", pid, i, colors[(pid-1)*nodes_per_block + i]);
	//printf("[rank %d] Send colors to master\n", pid);
	MPI_Send( &colors[(pid-1)*nodes_per_block], my_nodes, MPI_INT, 0, TAG, MPI_COMM_WORLD);
}


void greedy_parallel(int pid) {

	int nodes, node, i, tmp;
	int color, othernode;
	MPI_Status status;
	MPI_Request req;
	int * indexes, *colors, *madjacency;

	PblSet * psetcolors = (PblSet *) pblSetNewHashSet();
	

	MPI_Recv( &nodes, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status );

	//printf("[rank %d] NODES: %d\n", pid, nodes);

	srand( time(NULL) / pid );

	indexes = (int*) malloc(nodes*sizeof(int));
	colors =  (int*) malloc(nodes*sizeof(int));
	madjacency  = (int*) malloc(nodes*nodes*sizeof(int));

	MPI_Recv(madjacency, nodes*nodes, MPI_INT, 0, TAG, MPI_COMM_WORLD, &status);

	for (i = 0; i < nodes; i++)
		indexes[i] = i;

	for (i = nodes-1; i >= 0; i--) {
		tmp = rand_between(0, i);
		
		node = indexes[tmp];

		indexes[tmp] = indexes[i];

		pblSetClear(psetcolors);
		
		for (othernode = 0; othernode < nodes; othernode++) {
			if (othernode != node) {
				if (madjacency[POS(othernode, node, nodes)] == 1) {
					if (colors[othernode] != 0)
					/* the neighbor already has a color */
						pblSetAdd(psetcolors, (void *) colors[othernode]);
				}
				
			}
		}

		colors[node] = find_my_color(psetcolors, nodes);
	}

	//printf("[rank %d] colors = ", pid);
	//for (i = 0; i < nodes; i++)
	//	printf(" %d", colors[i]);
	//printf("\n");

	MPI_Send(colors, nodes, MPI_INT, 0, TAG, MPI_COMM_WORLD);
}




