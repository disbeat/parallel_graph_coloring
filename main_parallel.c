#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <pbl.h>
#include <time.h>
#include "common.h"
#include "parallel.h"


int main(int argc, char * argv[]) {
    
    if(argc <= 1) {
        printf("Wrong parameters!!\n You have to introduce the name of the algorithm you want to run!!\n");
        printf("Instructions:\n");
        printf("mpirun -np <n_processes> ./<exec_file> plassman_p < <input_file> \n");
        printf("mpirun -np <n_processes> ./<exec_file> block_partition < <input_file> \n");
        printf("mpirun -np <n_processes> ./<exec_file> plassman_v < <input_file> \n");
        printf("mpirun -np <n_processes> ./<exec_file> greedy_parallel < <input_file> \n");
        return 0;
    }
    else {
        
        if(strcmp("block_partition",argv[1]) == 0){
	
            if(block_partioning_algorithm(argc,argv) == -1) {
                printf("An error ocurred while using the block partion algorithm\n");
            }
    
        } 
        else if (strcmp("plassman_p",argv[1]) == 0) {
			if(plassman_algorithm_p_processors(argc,argv) == -1) {
				printf("An error ocurred while using the plassman algorithm\n");
			}
		}
		else if (strcmp("plassman_v",argv[1]) == 0) {
			if(plassman_algorithm_v_processors(argc,argv) == -1) {
				printf("An error ocurred while using the plassman algorithm v processors\n");
			}
		}else if (strcmp("greedy_parallel",argv[1]) == 0) {
			if(main_greedy_parallel(argc,argv) == -1) {
				printf("An error ocurred while using the greedy parallel algorithm\n");
			}
		}
        
    }
    return 0;
}



int block_partioning_algorithm(int argc, char * argv[]) {
    int * matrix, * colors, * conflicts;
	int nbrnodes, nbrcolors;
	int i, j, err, blocks;
	int  numtasks, rank, conflict_pos = 0, read_items;
	static int excl_ranks[1] = {0};
	clock_t ticks1, ticks2;
	MPI_Request * reqs;
	MPI_Status * status;

	MPI_Group orig_group, workers_group;
	MPI_Comm comm_workers;

	PblSet * psetcolors = (PblSet *) pblSetNewHashSet();
    ticks1 = clock();
	
	err = MPI_Init(&argc, &argv); /* Initialize MPI */
	if (err != MPI_SUCCESS) {
		printf("MPI_init failed!\n"); return 1;
	}

	err = MPI_Comm_size(MPI_COMM_WORLD, &numtasks);	/* Get nr of tasks */
	if (err != MPI_SUCCESS) {
		printf("MPI_Comm_size failed!\n"); return 1;
	}

	err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);    /* Get id of this process */
	if (err != MPI_SUCCESS) {
		printf("MPI_Comm_rank failed!\n"); return 1;
	}
	

	// create the workers group = all \ {master}
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

	MPI_Group_excl(orig_group, 1, excl_ranks, &workers_group);

	MPI_Comm_create(MPI_COMM_WORLD, workers_group, &comm_workers);


	reqs = (MPI_Request*)malloc((numtasks-1)*sizeof(MPI_Request));
	status = (MPI_Status*)malloc((numtasks-1)*sizeof(MPI_Status));

	if (rank == 0) {
		read_graph_to_adjacency_matrix(stdin, &matrix, &nbrnodes);
		
		colors = (int*)malloc(nbrnodes*sizeof(int));
		conflicts = (int*)malloc(nbrnodes*sizeof(int));

		//print_adjacency_matrix(0, matrix, nbrnodes, nbrnodes);
	    //printf("nbrnodes = %d\n", nbrnodes);
        
		blocks = floor(nbrnodes *1.0/ (numtasks-1));
		for (i = 1; i < numtasks-1; i++){
			MPI_Send( &nbrnodes, 1, MPI_INT, i, TAG, MPI_COMM_WORLD );
			MPI_Send( &matrix[(i-1) * nbrnodes * blocks], nbrnodes*blocks, MPI_INT, i, TAG, MPI_COMM_WORLD );
		}

		MPI_Send( &nbrnodes, 1, MPI_INT, numtasks-1, TAG, MPI_COMM_WORLD );
		MPI_Send( &matrix[(numtasks-2) * nbrnodes * blocks], nbrnodes*(nbrnodes-(numtasks-2)*blocks), MPI_INT, numtasks-1, TAG, MPI_COMM_WORLD );

		for (i = 1; i < numtasks-1; i++)
			MPI_Irecv( &colors[(i-1)*blocks], blocks, MPI_INT, i, TAG, MPI_COMM_WORLD, &reqs[i-1] );
		MPI_Irecv( &colors[(numtasks-2)*blocks], nbrnodes - (numtasks-2)*blocks, MPI_INT, i, TAG, MPI_COMM_WORLD, &reqs[numtasks-2] );
		MPI_Waitall(numtasks-1, reqs, status);
		
		for (i = 1; i < numtasks; i++) {
    		MPI_Recv( &conflicts[conflict_pos], nbrnodes, MPI_INT, i, TAG, MPI_COMM_WORLD, &status[i-1] );
    		MPI_Get_count(&status[i-1], MPI_INT, &read_items);

			conflict_pos += read_items;
		}

		// resolve conflicts

		for ( i = 0; i < conflict_pos; i++ ) {
			pblSetClear(psetcolors);
			for (j = 0; j < nbrnodes; j++) {
				if (j == conflicts[i]) 
					continue;

				if (matrix[POS(conflicts[i], j, nbrnodes)] == 1) {
					pblSetAdd(psetcolors, (void *) colors[j]);
				}
			}

			colors[conflicts[i]] = find_my_color(psetcolors, nbrnodes);
		}

		/*for ( i = 0; i < conflict_pos; i++ )
			printf("conflict index %d\n", conflicts[i]);*/

		nbrcolors = 0;
		// FIND MAX COLOR
		for (i = 1; i < nbrnodes; i++)
			if (colors[i] > nbrcolors)   
				nbrcolors = colors[i];

		printf("%d\n%d\n", nbrnodes, nbrcolors);
		for (i = 0; i < nbrnodes; i++) {
			printf("%d\n", colors[i]);
		}
        ticks2	= clock();
        printf("Time Elapsed = %lf\n", 1.0 * (ticks2 - ticks1) / CLOCKS_PER_SEC);

	} else {	
		block_partition(rank, numtasks-1, &comm_workers);
	}
	

	err = MPI_Finalize();	         /* Terminate MPI */
	if (err != MPI_SUCCESS) {
		printf("Error in MPI_Finalize!\n");
		return -1;
	}
	
    return 0;
    
}



int plassman_algorithm_p_processors(int argc, char * argv[]) {
	int * matrix, * colors, * conflicts;
	int nbrnodes, nbrcolors;
	int i, j, err, blocks;
	int  numtasks, rank, conflict_pos = 0, read_items;
	static int excl_ranks[1] = {0};
	clock_t ticks1, ticks2;
	MPI_Request * reqs;
	MPI_Status * status;

	MPI_Group orig_group, workers_group;
	MPI_Comm comm_workers;

	PblSet * psetcolors = (PblSet *) pblSetNewHashSet();
    ticks1 = clock();
	
	err = MPI_Init(&argc, &argv); /* Initialize MPI */
	if (err != MPI_SUCCESS) {
		printf("MPI_init failed!\n"); return 1;
	}

	err = MPI_Comm_size(MPI_COMM_WORLD, &numtasks);	/* Get nr of tasks */
	if (err != MPI_SUCCESS) {
		printf("MPI_Comm_size failed!\n"); return 1;
	}

	err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);    /* Get id of this process */
	if (err != MPI_SUCCESS) {
		printf("MPI_Comm_rank failed!\n"); return 1;
	}
	

	// create the workers group = all \ {master}
	MPI_Comm_group(MPI_COMM_WORLD, &orig_group);

	MPI_Group_excl(orig_group, 1, excl_ranks, &workers_group);

	MPI_Comm_create(MPI_COMM_WORLD, workers_group, &comm_workers);

	if (rank == 0) {
		read_graph_to_adjacency_matrix(stdin, &matrix, &nbrnodes);

		colors = (int*)malloc(nbrnodes*sizeof(int));
		conflicts = (int*)malloc(nbrnodes*sizeof(int));

		reqs = (MPI_Request*)malloc((numtasks-1)*sizeof(MPI_Request));
		status = (MPI_Status*)malloc((numtasks-1)*sizeof(MPI_Status));


		//print_adjacency_matrix(0, matrix, nbrnodes, nbrnodes);
	    //printf("nbrnodes = %d\n", nbrnodes);
        
		blocks = floor(nbrnodes *1.0/ (numtasks-1));
		for (i = 1; i < numtasks-1; i++){
			MPI_Send( &nbrnodes, 1, MPI_INT, i, TAG, MPI_COMM_WORLD );
			MPI_Send( &matrix[(i-1) * nbrnodes * blocks], nbrnodes*blocks, MPI_INT, i, TAG, MPI_COMM_WORLD );
		}

		MPI_Send( &nbrnodes, 1, MPI_INT, numtasks-1, TAG, MPI_COMM_WORLD );
		MPI_Send( &matrix[(numtasks-2) * nbrnodes * blocks], nbrnodes*(nbrnodes-(numtasks-2)*blocks), MPI_INT, numtasks-1, TAG, MPI_COMM_WORLD );

		for (i = 1; i < numtasks-1; i++)
			MPI_Irecv( &colors[(i-1)*blocks], blocks, MPI_INT, i, TAG, MPI_COMM_WORLD, &reqs[i-1] );
		MPI_Irecv( &colors[(numtasks-2)*blocks], nbrnodes - (numtasks-2)*blocks, MPI_INT, i, TAG, MPI_COMM_WORLD, &reqs[numtasks-2] );
		MPI_Waitall(numtasks-1, reqs, status);
		
	
		nbrcolors = 0;
		// FIND MAX COLOR
		for (i = 1; i < nbrnodes; i++)
			if (colors[i] > nbrcolors)   
				nbrcolors = colors[i];

		printf("%d\n%d\n", nbrnodes, nbrcolors);
		for (i = 0; i < nbrnodes; i++) {
			printf("%d\n", colors[i]);
		}
        ticks2	= clock();
        printf("Time Elapsed = %lf\n", 1.0 * (ticks2 - ticks1) / CLOCKS_PER_SEC);

	} else {	
		plassman_p_processors(rank, numtasks-1, &comm_workers);
	}
	

	err = MPI_Finalize();	         /* Terminate MPI */
	if (err != MPI_SUCCESS) {
		printf("Error in MPI_Finalize!\n");
		return -1;
	}
	
    return 0;
}



int plassman_algorithm_v_processors(int argc, char * argv[]) {
	int * matrix, * colors;
	int nbrnodes, nbrcolors;
	int i, err;
	int  numtasks, rank;
	clock_t ticks1, ticks2;

	MPI_Request * reqs;
	MPI_Status * status;

	ticks1 = clock();

	err = MPI_Init(&argc, &argv); /* Initialize MPI */
	if (err != MPI_SUCCESS) {
		printf("MPI_init failed!\n");
		return 1;
	}

	err = MPI_Comm_size(MPI_COMM_WORLD, &numtasks);	/* Get nr of tasks */
	if (err != MPI_SUCCESS) {
		printf("MPI_Comm_size failed!\n");
		return 1;
	}

	err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);    /* Get id of this process */
	if (err != MPI_SUCCESS) {
		printf("MPI_Comm_rank failed!\n");
		return 1;
	}

	reqs = (MPI_Request*)malloc((numtasks-1)*sizeof(MPI_Request));
	status = (MPI_Status*)malloc((numtasks-1)*sizeof(MPI_Status));

	if (rank == 0) {
		read_graph_to_adjacency_matrix(stdin, &matrix, &nbrnodes);
		colors = (int*)malloc(nbrnodes*sizeof(int));
	
		if (numtasks != nbrnodes+1) {
			printf("ERROR:  For that input you must execute %d processes but you did %d", nbrnodes+1, numtasks);
			err = MPI_Finalize();	         /* Terminate MPI */
			if (err != MPI_SUCCESS) {
				printf("Error in MPI_Finalize!\n");
				return 1;
			}
			return -1;
		}
	    //printf("nbrnodes = %d\n", nbrnodes);
        
		for (i = 1; i < numtasks; i++) {
			MPI_Send( &nbrnodes, 1, MPI_INT, i, TAG, MPI_COMM_WORLD );
			MPI_Send( &matrix[(i-1) * nbrnodes], nbrnodes, MPI_INT, i, TAG, MPI_COMM_WORLD );
		}

		for (i = 1; i < numtasks; i++)
			MPI_Irecv( &colors[i-1], 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &reqs[i-1] );
		MPI_Waitall(numtasks-1, reqs, status);
		
		nbrcolors = 0;
		// FIND MAX COLOR
		for (i = 0; i < numtasks-1; i++)
			if (colors[i] > nbrcolors)   
				nbrcolors = colors[i];

		printf("%d\n%d\n", nbrnodes, nbrcolors);
		for (i = 0; i < nbrnodes; i++)
			printf("%d\n", colors[i]);
		ticks2	= clock();
        printf("Time Elapsed = %lf\n", 1.0 * (ticks2 - ticks1) / CLOCKS_PER_SEC);

	}else{
		plassman_v_processors(rank);
	}

	err = MPI_Finalize();	         /* Terminate MPI */
	if (err != MPI_SUCCESS) {
		printf("Error in MPI_Finalize!\n");
		return 1;
	}

    return 0;
}
int main_greedy_parallel(int argc, char * argv[]) {
	int * matrix, * colors, * conflicts;
	int nbrnodes, nbrcolors,best_solution_pos,best_solution_max_colors;
	int i, j, err, blocks;
	int  numtasks, rank, conflict_pos = 0, read_items;
	static int excl_ranks[1] = {0};
	clock_t ticks1, ticks2;
	MPI_Request * reqs;
	MPI_Status * status;

	MPI_Group orig_group, workers_group;
	MPI_Comm comm_workers;

	PblSet * psetcolors = (PblSet *) pblSetNewHashSet();
    ticks1 = clock();
	
	

	colors = (int*)malloc(nbrnodes*sizeof(int));

	err = MPI_Init(&argc, &argv); /* Initialize MPI */
	if (err != MPI_SUCCESS) {
		printf("MPI_init failed!\n"); return 1;
	}

	err = MPI_Comm_size(MPI_COMM_WORLD, &numtasks);	/* Get nr of tasks */
	if (err != MPI_SUCCESS) {
		printf("MPI_Comm_size failed!\n"); return 1;
	}

	err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);    /* Get id of this process */
	if (err != MPI_SUCCESS) {
		printf("MPI_Comm_rank failed!\n"); return 1;
	}
	


	reqs = (MPI_Request*)malloc((numtasks-1)*sizeof(MPI_Request));
	status = (MPI_Status*)malloc((numtasks-1)*sizeof(MPI_Status));

	if (rank == 0) {
		read_graph_to_adjacency_matrix(stdin, &matrix, &nbrnodes);

		//print_adjacency_matrix(0, matrix, nbrnodes, nbrnodes);
	    //printf("nbrnodes = %d\n", nbrnodes);
	    
	   
        for(i = 1; i < numtasks; i++) {
            MPI_Send( &nbrnodes, 1, MPI_INT, i, TAG, MPI_COMM_WORLD );
    		MPI_Send( matrix, nbrnodes*nbrnodes, MPI_INT, i, TAG, MPI_COMM_WORLD );   
        }
        
        
        for (i = 1; i < numtasks; i++){
			MPI_Irecv( &matrix[(i-1)*nbrnodes], nbrnodes, MPI_INT, i, TAG, MPI_COMM_WORLD, &reqs[i-1] );
		}
		MPI_Waitall(numtasks-1, reqs, status);
        
		
        best_solution_pos = 0;
        best_solution_max_colors = 0;
		// FIND MAX COLOR and best solucion
        for(i = 0; i < numtasks-1; i++) {
			nbrcolors = 0;
		    for (j = 0; j < nbrnodes; j++) {
			    if (matrix[i*nbrnodes + j] > nbrcolors) {   
				    nbrcolors = matrix[i*nbrnodes + j];
			    }
			}
			
			if(best_solution_max_colors > nbrcolors || best_solution_max_colors == 0) {
                best_solution_pos = i * nbrnodes;
                best_solution_max_colors = nbrcolors;
			}
			
				    
		}

		printf("%d\n%d\n", nbrnodes, best_solution_max_colors);
        for (i = best_solution_pos; i < best_solution_pos+nbrnodes; i++) {
			printf("%d\n", matrix[i]);
		}
        ticks2	= clock();
        printf("Time Elapsed = %lf\n", 1.0 * (ticks2 - ticks1) / CLOCKS_PER_SEC);

	} else {	
		greedy_parallel(rank);
	}
	

	err = MPI_Finalize();	         /* Terminate MPI */
	if (err != MPI_SUCCESS) {
		printf("Error in MPI_Finalize!\n");
		return -1;
	}
	
    return 0;
    
}
