/*
 *  common.h
 *  P2
 *
 *  Created by Filipe Araujo on 4/28/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#define POS(x, y, side) ((x) + (y) * (side))



void error(char * msg);
void read_graph_to_adjacency_matrix(FILE * f, int ** pmatrix, int * pnodes);

