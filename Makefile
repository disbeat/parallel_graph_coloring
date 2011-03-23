CC = mpicc
LINKER_FLAGS=-Lpbl_1_04
CFLAGS=-I pbl_1_04
OBJS = main.o common.o greedy.o
OBJS_PAR = parallel.o main_parallel.o common.o
PROG = greedy_colouring
PROG_PAR = parallel_colouring


all: ${PROG} ${PROG_PAR}

parallel: ${PROG_PAR}

clean:	
	rm ${OBJS} parallel.o main_parallel.o

${PROG}:	${OBJS}
			${CC} -o $@ ${OBJS} $(LINKER_FLAGS) -lpbl
			
${PROG_PAR}:	${OBJS_PAR}
				${CC} -o $@ ${OBJS_PAR} $(LINKER_FLAGS) -lpbl
			
			

