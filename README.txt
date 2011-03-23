1 - Para compilar a nossa aplicação basta escrever make na pasta onde se encontram os ficheiros.
Serão gerados dois executáveis: um com uma versão sequencial, chamado greedy_colouring, e outro com versões paralelas, chamado parallel_colouring.

2 - a) Para correr a versão sequencial, basta na consola escrever:
		./greedy_colouring < [ficheiro_input]

2 - b) Para correr a versão paralela, basta escrever na consola:
		mpiexec -n [n_processors] ./parallel_colouring [nome_algoritmo] < [input_file]
		
		O parâmetro nome_algoritmo pode ter os seguintes valores:
			- block_partition
			- plassman_p
			- plassman_v
			- greedy_parallel