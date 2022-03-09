#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>

int main(int argc, char * argv[]) {
	//rank-ul procesului, numarul total de procese
	int rank, nProcesses;
	//lista workerilor fiecarui cluster
	int *clusters[3];
	//numarul de workeri pentru fiecare cluster
	int nWorkers[3];
	//parintele procesului(daca e cazul)
	int parent;
	int N, *V;
	int comErr;

	MPI_Init(&argc, &argv);
	MPI_Status status;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);

	comErr = atoi(argv[2]);

	//setare workers pentru fiecare coordonator
	if (rank <= 2) {
		//se deschide fisierul de input corespunzator coordonatorului
		std::string fileName = "cluster";
		fileName += (char)('0' + rank);
		fileName += ".txt";

		std::ifstream fin(fileName);

		//citim numarul de workers
		fin >> nWorkers[rank];

		clusters[rank] = (int*) malloc(nWorkers[rank] * sizeof(int));

		//se citesc workerii clusterului
		for (int i = 0; i < nWorkers[rank]; i++) {
			fin >> clusters[rank][i];
		}

		fin.close();
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//setare topologie pentru fiecare coordonator
	if (rank == 0) {
		if (!comErr) {
			//se trimite cluster-ul 0 catre coordonatorul 1
			MPI_Send(&nWorkers[0], 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("M(0,1)\n");

			MPI_Send(clusters[0], nWorkers[0], MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("M(0,1)\n");
		}

		//se trimite cluster-ul 0 catre coordonatorul 2
		MPI_Send(&nWorkers[0], 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
		printf("M(0,2)\n");

		MPI_Send(clusters[0], nWorkers[0], MPI_INT, 2, 0, MPI_COMM_WORLD);
		printf("M(0,2)\n");

		if (!comErr) {
			//se primeste cluster-ul 1 de la coordonatorul 1
			MPI_Recv(&nWorkers[1], 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
			clusters[1] = (int *) malloc(nWorkers[1] * sizeof(int));
			MPI_Recv(clusters[1], nWorkers[1], MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
		} else {
			//se primeste cluster-ul 1 de la coordonatorul 2
			MPI_Recv(&nWorkers[1], 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
			clusters[1] = (int *) malloc(nWorkers[1] * sizeof(int));
			MPI_Recv(clusters[1], nWorkers[1], MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
		}
		
		//se primeste cluster-ul 2 de la coordonatorul 2
		MPI_Recv(&nWorkers[2], 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
		clusters[2] = (int *) malloc(nWorkers[2] * sizeof(int));
		MPI_Recv(clusters[2], nWorkers[2], MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

		//se afiseaza topologia
		printf("0 -> 0:");

		int i;
		for(i = 0; i < nWorkers[0] - 1; i++) {
			printf("%d,", clusters[0][i]);
		}

		printf("%d 1:", clusters[0][i]);

		for(i = 0; i < nWorkers[1] - 1; i++) {
			printf("%d,", clusters[1][i]);
		}

		printf("%d 2:", clusters[1][i]);

		for (i = 0; i < nWorkers[2] - 1; i++) {
			printf("%d,", clusters[2][i]);
		}

		printf("%d\n", clusters[2][i]);
	}

	if (rank == 1) {
		if (!comErr) {
			//se primeste cluster-ul 0 de la coordonatorul 0
			MPI_Recv(&nWorkers[0], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			clusters[0] = (int *) malloc(nWorkers[0] * sizeof(int));
			MPI_Recv(clusters[0], nWorkers[0], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		} else {
			//se primeste cluster-ul 0 de la coordonatorul 2
			MPI_Recv(&nWorkers[0], 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
			clusters[0] = (int *) malloc(nWorkers[0] * sizeof(int));
			MPI_Recv(clusters[0], nWorkers[0], MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
		}

		if (!comErr) {
			//se trimite cluster-ul 1 catre coordonatorul 0
			MPI_Send(&nWorkers[1], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			printf("M(1,0)\n");

			MPI_Send(clusters[1], nWorkers[1], MPI_INT, 0, 0, MPI_COMM_WORLD);
			printf("M(1,0)\n");
		}

		//se trimite cluster-ul 1 catre coordonatorul 2
		MPI_Send(&nWorkers[1], 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
		printf("M(1,2)\n");

		MPI_Send(clusters[1], nWorkers[1], MPI_INT, 2, 0, MPI_COMM_WORLD);
		printf("M(1,2)\n");

		//se primeste cluster-ul 2 de la coordonatorul 2
		MPI_Recv(&nWorkers[2], 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
		clusters[2] = (int *) malloc(nWorkers[2] * sizeof(int));
		MPI_Recv(clusters[2], nWorkers[2], MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

		//se afiseaza topologia
		printf("1 -> 0:");

		int i;
		for(i = 0; i < nWorkers[0] - 1; i++) {
			printf("%d,", clusters[0][i]);
		}

		printf("%d 1:", clusters[0][i]);

		for(i = 0; i < nWorkers[1] - 1; i++) {
			printf("%d,", clusters[1][i]);
		}

		printf("%d 2:", clusters[1][i]);

		for (i = 0; i < nWorkers[2] - 1; i++) {
			printf("%d,", clusters[2][i]);
		}

		printf("%d\n", clusters[2][i]);
	}

	if (rank == 2) {
		//se primeste cluster-ul 0 de la coordonatorul 0
		MPI_Recv(&nWorkers[0], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		clusters[0] = (int *) malloc(nWorkers[0] * sizeof(int));
		MPI_Recv(clusters[0], nWorkers[0], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

		if (comErr) {
			//se trimite cluster-ul 0 catre coordonatorul 1
			MPI_Send(&nWorkers[0], 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("M(2,1)\n");

			MPI_Send(clusters[0], nWorkers[0], MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("M(2,1)\n");
		}

		//se primeste cluster-ul 1 de la coordonatorul 1
		MPI_Recv(&nWorkers[1], 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
		clusters[1] = (int *) malloc(nWorkers[1] * sizeof(int));
		MPI_Recv(clusters[1], nWorkers[1], MPI_INT, 1, 0, MPI_COMM_WORLD, &status);

		if (comErr) {
			//se trimite cluster-ul 1 catre coordonatorul 0
			MPI_Send(&nWorkers[1], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			printf("M(2,0)\n");

			MPI_Send(clusters[1], nWorkers[1], MPI_INT, 0, 0, MPI_COMM_WORLD);
			printf("M(2,0)\n");
		}

		//se trimite cluster-ul 2 catre coordonatorul 0
		MPI_Send(&nWorkers[2], 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		printf("M(%d,%d)\n", 2, 0);

		MPI_Send(clusters[2], nWorkers[2], MPI_INT, 0, 0, MPI_COMM_WORLD);
		printf("M(%d,%d)\n", 2, 0);

		//se trimite cluster-ul 2 catre coordonatorul 1
		MPI_Send(&nWorkers[2], 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		printf("M(%d,%d)\n", 2, 1);

		MPI_Send(clusters[2], nWorkers[2], MPI_INT, 1, 0, MPI_COMM_WORLD);
		printf("M(%d,%d)\n", 2, 1);

		//se afiseaza topologia
		printf("2 -> 0:");

		int i;
		for(i = 0; i < nWorkers[0] - 1; i++) {
			printf("%d,", clusters[0][i]);
		}

		printf("%d 1:", clusters[0][i]);

		for(i = 0; i < nWorkers[1] - 1; i++) {
			printf("%d,", clusters[1][i]);
		}

		printf("%d 2:", clusters[1][i]);

		for (i = 0; i < nWorkers[2] - 1; i++) {
			printf("%d,", clusters[2][i]);
		}

		printf("%d\n", clusters[2][i]);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//setare topologie pentru fiecare worker
	//coordonatorii trimit tuturor workerilor topologia
	if (rank <= 2) {
		for (int i = 0; i < nWorkers[rank]; i++) {
			//se trimit numarul de workeri pentru fiecare cluster
			MPI_Send(nWorkers, 3, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
			printf("M(%d,%d)\n", rank, clusters[rank][i]);

			//se trimit clusterele in ordine
			MPI_Send(clusters[0], nWorkers[0], MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
			printf("M(%d,%d)\n", rank, clusters[rank][i]);

			MPI_Send(clusters[1], nWorkers[1], MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
			printf("M(%d,%d)\n", rank, clusters[rank][i]);

			MPI_Send(clusters[2], nWorkers[2], MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
			printf("M(%d,%d)\n", rank, clusters[rank][i]);
		}
	}

	//workerii asteapta topologia de la coordonatori
	if (rank > 2) {
		//se astepata numarul de workeri pentru fiecare cluster
		MPI_Recv(nWorkers, 3, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		
		//se seteaza parintele
		parent = status.MPI_SOURCE;

		//se aloca memorie pentru fiecare cluster
		//se asteapta clusterele in ordine
		clusters[0] = (int *) malloc(nWorkers[0] * sizeof(int));
		MPI_Recv(clusters[0], nWorkers[0], MPI_INT, parent, 0, MPI_COMM_WORLD, &status);

		clusters[1] = (int *) malloc(nWorkers[1] * sizeof(int));
		MPI_Recv(clusters[1], nWorkers[1], MPI_INT, parent, 0, MPI_COMM_WORLD, &status);

		clusters[2] = (int *) malloc(nWorkers[2] * sizeof(int));
		MPI_Recv(clusters[2], nWorkers[2], MPI_INT, parent, 0, MPI_COMM_WORLD, &status);

		//se afiseaza topologia
		printf("%d -> 0:", rank);

		int i;
		for(i = 0; i < nWorkers[0] - 1; i++) {
			printf("%d,", clusters[0][i]);
		}

		printf("%d 1:", clusters[0][i]);

		for(i = 0; i < nWorkers[1] - 1; i++) {
			printf("%d,", clusters[1][i]);
		}

		printf("%d 2:", clusters[1][i]);

		for (i = 0; i < nWorkers[2] - 1; i++) {
			printf("%d,", clusters[2][i]);
		}

		printf("%d\n", clusters[2][i]);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//coordonatorul 0 citeste datele si le trimite mai departe
	if (rank == 0) {
		//citire dimensiune vector
		N = atoi(argv[1]);

		V = (int *) malloc(N * sizeof(int));

		//populare vector
		for (int i = 0; i < N; i++) {
			V[i] = i;
		}

		//calculare dimensiune calcul per worker
		int varsPerWorker = (int) ceil((float)N / (nProcesses - 3));

		//calculare inceput si dimensiune pentru fiecare cluster
		int start0 = 0;
		int total0 = varsPerWorker * nWorkers[0];

		int start1 = start0 + total0;
		int total1 = varsPerWorker * nWorkers[1];

		int start2 = start1 + total1;
		int total2 = N - start2;

		//trimitere date catre ceilalti coordonatori
		if (!comErr) {
			MPI_Send(&total1, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("M(0,1)\n");

			MPI_Send(V + start1, total1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("M(0,1)\n");
		} else {
			MPI_Send(&total1, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
			printf("M(0,2)\n");

			MPI_Send(V + start1, total1, MPI_INT, 2, 0, MPI_COMM_WORLD);
			printf("M(0,2)\n");
		}

		MPI_Send(&total2, 1, MPI_INT, 2, 0, MPI_COMM_WORLD);
		printf("M(0,2)\n");

		MPI_Send(V + start2, total2, MPI_INT, 2, 0, MPI_COMM_WORLD);
		printf("M(0,2)\n");

		//trimitere date catre workerii clusterului 0
		for (int i = 0; i < nWorkers[rank]; i++) {
			//calculare start interval pentru fiecare worker
			int startInterval = i * (double) total0 / nWorkers[rank];
			
			//trimitere lungime date
			MPI_Send(&varsPerWorker, 1, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
			printf("M(%d,%d)\n", rank, clusters[rank][i]);

			//trimitere date
			MPI_Send(V + startInterval, varsPerWorker, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
			printf("M(%d,%d)\n", rank, clusters[rank][i]);
		}

		for (int i = 0; i < nWorkers[rank]; i++) {
			//calculare start interval de primire pentru fiecare worker
			int startInterval = i * (double) total0 / nWorkers[rank];

			//asteptare date
			MPI_Recv(V + startInterval, varsPerWorker, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD, &status);
		}

		//asteptare date de la ceilalti coordonatori
		if (!comErr) {
			MPI_Recv(V + start1, total1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
		} else {
			MPI_Recv(V + start1, total1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
		}

		MPI_Recv(V + start2, total2, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
	}

	//ceilalti coordonatori se ocupa de share-uirea datelor la workerii lor
	if (!comErr) {
		if (rank == 1 || rank == 2) {
			//asteptare dimensiune vector
			MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			V = (int *) malloc(N * sizeof(int));

			//asteptare vector
			MPI_Recv(V, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			//trimitere date catre workerii clusterului
			for (int i = 0; i < nWorkers[rank]; i++) {
				//calculare start interval si dimensiune vector
				int startInterval = i * (double) N / nWorkers[rank];
				int endInterval = fmin((i + 1) * (double)N / nWorkers[rank], N);
				int diffInterval = endInterval - startInterval;

				//trimitere lungime date
				MPI_Send(&diffInterval, 1, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
				printf("M(%d,%d)\n", rank, clusters[rank][i]);

				//trimitere vector
				MPI_Send(V + startInterval, diffInterval, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
				printf("M(%d,%d)\n", rank, clusters[rank][i]);
			}

			//asteptare date de la workeri
			for (int i = 0; i < nWorkers[rank]; i++) {
				int startInterval = i * (double) N / nWorkers[rank];
				int endInterval = fmin((i + 1) * (double)N / nWorkers[rank], N);
				int diffInterval = endInterval - startInterval;

				MPI_Recv(V + startInterval, diffInterval, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD, &status);
			}

			//trimitere date catre cluster-ul 0
			MPI_Send(V, N, MPI_INT, 0, 0, MPI_COMM_WORLD);
			printf("M(%d,0)\n", rank);

			free(V);
		}
	} else {
		if (rank == 1) {
			//asteptare dimensiune vector
			MPI_Recv(&N, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

			V = (int *) malloc(N * sizeof(int));

			//asteptare vector
			MPI_Recv(V, N, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

			//trimitere date catre workerii clusterului
			for (int i = 0; i < nWorkers[rank]; i++) {
				//calculare start interval si dimensiune vector
				int startInterval = i * (double) N / nWorkers[rank];
				int endInterval = fmin((i + 1) * (double)N / nWorkers[rank], N);
				int diffInterval = endInterval - startInterval;

				//trimitere lungime date
				MPI_Send(&diffInterval, 1, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
				printf("M(%d,%d)\n", rank, clusters[rank][i]);

				//trimitere vector
				MPI_Send(V + startInterval, diffInterval, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
				printf("M(%d,%d)\n", rank, clusters[rank][i]);
			}

			//asteptare date de la workeri
			for (int i = 0; i < nWorkers[rank]; i++) {
				int startInterval = i * (double) N / nWorkers[rank];
				int endInterval = fmin((i + 1) * (double)N / nWorkers[rank], N);
				int diffInterval = endInterval - startInterval;

				MPI_Recv(V + startInterval, diffInterval, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD, &status);
			}

			//trimitere date catre cluster-ul 2
			MPI_Send(V, N, MPI_INT, 2, 0, MPI_COMM_WORLD);
			printf("M(1,2)\n");

			free(V);
		}

		if (rank == 2) {
			//asteptare vector coordonator 1
			int N1;
			MPI_Recv(&N1, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			int *V1 = (int *) malloc(N1 * sizeof(int));

			MPI_Recv(V1, N1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			MPI_Send(&N1, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("M(2,1)\n");

			MPI_Send(V1, N1, MPI_INT, 1, 0, MPI_COMM_WORLD);
			printf("M(2,1)\n");

			//asteptare dimensiune vector
			MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			V = (int *) malloc(N * sizeof(int));

			//asteptare vector
			MPI_Recv(V, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

			//trimitere date catre workerii clusterului
			for (int i = 0; i < nWorkers[rank]; i++) {
				//calculare start interval si dimensiune vector
				int startInterval = i * (double) N / nWorkers[rank];
				int endInterval = fmin((i + 1) * (double)N / nWorkers[rank], N);
				int diffInterval = endInterval - startInterval;

				//trimitere lungime date
				MPI_Send(&diffInterval, 1, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
				printf("M(%d,%d)\n", rank, clusters[rank][i]);

				//trimitere vector
				MPI_Send(V + startInterval, diffInterval, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD);
				printf("M(%d,%d)\n", rank, clusters[rank][i]);
			}

			//asteptare date de la workeri
			for (int i = 0; i < nWorkers[rank]; i++) {
				int startInterval = i * (double) N / nWorkers[rank];
				int endInterval = fmin((i + 1) * (double)N / nWorkers[rank], N);
				int diffInterval = endInterval - startInterval;

				MPI_Recv(V + startInterval, diffInterval, MPI_INT, clusters[rank][i], 0, MPI_COMM_WORLD, &status);
			}

			//asteptare date de la coordonatorul 1
			MPI_Recv(V1, N1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);

			//trimitere date la coordonatorul 0
			MPI_Send(V1, N1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			printf("M(2,0)\n");

			//trimitere date catre cluster-ul 0
			MPI_Send(V, N, MPI_INT, 0, 0, MPI_COMM_WORLD);
			printf("M(2,0)\n");

			free(V);
			free(V1);
		}
	}

	//procesare date de catre workeri
	if (rank > 2) {
		//asteptare lungime vector
		MPI_Recv(&N, 1, MPI_INT, parent, 0, MPI_COMM_WORLD, &status);

		V = (int *) malloc(N * sizeof(int));

		//asteptare vector
		MPI_Recv(V, N, MPI_INT, parent, 0, MPI_COMM_WORLD, &status);

		//procesare date
		for (int i = 0; i < N; i++) {
			V[i] *= 2;
		}

		//trimitere date catre coordonator
		MPI_Send(V, N, MPI_INT, parent, 0, MPI_COMM_WORLD);
		printf("M(%d,%d)\n", rank, parent);

		free(V);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//afisare rezultat
	if (rank == 0) {
		printf("Rezultat:");

		for (int i = 0; i < N; i++) {
			printf(" %d", V[i]);
		}

		printf("\n");
	}

	//eliberare memorie
	{
		free(clusters[0]);
		free(clusters[1]);
		free(clusters[2]);

		if (rank == 0) {
			free(V);
		}
	}

	MPI_Finalize();
	return 0;
}
