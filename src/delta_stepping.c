/*
 ============================================================================
 Name        : delta_stepping.c
 Author      : stam
 Version     :
 Copyright   : Your copyright notice
 Description : Calculate Pi in MPI
 ============================================================================
 */
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#define NEWARRAY(x, n) do { (x) = calloc((n), sizeof *(x)); } while (0)
typedef struct NODE{
	int nodeId;
	int *edgeList;
} NODE;

NODE* New_grapg(int numberOfNodes){
	NODE *newGraph= malloc(sizeof(NODE)*numberOfNodes);
	    assert(newGraph != NULL);
	return newGraph;
}

int 
main(int argc, char *argv[])
{
	int			my_rank;		/* rank of process */
	int			num_procs;		/* number of processes */
	int			source;			/* rank of sender */
	int			dest = 0;		/* rank of receiver */
	int			tag = 0;		/* tag for messages */
	char		message[100];	/* storage for message */
	MPI_Status	status ;		/* return status for receive */
	char		*nameOfGraphFile;
	char line[80];
	long elapsed_seconds;
	FILE *fr;            /* declare the file pointer */
	int sizeOfGraph;
	NODE *graph;

	/* start up MPI */


	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 
	if(argc < 2){
		printf("No Arguments \n Arguments:\n name of graph file \nintial delta value ");
		exit(1);
    }
	nameOfGraphFile=argv[1];
	double delta=atof(argv[2]);
	if(nameOfGraphFile!=NULL && nameOfGraphFile!=""){
		fr = fopen (nameOfGraphFile, "r");  /* open the file for reading */
		   /* elapsed.dta is the name of the file */
		   /* "rt" means open the file for reading text */
		/* read the graph */
		if(fr==NULL){
			 printf("error: the file does not exists \n");
			 exit(1);
		}
		while(fgets(line, 80, fr) != NULL)
		{
			 /* get a line, up to 80 chars from fr.  done if NULL */
			 sscanf (line, "%ld", &elapsed_seconds);
			 /* convert the string to a double int */
			 sizeOfGraph=atoi(line);
			 if(sizeOfGraph<1){
				 printf("error: The size of the graph is  < 1 \n");
				 exit(1);
			 }
			 graph=New_grapg(sizeOfGraph)
			 if(atoi(line)!=-1){

			 }
			 printf ("line is %s\n", line);
		}
	}


	/* shut down MPI */
	MPI_Finalize(); 

	return 0;
}
