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
#include <assert.h>
#define NEWARRAY(x, n) do { (x) = calloc((n), sizeof *(x)); } while (0)
#define RESIZEARRAY(x, n) do { (x) = realloc((n), sizeof *(x)); } while (0)


//the graph will be represented as an array of nodes
typedef struct NODE{
	int nodeId;//the id of the node
	int *edgeList;//the array containing the adj node ids
	double *weigthsList;//the weight to each edge
	int maxEdges;//to dynamically allovate new size for the edge list
	int numberOfEdges;
} NODE;

typedef struct BUCKETELEMENT{
	double* tentDist;
	int* nodeIds;
	int maxNodes;
	int numberOfNodes;
} BUCKETELEMENT;

NODE* New_grapg(int numberOfNodes){
	NODE *newGraph= calloc(numberOfNodes, sizeof(NODE));
	    assert(newGraph != NULL);

	newGraph->numberOfEdges=0;
	return newGraph;
}

BUCKETELEMENT* New_Bucket(int numberOfBuckets){
	BUCKETELEMENT *newBucket = calloc(numberOfBuckets, sizeof(BUCKETELEMENT));
		assert(newBucket!=NULL);

	newBucket->numberOfNodes=0;
	return newBucket;
}

/* rangefrom inclusive end not inclusive */
int get_randomNumber(int rangefrom,int rangeTo)
{
	return (rand() % rangeTo + rangefrom);
}

void pushEdge(NODE *anode,int edgeId,double edgeWeight){
	if(anode->edgeList==NULL){
		NEWARRAY(anode->edgeList,10);
		anode->maxEdges=10;
	}
	/*we dynamically allocate space for the edges */
	if(anode->maxEdges==anode->numberOfEdges){
		anode->maxEdges=anode->maxEdges*2;
		RESIZEARRAY(anode->edgeList,anode->maxEdges);
	}
	*(anode->edgeList+anode->numberOfEdges)=edgeId;
	*(anode->weigthsList+anode->numberOfEdges)=edgeWeight;
	anode->numberOfEdges++;
}

void pushDistanceAndNode(BUCKETELEMENT *bucketarray, int nodeId, double distance){

	if(bucketarray->nodeIds==NULL){
		NEWARRAY(bucketarray->nodeIds,5);
		NEWARRAY(bucketarray->tentDist,5);
		bucketarray->maxNodes=5;
	}

	/*we dynamically allocate space for the edges */
	if(bucketarray->maxNodes==bucketarray->numberOfNodes){
		bucketarray->maxNodes=bucketarray->maxNodes*2;
		RESIZEARRAY(bucketarray->nodeIds,bucketarray->maxNodes);
		RESIZEARRAY(bucketarray->tentDist,bucketarray->maxNodes);
	}
	*(bucketarray->nodeIds+bucketarray->numberOfNodes)=nodeId;
	*(bucketarray->tentDist+bucketarray->numberOfNodes)=distance;
	bucketarray->numberOfNodes++;
}

void popDistanceAndNode(BUCKETELEMENT *bucketarray, int nodeId, double distance){

	if(bucketarray->nodeIds==NULL){
		return;
	}

	int* tempNodeIds = calloc((bucketarray->maxNodes - 1), sizeof(int)); // allocate an array with a size 1 less than the current one
	double* tempDistances = calloc((bucketarray->maxNodes - 1), sizeof(double)); // allocate an array with a size 1 less than the current one

	int i=0;
	int j=0;
	for(i=0; i<bucketarray->numberOfNodes; i++){
		if(*(bucketarray->nodeIds+i)!=nodeId){

			*(tempNodeIds+j)=nodeId;
			*(tempDistances+j)=distance;
			j++;
		}
	}
	memcpy(bucketarray->nodeIds, tempNodeIds, ((bucketarray->maxNodes - 1)) * sizeof(int));
	memcpy(bucketarray->tentDist, tempDistances, ((bucketarray->maxNodes - 1)) * sizeof(int));
	bucketarray->numberOfNodes--;

	int N = bucketarray->numberOfNodes;
	int totalSize = bucketarray->maxNodes;

	if (N > 0 && N == totalSize/4){
		bucketarray->maxNodes = totalSize/2;
	}

}


void readEdges(NODE* graph,const char *line){
	char *token;
	short counter=0;
	int startNode,endNode;
	double weight;
	const char s[2] = " ";
	int startNodeId,endNodeId;
	/* get the first token */
	token = strtok(line, s);

	/* walk through other tokens */
	 while( token != NULL )
	 {
		if(counter==0){
			startNode=atoi(token);
		}else if(counter==1){
			endNode=atoi(token);
		}else{
			weight=atof(token);
		}
	    token = strtok(NULL, s);
	    counter++;
	 }
	 pushEdge(graph+startNode-1,endNode,weight);
}

int 
main(int argc, char *argv[])
{
	int			my_rank;		/* rank of process */
	int			num_procs=0;		/* number of processes */
	int			source;			/* rank of sender */
	int			dest = 0;		/* rank of receiver */
	int			tag = 0;		/* tag for messages */
	char		message[100];	/* storage for message */
	MPI_Status	status ;		/* return status for receive */
	char		*nameOfGraphFile;
	char line[80];
	long elapsed_seconds;
	int sizeOfScatteredList;
	FILE *fr;            /* declare the file pointer */
	int sizeOfGraph;
	NODE *graph;
	int *ind;
	int *bcastBuff;
	int maxSize;
	NODE *localScatteredReceiveList;
	NODE *globalScatteredList;

	/* start up MPI */
	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 


	typedef struct NODE{
		int nodeId;//the id of the node
		int *edgeList;//the array containing the adj node ids
		double *weigthsList;//the weight to each edge
		int maxEdges;//to dynamically allovate new size for the edge list
		int numberOfEdges;
	} NODE;


	/* Declaration of new datatype */
	const int nitems=5;
		int          blocklengths[3] = {1,1,1,1,1};
		MPI_Datatype types[3] = {MPI_INT, MPI_INT,MPI_DOUBLE,MPI_INT,MPI_INT};
		MPI_Datatype mpi_local_counts_type;
		MPI_Aint     offsets[5];
		MPI_Aint 	 addr[3];
		MPI_Get_address(globalScatteredList, &addr[0]);
		MPI_Get_address(&globalScatteredList->nodeId, &addr[1]);
		MPI_Get_address(&globalScatteredList->edgeList, &addr[2]);
		MPI_Get_address(&globalScatteredList->weigthsList, &addr[3]);
		MPI_Get_address(&globalScatteredList->maxEdges, &addr[4]);
		MPI_Get_address(&globalScatteredList->numberOfEdges, &addr[5]);

		offsets[0] = addr[1] - addr[0];
		offsets[1] = addr[2] - addr[0];
		offsets[2] = addr[3] - addr[0];
		offsets[3] = addr[4] - addr[0];
		offsets[4] = addr[5] - addr[0];
		MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_local_counts_type);
		MPI_Type_commit(&mpi_local_counts_type);

	NEWARRAY(ind,sizeOfGraph-1);
	if(my_rank==0){
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
				 graph=New_grapg(sizeOfGraph);
				 if(atoi(line)!=-1){
					 readEdges(&graph,line);
				 }
			}
		}
		srand(clock(NULL));
		int counter=0;
		int countsPerProcess[num_procs]={0};
		for(counter=0;counter<sizeOfGraph;counter++){
			ind[counter]=get_randomNumber(0,num_procs);//generate the randomly assignments of nodes to processes
			countsPerProcess[ind[counter]]++;
		}
		maxSize=0;
		for(counter=0;counter<num_procs;counter++){
			if (maxSize<countsPerProcess[counter]){
				maxSize=countsPerProcess[counter];
			}
		}
		sizeOfScatteredList= maxSize*num_procs;
		NEWARRAY(globalScatteredList,sizeOfScatteredList);
		int tempCountsPerProcess[num_procs]={0};///te
		for(counter=0;counter<sizeOfScatteredList;counter++){
			globalScatteredList+maxSize*(ind[counter]-1)=graph+tempCountsPerProcess[ind[counter]];
			tempCountsPerProcess[ind[counter]]++;
		}
	}
	MPI_Bcast(ind,sizeOfGraph-1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(sizeOfScatteredList,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(maxSize,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Scatter(globalScatteredList,maxSize,mpi_local_counts_type,localScatteredReceiveList,maxSize,mpi_local_counts_type,0,MPI_COMM_WORLD);

	///send the buckets ??or send the nodes??
	if(sizeOfGraph%num_procs==0){
		localScatteredReceiveList=New_grapg(sizeOfGraph/num_procs);
	}else{
		//to do
		printf("sizeOfGraph%num_procs==0 ");
		exit(1);
	}


	int ii=0;
	int localCounter=0;
	for(int ii=1;ii<num_procs;ii++){
		if(my_rank==ind[ii]){

		}
	}



	/* shut down MPI */
	MPI_Finalize(); 

	return 0;
}
