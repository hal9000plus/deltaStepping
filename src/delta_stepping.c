/*
 ============================================================================
 Name        : delta_stepping.c
 Author      : stam
 Version     :
 Copyright   : Your copyright notice
 Description : Calculate Pi in MPI
 ============================================================================
 */
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "mpi.h"
#include <time.h>
#include "limits.h"
#include <stddef.h>
#define NEWARRAY(x, n) do { (x) = calloc((n), sizeof *(x)); } while (0)
#define RESIZEARRAY(x, n) do { (x) = realloc((x), sizeof *(x) * (n) ); } while (0)
float precision = 0.00001;

//the graph will be represented as an array of nodes
typedef struct NODE{
	int nodeId;//the id of the node
	int *edgeList;//the array containing the adj node ids
	double *weigthsList;//the weight to each edge
	int maxEdges;//to dynamically allovate new size for the edge list
	int numberOfEdges;
} NODE;

typedef struct BUCKETELEMENT{
	double* tentDist;//-1.0 denotes infinite distance  //-2.0 just initialized
	int* nodeIds;
	int* bucket2Index;//helper to find the element in constant time
	int maxNodes;
	int numberOfNodes;
	int index;
} BUCKETELEMENT;

typedef struct REQUESTS{
	int *targetId;
	int numberOfelements;
	double *weightSourceToTarget;
	double tendDistanceOfSource;
}REQUESTS;

NODE* New_grapg(int numberOfNodes){
	NODE *newGraph= calloc(numberOfNodes, sizeof(NODE));
	    assert(newGraph != NULL);
	int i=0;
	for(i=0;i<numberOfNodes;i++){
		newGraph->numberOfEdges=0;
	}

	return newGraph;
}

BUCKETELEMENT* New_Bucket(int numberOfBuckets){
	BUCKETELEMENT *newBucket = calloc(numberOfBuckets, sizeof(BUCKETELEMENT));
		assert(newBucket!=NULL);
	int i=0;
	for(i=0;i<numberOfBuckets;i++){
		(newBucket+i)->index=-2;//-2 just initialized
		(newBucket+i)->numberOfNodes=0;
	}
	return newBucket;
}

BUCKETELEMENT* resizeBucket(BUCKETELEMENT* bucketArray,int newSize,int oldSize){
	RESIZEARRAY(bucketArray,newSize);
		assert(bucketArray!=NULL);
	int i=0;
	if(oldSize<newSize){
		for(i=oldSize;i<newSize;i++){
			(bucketArray+i)->index=-2;//-2 just initialized
			(bucketArray+i)->numberOfNodes=0;
		}
	}
	return bucketArray;
}

void buildRequests(NODE aNode,REQUESTS *lightEdgesRequests,double delta,double tendDistSource){
	int neigboursSize=aNode.numberOfEdges;
	int kk=0;
	int innerCounter=0;
	if(lightEdgesRequests->targetId==NULL||lightEdgesRequests->weightSourceToTarget==NULL){
		NEWARRAY(lightEdgesRequests->targetId,neigboursSize);//at max is neigboursSize
		NEWARRAY(lightEdgesRequests->weightSourceToTarget,neigboursSize);
		lightEdgesRequests->numberOfelements=0;
	}else{
		RESIZE(lightEdgesRequests->targetId,lightEdgesRequests->numberOfelements+neigboursSize);
		RESIZE(lightEdgesRequests->weightSourceToTarget,lightEdgesRequests->numberOfelements+neigboursSize);
	}
	lightEdgesRequests->tendDistanceOfSource=tendDistSource;
	for(kk=0 ; kk < neigboursSize ; kk++ ){
		if( *((aNode.weigthsList)+kk) <= delta ){
			*(lightEdgesRequests->(targetId+innerCounter+lightEdgesRequests->numberOfelements))=*(aNode.edgeList+kk);
			*(lightEdgesRequests->(weightSourceToTarget+innerCounter+lightEdgesRequests->numberOfelements)=*((aNode.weigthsList)+kk);
			lightEdgesRequests->numberOfelements++;
			innerCounter++;
		}
	}
}

void buildHeavyRequests(NODE aNode,REQUESTS *heaveEdgesRequests,double delta,double tendDistSource){
		int neigboursSize=aNode.numberOfEdges;
		int kk=0;
		int innerCounter=0;
		if(heaveEdgesRequests->targetId==NULL||heaveEdgesRequests->weightSourceToTarget==NULL){
			NEWARRAY(heaveEdgesRequests->targetId,neigboursSize);//at max is neigboursSize
			NEWARRAY(heaveEdgesRequests->weightSourceToTarget,neigboursSize);
			heaveEdgesRequests->numberOfelements=0;
		}else{
			RESIZE(heaveEdgesRequests->targetId,heaveEdgesRequests->numberOfelements+neigboursSize);
			RESIZE(heaveEdgesRequests->weightSourceToTarget,heaveEdgesRequests->numberOfelements+neigboursSize);
		}
		heaveEdgesRequests->tendDistanceOfSource=tendDistSource;
		for(kk=0 ; kk < neigboursSize ; kk++ ){
			if( *((aNode.weigthsList)+kk) >= delta ){
				*(lightEdgesRequests->(targetId+innerCounter+lightEdgesRequests->numberOfelements))=*(aNode.edgeList+kk);
				*(lightEdgesRequests->(weightSourceToTarget+innerCounter+lightEdgesRequests->numberOfelements)=*((aNode.weigthsList)+kk);
				lightEdgesRequests->(numberOfelements++);
				innerCounter++;
			}
		}
	}


int isBucketEmpty(BUCKETELEMENT *aBucketElement,int *removedNodes,int size){
	if(aBucketElement->numberOfNodes==0 || aBucketElement==NULL){
		return 1;
	}else{
		int i=0;
		int notFound=0;
		for(i=0;i<aBucketElement->numberOfNodes;i++){
			if(!containtsElem(aBucketElement->nodeIds,removedNodes,*removedNodes)){
				notFound=1;
				break;
			}
		}
		return !notFound;
	}
}
int isBucketEmpty(BUCKETELEMENT *aBucketElement){
	if(aBucketElement->numberOfNodes==0 || aBucketElement==NULL){
		return 1;
	}else{
		return 0;
	}
}

int containtsElem(int searchNode,int *searchlist,int size){
	int i;
	for(i=1;i<size;i++){
		if(searchlist[i]==searchNode){
			return 1;
		}else{
			return 0;
		}
	}
}

int pushToRemovedNodes(int nodeId,int pos,int arr,int *sizeOfRemoveArray){
	if(*(arr)==*sizeOfRemoveArray-1){
		*sizeOfRemoveArray=2*(*sizeOfRemoveArray);
		RESIZEARRAY(arr,sizeOfRemoveArray);
	}else{
		*(arr+pos)=nodeId;
		(*arr)++;
	}
}

void findnextNonEmptyBucket(BUCKETELEMENT *localBucketArray,int currentBucketIndex,int *localMinBuckInd,int localBucketSize){
	int ii=0;
	int found=0;
	for(ii=currentBucketIndex ; ii < localBucketSize ;ii++){
		if( !isBucketEmpty((localBucketArray+ii)) ){
			found=1;
			break;
		}
	}
	if(found){
		*localMinBuckInd=ii;
	}else{
		*localMinBuckInd=0;///only the infiniteBucket is Non Empty
	}
}



/* rangefrom inclusive end not inclusive */
int get_randomNumber(int rangefrom,int rangeTo)
{
	return (rand() % rangeTo + rangefrom);
}

void pushEdge(NODE *anode,int edgeId,double edgeWeight){
	if(anode->edgeList==NULL){
		NEWARRAY(anode->edgeList,10);
		NEWARRAY(anode->weigthsList,10);
		anode->maxEdges=10;
	}
	/*we dynamically allocate space for the edges */
	if(anode->maxEdges==anode->numberOfEdges){
		anode->maxEdges=anode->maxEdges*2;
		RESIZEARRAY((anode->edgeList),anode->maxEdges);
		RESIZEARRAY((anode->weigthsList),anode->maxEdges);
	}
	*(anode->edgeList+anode->numberOfEdges)=edgeId;
	*(anode->weigthsList+anode->numberOfEdges)=edgeWeight;
	anode->numberOfEdges++;
}

void pushDistanceAndNode(BUCKETELEMENT *bucketarray, int nodeId, double distance,int index,int *localBucketSize){

	if(bucketarray->nodeIds==NULL){
		NEWARRAY(bucketarray->nodeIds,5);
		NEWARRAY(bucketarray->tentDist,5);
		NEWARRAY(bucketarray->bucket2Index,5);
		bucketarray->maxNodes=5;
	}

	/*we dynamically allocate space for the edges */
	if(bucketarray->maxNodes==bucketarray->numberOfNodes){
		bucketarray->maxNodes=bucketarray->maxNodes*2;
		RESIZEARRAY((bucketarray->nodeIds),bucketarray->maxNodes);
		RESIZEARRAY((bucketarray->tentDist),bucketarray->maxNodes);
		RESIZEARRAY((bucketarray->bucket2Index),bucketarray->maxNodes);
	}
	*(bucketarray->nodeIds+bucketarray->numberOfNodes)=nodeId;
	*(bucketarray->tentDist+bucketarray->numberOfNodes)=distance;
	*(bucketarray->bucket2Index+bucketarray->numberOfNodes)=index;
	bucketarray->numberOfNodes++;
}

void popDistanceAndNode(BUCKETELEMENT *bucketarray, int nodeId){
	if(bucketarray->nodeIds==NULL){
		return;
	}
	short foundFlag=0;
	int* tempNodeIds = calloc((bucketarray->maxNodes - 1), sizeof(int)); // allocate an array with a size 1 less than the current one
	double* tempDistances = calloc((bucketarray->maxNodes - 1), sizeof(double)); // allocate an array with a size 1 less than the current one
	int* tempBucket2Index = calloc((bucketarray->maxNodes - 1), sizeof(double));
	int i=0;
	int j=0;
	for(i=0; i<bucketarray->numberOfNodes; i++){
		if(*(bucketarray->nodeIds+i)!=nodeId){
			*(tempNodeIds+j)=*(bucketarray->nodeIds+i);
			*(tempDistances+j)=*(bucketarray->tentDist+i);
			*(tempBucket2Index+j)=*(bucketarray->bucket2Index+i);
			j++;
		}else{
			foundFlag=1;
		}
	}
	if(foundFlag){
		memcpy(bucketarray->nodeIds, tempNodeIds, ((bucketarray->maxNodes - 1)) * sizeof(int));
		memcpy(bucketarray->tentDist, tempDistances, ((bucketarray->maxNodes - 1)) * sizeof(int));
		bucketarray->numberOfNodes--;

		int N = bucketarray->numberOfNodes;
		int totalSize = bucketarray->maxNodes;

		if (N > 0 && N <= totalSize/4){
			bucketarray->maxNodes = totalSize/2;
			RESIZEARRAY((bucketarray->nodeIds),bucketarray->maxNodes);
			RESIZEARRAY((bucketarray->tentDist),bucketarray->maxNodes);
		}
	}
}

double findTargetDistance(BUCKETELEMENT *locBucketArray,int nodeId,int sizeOfBucketArray,int *targetNodeBucketIndex,int *targetNodeIndexInBucket){
	int ii=0;
	int ee=0;
	for(ii=0; ii<sizeOfBucketArray ;ii++ ){
		for(ee=0 ; ee < (locBucketArray+ii)->numberOfNodes ; ee++){
			int currentNode=*((locBucketArray+ii)->nodeIds+ee);
			if(currentNode==nodeId){
				*targetNodeBucketIndex=ii;
				*targetNodeIndexInBucket=ee;
				return (locBucketArray+ii)->tentDist;
			}
		}
	}

}

int findNextBucketIndex(int *globalMinBuffer,int size){
	int ii=0;
	int tempMin=size;
	int allmin=1;
	for(ii=0;ii < size; ii++){
		if(*globalMinBuffer!=0 && (*globalMinBuffer)<tempMin ){
			tempMin=*globalMinBuffer;
			allmin=0;
		}
	}
	if(allmin){
		return 0;
	}else{
		return tempMin;
	}
}



void readEdges(NODE* graph,const char *line){
	char *token;
	short counter=0;
	int startNode,endNode;
	double weight;
	const char s[2] = " ";
	int startNodeId,endNodeId;
	//printf("reading line %s \n",line);
	/* get the first token */
	token = strtok((const char *)line, s);

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


	/* start up MPI */
	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 

	/*declarations and initializations*/
	int tempCountsPerProcess[num_procs-1];
	char		*nameOfGraphFile;
	char line[80];
	long elapsed_seconds;
	int sizeOfScatteredList;
	FILE *fr;            /* declare the file pointer */
	int *sizeOfGraph = NULL;
	NODE *graph= NULL;
	int *ind= NULL;
	int maxSize;
	NODE *localScatteredReceiveList;
	NODE *globalScatteredList;
	BUCKETELEMENT *localBucketArray;
	NEWARRAY(sizeOfGraph,1);
	memset( tempCountsPerProcess, 0, num_procs*sizeof(int) );//for c99 compiler
	int sourceNodeId;
	int currentBucketIndex;
	int localMinBuckInd;
	int *globalMinBuffer;
	REQUESTS *requestsArray;
	int *removed;
	double delta;
	int *gatherCountsSendBuffer=NULL;
	int *gatherCountsRcvBuffer=NULL;
	int gatherMax=0;
	int localBucketSize=0;
	int terminateFlag=0;
	int stillWorking;
	int *stillWorkingArray=NULL;
	NEWARRAY(stillWorkingArray,num_procs-1);
	/*end of declarations and initializations*/

	if(my_rank==0){
		if(argc < 2){
			printf("No Arguments \n Arguments:\n name of graph file \nintial delta value \n");
			exit(1);
		}
		nameOfGraphFile=argv[1];
		delta=atof(argv[2]);
		sourceNodeId=atoi(argv[3]);
		if(nameOfGraphFile!=NULL){
			fr = fopen (nameOfGraphFile, "r");  /* open the file for reading */
			   /* elapsed.dta is the name of the file */
			   /* "rt" means open the file for reading text */
			/* read the graph */
			if(fr==NULL){
				 printf("error: the file does not exists \n");
				 exit(1);
			}
			fgets(line, 80, fr);
			sscanf (line, "%ld", &elapsed_seconds);
			*sizeOfGraph=atoi(line);
			printf("first line %s \n",line);
			printf("the size of the graph %d \n",*sizeOfGraph);
			while(fgets(line, 80, fr) != NULL)
			{
				 /* get a line, up to 80 chars from fr.  done if NULL */
				 sscanf (line, "%ld", &elapsed_seconds);
				 /* convert the string to a double int */
				 if(*sizeOfGraph<1){
					 printf("error: The size of the graph is  < 1 \n");
					 exit(1);
				 }
				 graph=New_grapg(*sizeOfGraph);
				 if(atoi(line)!=-1){
					 readEdges(graph,line);
				 }
			}
		}

		srand(clock());
		int counter=0;
		int countsPerProcess[num_procs-1];
		memset( countsPerProcess, 0, (num_procs-1)*sizeof(int) );//for c99 compiler
		NEWARRAY(ind,*(sizeOfGraph)-1);
		for(counter=0;counter<*sizeOfGraph;counter++){
			ind[counter]=get_randomNumber(0,num_procs);//generate the randomly assignments of nodes to processes u
			countsPerProcess[ind[counter]]++;
		}
		maxSize=0;
		for(counter=0;counter<num_procs;counter++){
			if (maxSize<countsPerProcess[counter]){
				maxSize=countsPerProcess[counter];
			}
		}
		sizeOfScatteredList= maxSize*num_procs;
		globalScatteredList=New_grapg(sizeOfScatteredList);
		for(counter=0;counter<*sizeOfGraph;counter++){
			*(globalScatteredList+maxSize*(ind[counter])+tempCountsPerProcess[(ind[counter])])=*(graph+counter);
			tempCountsPerProcess[(ind[counter])]++;
		}
		printf("in here \n");
	}
	/* Declaration of new datatype */
	const int nitems=5;
	int          blocklengths[5] = {1,1,1,1,1};
	MPI_Datatype types[5] = {MPI_INT, MPI_INT,MPI_DOUBLE,MPI_INT,MPI_INT};
	MPI_Datatype mpi_local_counts_type;
	MPI_Aint     offsets[5];
	MPI_Aint 	 addr[5];
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
	/* end of new datatype */

	/* Declaration of new datatype */
	const int nitems2=3;
	int          blocklengths2[3] = {1,1,1};
	MPI_Datatype types2[3] = {MPI_INT, MPI_INT,MPI_DOUBLE};
	MPI_Datatype mpi_request_type;
	MPI_Aint     offsets2[3];
	MPI_Aint 	 addr2[3];
	MPI_Get_address(requestsArray, &addr2[0]);
	MPI_Get_address(&requestsArray->targetId, &addr2[1]);
	MPI_Get_address(&requestsArray->processIdOfTarget, &addr2[2]);
	MPI_Get_address(&requestsArray->weightSourceToTarget, &addr2[3]);
	offsets2[0] = addr2[1] - addr2[0];
	offsets2[1] = addr2[2] - addr2[0];
	offsets2[2] = addr2[3] - addr2[0];
	MPI_Type_create_struct(nitems2, blocklengths2, offsets2, types2, &mpi_request_type);
	MPI_Type_commit(&mpi_request_type);
		/* end of new datatype */
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(sizeOfGraph,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
		assert(sizeOfGraph != NULL);
	if(my_rank!=0){
		NEWARRAY(ind,*(sizeOfGraph)-1);
	}
		assert(ind  != NULL);
	MPI_Bcast(tempCountsPerProcess,num_procs-1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(ind,(*sizeOfGraph)-1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&delta,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sizeOfScatteredList,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&maxSize,1,MPI_INT,0,MPI_COMM_WORLD);
	localScatteredReceiveList=New_grapg(maxSize);
	MPI_Scatter(globalScatteredList,maxSize,mpi_local_counts_type,localScatteredReceiveList,maxSize,mpi_local_counts_type,0,MPI_COMM_WORLD);
	localBucketArray=New_Bucket(2);//initial we start with only the 0 bucket and the infinite distance bucket
	localBucketSize=2;
	(localBucketArray)->index=-1;//go to the bucket that holds infinite
	(localBucketArray)->numberOfNodes=0;
	(localBucketArray+1)->index=0;//go to the bucket that holds source
	(localBucketArray+1)->numberOfNodes=0;

	/* fill the all the local nodes to the infinite bucket*/
	int ll=0;
	for(ll=0;ll<tempCountsPerProcess[my_rank];ll++){
		if(((localScatteredReceiveList+ll)->nodeId)==sourceNodeId){
			printf("found source node on process %d",my_rank);
			pushDistanceAndNode(localBucketArray+1,(localScatteredReceiveList+ll)->nodeId,0.0);
		}else{
			pushDistanceAndNode(localBucketArray,(localScatteredReceiveList+ll)->nodeId,-1.0);
		}
	}
	currentBucketIndex=1;//the infinite bucket is the bucket array element 0;
	/*create requests */

	while(terminateFlag){

		int *removedNodesList;//the 0 index is used for the size of the elements stored in array
		int sizeOfremoveArray=localBucketArray->numberOfNodes+1;
		NEWARRAY(removedNodesList,sizeOfremoveArray);
		*removedNodesList=0;
		int isHeavyEgdesProcessed=-1;
		int isBucketEmptyFlag=isBucketEmpty( (localBucketArray+currentBucketIndex),removedNodesList );

		while( !isBucketEmptyFlag || !isHeavyEgdesProcessed  ){

			/*allocate variables */
			REQUESTS *requestsToSendList;
			REQUESTS *requestsToReceive;
			NEWARRAY(requestsToSendList,1);
			requestsToSendList->targetId=NULL;
			requestsToSendList->weightSourceToTarget=NULL;
			NEWARRAY(requestsToReceive,num_procs-1);///
			requestsToReceive->targetId=NULL;
			requestsToReceive->weightSourceToTarget=NULL;
			/*end of allocate variables */

			/* iterate all the nodes the bucket */
			for(ll=0;ll<(localBucketArray+currentBucketIndex)->numberOfNodes;ll++){
				//iterate all the edges
				int indexToInsert=0;
				int nodeToProcess=*((localBucketArray+currentBucketIndex)->nodeIds+ll);
				int indexInScatteredList=*((localBucketArray+currentBucketIndex)->bucket2Index+ll);
				double sourceTentDist=*((localBucketArray+currentBucketIndex)->tentDist+ll);
				if(!containtsElem(nodeToProcess,removedNodesList,*removedNodesList)){//if itis not inside the removed list
					buildRequests(nodeToProcess,requestsToSendList,delta,sourceTentDist);
					pushToRemovedNodes(nodeToProcess,indexToInsert+1,removedNodesList,&sizeOfremoveArray);
					indexToInsert++;
				}else if(isHeavyEgdesProcessed==0){
					buildHeavyRequests(nodeToProcess,requestsToSendList,delta,sourceTentDist);
				}
			}
			/*end of iterate all the nodes the bucket */

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allgather(requestsToSendList,1,mpi_request_type,requestsToReceive,1,mpi_request_type,MPI_COMM_WORLD);

			int tempOuterCounter=0;
			int tempInnerCounter=0;

			/* all processes perform Relaxations */
				for( tempOuterCounter=0 ; tempOuterCounter < num_procs-1; tempOuterCounter++){
					for(tempInnerCounter=0 ; tempInnerCounter < (requestsToReceive+tempOuterCounter)->numberOfelements ; tempInnerCounter){
						if((requestsToReceive+tempOuterCounter)->numberOfelements!=0){
							int currentTargeNodeId=*( (requestsToReceive+tempOuterCounter)->targetId+tempInnerCounter )
								if( ind[currentTargeNodeId]==my_rank ){//if true the the node must be someware in the bucket O(n) time for searching ///or improvement we could use a hashSet ..
								int targetNodeBucketIndex=0;
								int targeNodeIndexInsideBucket=0;
								double tentDistSource=(requestsToReceive+tempOuterCounter)->tendDistanceOfSource;
								double tentDistTarget=findTargetDistance(localBucketArray,currentTargeNodeId,localBucketSize,&targetNodeBucketIndex,&targeNodeIndexInsideBucket);
								double sourceToTargetWeight=*((requestsToReceive+tempOuterCounter)->weightSourceToTarget+tempInnerCounter);

								if( (tentDistSource+sourceToTargetWeight < tentDistTarget) || ( tentDistTarget == ( -0.1 ) ) ){
									popDistanceAndNode(localBucketArray+targetNodeBucketIndex);
									int arrayIndexToInsertNeWdDist=( (tentDistSource+sourceToTargetWeight)/(int)delta)+1;
									if( (arrayIndexToInsertNeWdDist+1) >= localBucketSize){
										RESIZEARRAY(localBucketArray,(arrayIndexToInsertNeWdDist+2));
										localBucketSize=(arrayIndexToInsertNeWdDist+2);
										pushDistanceAndNode(localBucketArray+arrayIndexToInsertNeWdDist+1,currentTargeNodeId,tentDistSource+sourceToTargetWeight,arrayIndexToInsertNeWdDist)//change the bucket the node is
									}
								}
							}
						}
					}
				}

				/* deallocate space for dynamically allocated variables */
				free(requestsToSendList->targetId);
				free(requestsToSendList->weightSourceToTarget);
				free(requestsToReceive->weightSourceToTarget);
				free(requestsToReceive->targetId);
				free(requestsToSendList);
				free(requestsToReceive);
				isBucketEmptyFlag=isBucketEmpty( (localBucketArray+currentBucketIndex),removedNodesList );
				if(isBucketEmptyFlag){
					isHeavyEgdesProcessed++;
				}
				MPI_Barrier(MPI_COMM_WORLD);
				stillWorking=(!isBucketEmptyFlag || !isHeavyEgdesProcessed);
				MPI_Allgather(&stillWorking,1,MPI_INT,stillWorkingArray,1,MPI_INT,MPI_COMM_WORLD);

		}

		free(removedNodesList);
		int ii=0;
		terminateFlag=1;
		for(ii=0;ii<num_procs-1;ii++){
			if(*(stillWorkingArray+ii)!=0){
				terminateFlag=0;
				break;
			}
		}

		if(isBucketEmptyFlag && !terminateFlag && !stillWorking){//other processes still work  but current has no more work
			/*allocate variables */
			REQUESTS *requestsToSendList;
			REQUESTS *requestsToReceive;
			NEWARRAY(requestsToSendList,1);
			requestsToSendList->targetId=NULL;
			requestsToSendList->weightSourceToTarget=NULL;
			requestsToSendList->numberOfelements=0;
			NEWARRAY(requestsToReceive,num_procs-1);///
			requestsToReceive->targetId=NULL;
			requestsToReceive->weightSourceToTarget=NULL;
			/*end of allocate variables */

			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allgather(requestsToSendList,1,mpi_request_type,requestsToReceive,1,mpi_request_type,MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			stillWorking=(!isBucketEmptyFlag || !isHeavyEgdesProcessed);

			free(requestsToSendList->targetId);
			free(requestsToSendList->weightSourceToTarget);
			free(requestsToReceive->weightSourceToTarget);
			free(requestsToReceive->targetId);
			free(requestsToSendList);
			free(requestsToReceive);
			/* deallocate space for dynamically allocated variables */

		}else if(isBucketEmptyFlag && terminateFlag && !stillWorking){//all processes finished processing their bucket check if there is a bucket to process

			findnextNonEmptyBucket(localBucketArray,currentBucketIndex,&localMinBuckInd);
			MPI_Allgather(&localMinBuckInd,1,MPI_INT,globalMinBuffer,1,MPI_INT,MPI_COMM_WORLD)
			currentBucketIndex=findNextBucketIndex(globalMinBuffer,num_procs-1);
			if(currentBucketIndex==0){
				terminateFlag=1;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Gather(&localMinBuckInd,1,MPI_INT,globalMinBuffer,1,MPI_INT,MPI_COMM_WORLD);


	/* shut down MPI */
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize(); 

	return 0;
}
