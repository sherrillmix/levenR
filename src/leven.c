//R CMD SHLIB leven.c
//#include <string>
//#include <algorithm>
//#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <stdint.h>

#define forR 1
#if forR==1
	#include <R.h>
	#define errorMessage(A,B) error(A)
	#define warningMessage(A) warning(A)
	//printf generates errors if called from threads. just suppress debug
	#define printMessage(...) Rprintf(__VA_ARGS__)
	//#define printMessage(...)
#else
	#define errorMessage(A,B) printf(A);exit(B)
	#define warningMessage(A) printf(A)
	#define printMessage(...) printf(__VA_ARGS__)
#endif

void combineAligns(char **queries, char **refs, char **out, int *nReadsIn, int *readLengths,int *maxOutLengthIn){
	int nReads=nReadsIn[0];
	int maxOutLength=maxOutLengthIn[0];
	int *readIndexes= malloc(sizeof(int)*(nReads));
	int outIndex=0;
	int maxRemaining=0;
	int ii;
	int anyGap=0;
	char thisBase;
	//get ready
	for(ii=0;ii<nReads;ii++){
		//find max
		if(maxRemaining<readLengths[ii])maxRemaining=readLengths[ii];
		//point pointers to start of read
		readIndexes[ii]=0;
	}
	while(maxRemaining>0){
		maxRemaining=0;
		anyGap=0;
		for(ii=0;ii<nReads;ii++){
			if(readIndexes[ii]<readLengths[ii]&&refs[ii][readIndexes[ii]]=='-'){
				anyGap=1;
				break;
			}
		}
		for(ii=0;ii<nReads;ii++){
			if(readIndexes[ii]>=readLengths[ii]){
				thisBase='-';
			}else if(!anyGap || refs[ii][readIndexes[ii]]=='-'){
				thisBase=queries[ii][readIndexes[ii]];
				readIndexes[ii]++;
			}else{
				thisBase='-';
			}
			out[ii][outIndex]=thisBase;
			if(maxRemaining<readLengths[ii]-readIndexes[ii])maxRemaining=readLengths[ii]-readIndexes[ii];
		}
		outIndex++;
		if(outIndex>maxOutLength){errorMessage("Out read length longer than anticipated. Adjust C code",11);}
	}
	//wrap up
	for(ii=0;ii<nReads;ii++){
		out[ii][outIndex]='\0';
	}

	free(readIndexes);
}

int countMismatch(char *s1, char *s2, int length, int s2Offset, int cutoff, int *weights) {
	int answer=0;
	for(unsigned int ii = 0; ii < length; ii++){
		if(s1[ii]!=s2[ii+s2Offset]) answer+=weights[ii];
		if(answer >= cutoff)return(answer);
	}
	return(answer);
}

void bestMismatch(int *ans, char **s1, char **s2, int *weights){
	unsigned int s1Length=strlen(s1[0]);
	unsigned int s2Length=strlen(s2[0]);
	unsigned int lastPos=s2Length-s1Length;
	int tmp;
	//INT_MAX doesn't appear to be defined
	ans[0]=99999;
	for(int ii = 0; ii <= lastPos; ii++){
		tmp=countMismatch(s1[0],s2[0],s1Length,ii,ans[0],weights);
		if(tmp < ans[0]){
			ans[0]=tmp;
			ans[1]=ii+1;
		}
	}
}

void levenAll(int *answer, char **s1, char **s2, int *homoLimit, int *prepend, int *append, int *debug, int *isAlign, char **align1, char **align2) {
	unsigned int ii,jj;//counters
	unsigned int cost_del = 1;
	unsigned int cost_ins = 1;
	const unsigned int cost_sub = 1;
	unsigned int n1 = strlen(s1[0]);
	unsigned int n2 = strlen(s2[0]);
	//printMessage("%d %d\n",n1,n2);
	unsigned int* endCol = malloc(sizeof(unsigned int)*(n1+1));
	unsigned int* lastRow = malloc(sizeof(unsigned int)*(n2+1));
	unsigned int* thisRow = malloc(sizeof(unsigned int)*(n2+1));
	unsigned int* cost_dels = malloc(sizeof(unsigned int)*(n1)); //store costs derived from homopolymers
	unsigned int* cost_inss = malloc(sizeof(unsigned int)*(n2)); //store costs derived from homopolymers
	unsigned int* tmpRow;
	unsigned int min;
	unsigned int **array;
	uint8_t **trace; //bit field: 1=goes down,2=goes right,4=goes diagonal,8=comes from up, 16=comes from left, 32=comes from diagonal, 64=on final path
	unsigned int coords[3],newCoords[2]; //for tracing back
	unsigned int d_del,d_ins,d_sub;
	if(isAlign[0]){
		//Make giant array
		array=malloc(sizeof(unsigned int *)*(n1+1));
		for(ii = 0;ii<n1+1;ii++)array[ii]=malloc(sizeof(unsigned int)*(n2+1));
		//Make traceback array for prettier alignments
		trace=malloc(sizeof(uint8_t *)*(n1+1));
		for(ii = 0;ii<n1+1;ii++)trace[ii]=malloc(sizeof(uint8_t)*(n2+1));
		for(ii = 0;ii<n1+1;ii++){
			for(jj = 0;jj<n2+1;jj++){
				trace[ii][jj]=0;
				if(ii==0){
					if(jj!=n2)trace[ii][jj] |=2;
					if(jj!=0)trace[ii][jj] |=16;
				}
				if(jj==0){
					if(ii!=n1)trace[ii][jj] |=1;
					if(ii!=0)trace[ii][jj] |=8;
				}
			}
		}
	}

	if(*homoLimit){
		for(ii = 0; ii < n1; ii++){
			cost_del=0;
			if(ii>=*homoLimit){
				for(unsigned int step = 1; step <= *homoLimit; step++){
					//printMessage("|%c %c|",s1[0][ii-1],s1[0][ii-1-step]);
					if(s1[0][ii]!=s1[0][ii-step]){
						cost_del = 1;
						break;
					}
				}
			}else cost_del=1;
			cost_dels[ii]=cost_del;
		}
		for(jj = 0; jj < n2; jj++){
			cost_ins=0;
			if(jj>=*homoLimit){
				for(unsigned int step = 1; step <= *homoLimit; step++){
					if(s2[0][jj]!=s2[0][jj-step]){
						cost_ins = 1;
						break;	
					}
				}
			}else cost_ins=1;
			cost_inss[jj]=cost_ins;
		}
	}else{
		for(ii = 0; ii < n1; ii++) cost_dels[ii]=1;
		for(jj = 0; jj < n2; jj++) cost_inss[jj]=1;
	}

	lastRow[0] = 0;
	//set first row to 0s if prepend on second string
	if(prepend[1]) for(jj = 1; jj <= n2; jj++) lastRow[jj] = 0;
	else for(jj = 1; jj <= n2; jj++) lastRow[jj] = lastRow[jj-1] + cost_inss[jj-1];
	if(isAlign[0]){
		for(jj=0;jj<=n2;jj++)array[0][jj]=lastRow[jj];
	}
	
	endCol[0]=lastRow[n2];

	for(ii = 1; ii <= n1; ii++){
		cost_del=cost_dels[ii-1];
		//Set first column to 0s if prepend on first string
		if(prepend[0]) thisRow[0] = lastRow[0];
		else thisRow[0] = lastRow[0] + cost_del;
		for(jj = 1; jj <= n2; jj++){
			//printMessage("S1[%d]:%c S2[%d]:%c\n",ii-1,s1[0][ii-1],jj-1,s2[0][jj-1]);
			cost_ins=cost_inss[jj-1];
			if(ii==n1&append[1])cost_ins=0;
			if(jj==n2&append[0])cost_del=0;
			d_del = lastRow[jj] + cost_del;
			d_ins = thisRow[jj-1] + cost_ins;
			//printMessage("%d %d",ii,jj);
			d_sub = lastRow[jj-1] + ( s1[0][ii-1] == s2[0][jj-1] ? 0 : cost_sub );
			thisRow[jj]=d_ins;
			if (d_del < thisRow[jj])thisRow[jj]=d_del;
			if (d_sub < thisRow[jj])thisRow[jj]=d_sub;	
			if(isAlign[0]){
				//printMessage("%d,%d:%d,%d,%d\n",ii,jj,d_del==thisRow[jj],d_ins==thisRow[jj],d_sub==thisRow[jj]);
				if(d_del==thisRow[jj]){
					trace[ii-1][jj]|= 1;
					trace[ii][jj]|=8;
				}
				if(d_ins==thisRow[jj]){
					trace[ii][jj-1]|= 2;
					trace[ii][jj]|=16;
				}
				if(d_sub==thisRow[jj]){
					trace[ii-1][jj-1]|= 4;
					trace[ii][jj]|=32;
				}
			}
		}
		//keep track of the final column for substringing
		endCol[ii]=thisRow[n2];
		//store this row if we're aligning (could store directly instead)
		if(isAlign[0]){
			for(jj=0;jj<=n2;jj++)array[ii][jj]=thisRow[jj];
		}
		//Switch the pointers around
		tmpRow = lastRow;
		lastRow = thisRow;
		thisRow = tmpRow;
		if(*debug){for(int printer = 0; printer <=n2; ++printer) printMessage("%d ",thisRow[printer]); printMessage("\n");}
	}
	if(*debug){for(int printer = 0; printer <=n2; ++printer) printMessage("%d ",lastRow[printer]); printMessage("\n");}
	*answer = lastRow[n2];
	if(isAlign[0]){
		//traceback to start storing valid paths, then descend following lowest values preferring diagonals
		//traceback
		unsigned int rightCoord=n2;
		unsigned int newRightCoord;
		unsigned int leftCoord=0;
		unsigned int newLeftCoord;
		trace[n1][n2] |=64;
		//trace bit field: 1=goes down,2=goes right,4=goes diagonal,8=comes from up, 16=comes from left, 32=comes from diagonal, 64=on final path
		//mark what nodes are on true path with a 64 in trace[ii][jj]
		for(ii=n1;ii+1>=0+1;ii--){//careful about wrapping around here
			newRightCoord=0;
			newLeftCoord=n2;
			for(jj=rightCoord;jj+1>=leftCoord+1;jj--){//careful about wrapping here
				if((trace[ii][jj] & 64)==0)continue;//current node not on true path so continue
				//if(debug*)printMessage("%d,%d=%d: %d-%d\ bits:%d\n",ii,jj,trace[ii][jj]&64,rightCoord,leftCoord,trace[ii][jj]);
				//up
				if(trace[ii][jj] & 8){
					if(jj>=newRightCoord)newRightCoord=jj;
					if(jj<=newLeftCoord)newLeftCoord=jj;
					trace[ii-1][jj] |= 64;
				}
				//left
				if(trace[ii][jj] & 16){
					trace[ii][jj-1] |= 64;
					if(leftCoord==jj)leftCoord--;
				}
				//diagonal
				if(trace[ii][jj] & 32){
					if(jj-1>=newRightCoord)newRightCoord=jj-1;
					if(jj-1<=newLeftCoord)newLeftCoord=jj-1;
					trace[ii-1][jj-1] |= 64;
				}
			}
			leftCoord=newLeftCoord;
			rightCoord=newRightCoord;
		}
		if(*debug){
			printMessage("Trace matrix\n");
			for(ii=0;ii<n1+1;ii++){
				for(jj=0;jj<n2+1;jj++){
					printMessage("%d ",trace[ii][jj]);
				}
				printMessage("\n");
			}
		}
		//walk down trace from top left to bottom right
		coords[0]=0;
		coords[1]=0;
		coords[2]=0;

		//trace bit field: 1=goes down,2=goes right,4=goes diagonal,8=comes from up, 16=comes from left, 32=comes from diagonal, 64=on final path
		while(coords[0]<n1||coords[1]<n2){
			//find min of three choices
			min=UINT_MAX;

			//go diagonal
			if(trace[coords[0]][coords[1]]&4 && trace[coords[0]+1][coords[1]+1]&64 && array[coords[0]+1][coords[1]+1]<min ){
				min=array[coords[0]+1][coords[1]+1];
				newCoords[0]=coords[0]+1;
				newCoords[1]=coords[1]+1;
			}
			//go down
			if(trace[coords[0]][coords[1]]&1 && trace[coords[0]+1][coords[1]]&64 && array[coords[0]+1][coords[1]]<min ){
				min=array[coords[0]+1][coords[1]];
				newCoords[0]=coords[0]+1;
				newCoords[1]=coords[1];
			}
			//go right
			if(trace[coords[0]][coords[1]]&2 && trace[coords[0]][coords[1]+1]&64 && array[coords[0]][coords[1]+1]<min ){
				min=array[coords[0]][coords[1]+1];
				newCoords[0]=coords[0];
				newCoords[1]=coords[1]+1;
			}
			if(*debug)printMessage("%d:%d,%d=>%d,%d",coords[2],coords[0],coords[1],newCoords[0],newCoords[1]);
			if(newCoords[0]==coords[0]&&coords[1]==newCoords[1])break;
			if(newCoords[0]==coords[0]){
				align1[0][coords[2]]='-';	
			}else{
				align1[0][coords[2]]=s1[0][coords[0]];
			}
			coords[0]=newCoords[0];
			if(newCoords[1]==coords[1]){
				align2[0][coords[2]]='-';	
			}else{
				align2[0][coords[2]]=s2[0][coords[1]];
			}
			coords[1]=newCoords[1];
			coords[2]++;
			if(*debug){
				align1[0][coords[2]]='\0';
				align2[0][coords[2]]='\0';
				printMessage("%s = %s\n",align1[0],align2[0]);
			}
		}
		align1[0][coords[2]]='\0';
		align2[0][coords[2]]='\0';
		if(*debug)printMessage("%s\n%s\n",align1[0],align2[0]);

		if(*debug){
			for(ii=0;ii<n1+1;ii++){
				for(jj=0;jj<n2+1;jj++){
					printMessage("%d ",array[ii][jj]);
				}
				printMessage("\n");
			}
		}
	}
	//printMessage(" %d ",*answer);
	free(endCol);
	free(lastRow);
	free(thisRow);
	free(cost_dels);
	free(cost_inss);
	if(isAlign[0]){
		for(ii = 0;ii<n1+1;ii++)free(array[ii]);
		free(array);
		for(ii = 0;ii<n1+1;ii++)free(trace[ii]);
		free(trace);
	}
}

struct levenArgs{
	int *answer;
	char **s1;
	char **s2;
	int *homoLimit;
	int *prepend;
	int *append;
	int *debug;
	int *isAlign;
	char **align1;
	char **align2;
	int *nStrings;
	int iRange[2];
	int jRange[2];
};

void *levenAllPar(void *levenArgs){
	struct levenArgs *args=(struct levenArgs *)levenArgs;
	if(*args->debug)printMessage("Thread started. iRange: %d-%d, %d-%d\n",args->iRange[0],args->iRange[1],args->jRange[0],args->jRange[1]);
	unsigned int ii, jj;
	int *ansPoint;
	char **align1Point, **align2Point;
	for(ii=args->iRange[0];ii<=args->iRange[1];ii++){
		for(jj=args->jRange[0];jj<=args->jRange[1];jj++){
			ansPoint=&(args->answer[ii+args->nStrings[0]*jj]);
			align1Point=&(args->align1[ii+args->nStrings[0]*jj]);
			align2Point=&(args->align2[ii+args->nStrings[0]*jj]);
			levenAll(ansPoint,&args->s1[ii],&args->s2[jj],args->homoLimit,args->prepend,args->append,args->debug,args->isAlign,align1Point,align2Point);
		}
	}
	if(*args->debug)printMessage("Thread ended. iRange: %d-%d\n",args->iRange[0],args->iRange[1]);
}
void parallelLeven(int *answer, char **s1, char **s2, int *nStrings, int *homoLimit, int *prepend, int *append, int *debug, int *nThread, int *isAlign, char **align1,char **align2){
	//just split on ii  for now
	unsigned int ii;//counter
	int nextI=0;//count up when dividing strings
	int nThreads=nThread[0];
	int stepSize;
	if(nStrings[0]<nThreads)nThreads=nStrings[0];
	//threads
	pthread_t *threads=(pthread_t *)malloc(sizeof(pthread_t)*nThreads);
	struct levenArgs **args=(struct levenArgs **)malloc(sizeof(struct levenArgs*)*nThreads);
	if(*debug)printMessage("Starting %d threads",nThreads);
	for(ii=0;ii<nThreads;ii++){
		args[ii]=(struct levenArgs *)malloc(sizeof(struct levenArgs));
		args[ii]->answer=answer; args[ii]->s1=s1;args[ii]->s2=s2;args[ii]->homoLimit=homoLimit;args[ii]->prepend=prepend;args[ii]->append=append;args[ii]->debug=debug;args[ii]->isAlign=isAlign;args[ii]->align1=align1;args[ii]->align2=align2;args[ii]->nStrings=nStrings;
		
		stepSize=(nStrings[0]-nextI)/(nThreads-ii); //integer division
		args[ii]->iRange[0]=nextI;
		args[ii]->iRange[1]=nextI+stepSize-1;
		if(*debug){printMessage("irange[%d]: %d-%d \n",ii,args[ii]->iRange[0],args[ii]->iRange[1]);}
		nextI+=stepSize;
		//not doing anything smart with jj now
		args[ii]->jRange[0]=0;
		args[ii]->jRange[1]=nStrings[1]-1;
		if(pthread_create(&threads[ii],NULL,levenAllPar,args[ii])){errorMessage("Couldn't create thread",9);}
	}
	for(ii=0;ii<nThreads;ii++){
		if(pthread_join(threads[ii],NULL)){errorMessage("Couldn't join thread",10);}
	}

	for(ii=0;ii<nThreads;ii++)free(args[ii]);
	free(args);
	free(threads);
}

