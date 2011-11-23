//R CMD SHLIB leven.c
//#include <string>
//#include <algorithm>
//#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#define forR 1
#if forR==1
	#include <R.h>
	#define errorMessage(A,B) error(A)
	#define warningMessage(A) warning(A)
#else
	#define errorMessage(A,B) printf(A);exit(B)
	#define warningMessage(A) printf(A)
#endif


int countMismatch(char *s1, char *s2, int length, int s2Offset, int cutoff, int *weights) {
	int answer=0;
	for(unsigned int i = 0; i < length; i++){
		if(s1[i]!=s2[i+s2Offset]) answer+=weights[i];
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
	for(int i = 0; i <= lastPos; i++){
		tmp=countMismatch(s1[0],s2[0],s1Length,i,ans[0],weights);
		if(tmp < ans[0]){
			ans[0]=tmp;
			ans[1]=i+1;
		}
	}
}

void revString(char *string, int nChar){
	if(nChar<2)return;
	int i;
	char tmp;
	for(i=0;i<nChar/2;i++){ //using integer division to round down
		tmp=string[i];
		string[i]=string[nChar-1-i];
		string[nChar-1-i]=tmp;
	}
}

void levenAll(int *answer, char **s1, char **s2, int *homoLimit, int *prepend, int *append, int *debug, int *isAlign, char **align) {
	unsigned int i,j;//counters
	unsigned int cost_del = 1;
	unsigned int cost_ins = 1;
	const unsigned int cost_sub = 1;
	unsigned int n1 = strlen(s1[0]);
	unsigned int n2 = strlen(s2[0]);
	//printf("%d %d\n",n1,n2);
	unsigned int* endCol = malloc(sizeof(unsigned int)*(n1+1));
	unsigned int* lastRow = malloc(sizeof(unsigned int)*(n2+1));
	unsigned int* thisRow = malloc(sizeof(unsigned int)*(n2+1));
	unsigned int* cost_dels = malloc(sizeof(unsigned int)*(n1)); //store costs derived from homopolymers
	unsigned int* cost_inss = malloc(sizeof(unsigned int)*(n2)); //store costs derived from homopolymers
	unsigned int* tmpRow;
	unsigned int min;
	unsigned int **array;
	unsigned int appendTracker; //keep track of location for alignment when we're doing ends free
	unsigned int coords[3],newCoords[2],upCoord,leftCoord; //for tracing back
	if(isAlign[0]){
		//Make giant array
		array=malloc(sizeof(unsigned int *)*(n1+1));
		for(i = 0;i<n1+1;i++)array[i]=malloc(sizeof(unsigned int)*(n2+1));
	}

	if(*homoLimit){
		for(i = 0; i < n1; i++){
			cost_del=0;
			if(i>=*homoLimit){
				for(unsigned int step = 1; step <= *homoLimit; step++){
					//printf("|%c %c|",s1[0][i-1],s1[0][i-1-step]);
					if(s1[0][i]!=s1[0][i-step]){
						cost_del = 1;
						break;
					}
				}
			}else cost_del=1;
			cost_dels[i]=cost_del;
		}
		for(j = 0; j < n2; j++){
			cost_ins=0;
			if(j>=*homoLimit){
				for(unsigned int step = 1; step <= *homoLimit; step++){
					if(s2[0][j]!=s2[0][j-step]){
						cost_ins = 1;
						break;	
					}
				}
			}else cost_ins=1;
			cost_inss[j]=cost_ins;
		}
	}else{
		for(i = 0; i < n1; i++) cost_dels[i]=1;
		for(j = 0; j < n2; j++) cost_inss[j]=1;
	}

	lastRow[0] = 0;
	//set first row to 0s if append ok on first string
	if(prepend[1]) for(j = 1; j <= n2; j++) lastRow[j] = 0;
	else for(j = 1; j <= n2; j++) lastRow[j] = lastRow[j-1] + cost_inss[j-1];
	if(isAlign[0]){
		for(j=0;j<=n2;j++)array[0][j]=lastRow[j];
	}
	
	endCol[0]=lastRow[n2];

	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	for(i = 1; i <= n1; i++){
		cost_del=cost_dels[i-1];
		//Set first column to 0s if prepend ok on first string
		if(prepend[0]) thisRow[0] = lastRow[0];
		else thisRow[0] = lastRow[0] + cost_del;
		for(j = 1; j <= n2; j++){
			//printf("S1[%d]:%c S2[%d]:%c\n",i-1,s1[0][i-1],j-1,s2[0][j-1]);
			cost_ins=cost_inss[j-1];
			unsigned int d_del = lastRow[j] + cost_del;
			unsigned int d_ins = thisRow[j-1] + cost_ins;
			//printf("%d %d",i,j);
			unsigned int d_sub = lastRow[j-1] + ( s1[0][i-1] == s2[0][j-1] ? 0 : cost_sub );
			thisRow[j]=d_ins;
			if (d_del < thisRow[j])thisRow[j]=d_del;
			if (d_sub < thisRow[j])thisRow[j]=d_sub;	
		}
		//keep track of the final column for substringing
		endCol[i]=thisRow[n2];
		//store this row if we're aligning (could store directly instead)
		if(isAlign[0]){
			for(j=0;j<=n2;j++)array[i][j]=thisRow[j];
		}
		//Switch the pointers around
		tmpRow = lastRow;
		lastRow = thisRow;
		thisRow = tmpRow;
		if(*debug){for(int printer = 0; printer <=n2; ++printer) printf("%d ",thisRow[printer]); printf("\n");}
	}
	if(*debug){for(int printer = 0; printer <=n2; ++printer) printf("%d ",lastRow[printer]); printf("\n");}
	*answer = lastRow[n2];
	if(append[1]){
		min=*answer;
		appendTracker=n2;
		for(j = 0; j <= n2; ++j ){
			if(lastRow[j]<min){
				min=lastRow[j];
				appendTracker=j;
			}
		}
		*answer = min;
	}
	if(append[0]){
		min=*answer;
		appendTracker=n2;
		if(*debug)printf("Last column: ");
		for(i = 0; i <= n1; ++i ){
			if(*debug)printf(" %d ",endCol[i]);
			if(endCol[i]<min){
				min=endCol[i];
				appendTracker=i;
			}
		}
		*answer = min;
	}
	if(isAlign[0]){
		coords[0]=n1;
		coords[1]=n2;
		coords[2]=0;
		if(append[0]){
			for(i=n1;i>appendTracker;i--){
				align[0][coords[2]]=s1[0][coords[0]];
				align[1][coords[2]]='-';	
				coords[0]--;
				coords[2]++;
				if(i==0)exit(1);//if this happens we're going to wrap around (shouldn't happen)
			}
		}
		if(append[1]){
			for(j=n2;j>appendTracker;j--){
				align[0][coords[2]]='-';	
				align[1][coords[2]]=s2[0][coords[1]];
				coords[1]--;
				coords[2]++;
				if(j==0)exit(1);//if this happens we're going to wrap around (shouldn't happen)
			}
		}
		//could reprogram to look for least gappy alignment?
		while(coords[0]!=0||coords[1]!=0){
			//deal with edge of array
			upCoord=coords[0]<1?0:coords[0]-1;
			leftCoord=coords[1]<1?0:coords[1]-1;
			//find min of three choices
			newCoords[0]=upCoord;newCoords[1]=leftCoord;
			min=array[upCoord][leftCoord];
			if(array[upCoord][coords[1]]<min){
				min=array[upCoord][coords[1]];
				newCoords[0]=upCoord;
				newCoords[1]=coords[1];
			}
			if(array[coords[0]][leftCoord]<min){
				min=array[coords[0]][leftCoord];
				newCoords[0]=coords[0];
				newCoords[1]=leftCoord;
			}
			if(newCoords[0]==coords[0]){
				align[0][coords[2]]='-';	
			}else{
				align[0][coords[2]]=s1[0][coords[0]-1];
				coords[0]--;
			}
			if(newCoords[1]==coords[1]){
				align[1][coords[2]]='-';	
			}else{
				align[1][coords[2]]=s2[0][coords[1]-1];
				coords[1]--;
			}
			//printf("New %d %d left %d up %d coords[2] %d\n",coords[0],coords[1],leftCoord,upCoord,coords[2]);
			coords[2]++;
		}
		align[0][coords[2]]='\0';
		align[1][coords[2]]='\0';
		//reverse these strings
		revString(align[0],coords[2]);
		revString(align[1],coords[2]);
		if(*debug)printf("%s\n%s\n",align[0],align[1]);

		if(*debug){
			for(i=0;i<n1+1;i++){
				for(j=0;j<n2+1;j++){
				printf("%d ",array[i][j]);
				}
				printf("\n");
			}
		}
			
	}
	//printf(" %d ",*answer);
	free(endCol);
	free(lastRow);
	free(thisRow);
	free(cost_dels);
	free(cost_inss);
	if(isAlign[0]){
		for(i = 0;i<n1+1;i++)free(array[i]);
		free(array);
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
	char **align;
	int *nStrings;
	int iRange[2];
	int jRange[2];
};

void *levenAllPar(void *levenArgs){
	struct levenArgs *args=(struct levenArgs *)levenArgs;
	if(*args->debug)printf("Thread started. iRange: %d-%d, %d-%d\n",args->iRange[0],args->iRange[1],args->jRange[0],args->jRange[1]);
	printf("Thread started. iRange: %d-%d, %d-%d\n",args->iRange[0],args->iRange[1],args->jRange[0],args->jRange[1]);
	unsigned int i, j;
	int *ansPoint;
	for(i=args->iRange[0];i<=args->iRange[1];i++){
		for(j=args->jRange[0];j<=args->jRange[1];j++){
			ansPoint=&(args->answer[i+args->nStrings[0]*j]);
			levenAll(ansPoint,&args->s1[i],&args->s2[j],args->homoLimit,args->prepend,args->append,args->debug,args->isAlign,args->align);
		}
	}
	if(*args->debug)printf("Thread ended. iRange: %d-%d\n",args->iRange[0],args->iRange[1]);
	printf("Thread ended. iRange: %d-%d\n",args->iRange[0],args->iRange[1]);
}
void parallelLeven(int *answer, char **s1, char **s2, int *nStrings, int *homoLimit, int *prepend, int *append, int *debug, int *nThread){
	//just split on i  for now
	unsigned int i;//counter
	int nextI=0;//count up when dividing strings
	int nThreads=nThread[0];
	int stepSize;
	if(nStrings[0]<nThreads)nThreads=nStrings[0];
	int isAlign[1]={0};
	char **align=NULL;
	//threads
	pthread_t *threads=(pthread_t *)malloc(sizeof(pthread_t)*nThreads);
	struct levenArgs **args=(struct levenArgs **)malloc(sizeof(struct levenArgs*)*nThreads);
	if(*debug)printf("Starting %d threads",nThreads);
	for(i=0;i<nThreads;i++){
		args[i]=(struct levenArgs *)malloc(sizeof(struct levenArgs));
		args[i]->answer=answer; args[i]->s1=s1;args[i]->s2=s2;args[i]->homoLimit=homoLimit;args[i]->prepend=prepend;args[i]->append=append;args[i]->debug=debug;args[i]->isAlign=isAlign;args[i]->align=align;args[i]->nStrings=nStrings;
		
		stepSize=(nStrings[0]-nextI)/(nThreads-i); //integer division
		args[i]->iRange[0]=nextI;
		args[i]->iRange[1]=nextI+stepSize-1;
		if(*debug){printf("irange[%d]: %d-%d \n",i,args[i]->iRange[0],args[i]->iRange[1]);}
		nextI+=stepSize;
		//not doing anything smart with j now
		args[i]->jRange[0]=0;
		args[i]->jRange[1]=nStrings[1]-1;
		if(pthread_create(&threads[i],NULL,levenAllPar,args[i])){errorMessage("Couldn't create thread",9);}
	}
	for(i=0;i<nThreads;i++){
		if(pthread_join(threads[i],NULL)){errorMessage("Couldn't join thread",10);}
	}

	for(i=0;i<nThreads;i++)free(args[i]);
	free(args);
	free(threads);



}

