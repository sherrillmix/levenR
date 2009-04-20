//#include <string>
//#include <algorithm>
//#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
void leven(int *answer, char **s1, char **s2) {
	const unsigned int cost_del = 1;
	const unsigned int cost_ins = 1;
	const unsigned int cost_sub = 1;
	unsigned int n1 = strlen(s1[0]);// s1.length();
	unsigned int n2 = strlen(s2[0]);
	//printf("%d %d\n",n1,n2);
	unsigned int* p = malloc(sizeof(unsigned int)*(n2+1));
	unsigned int* q = malloc(sizeof(unsigned int)*(n2+1));//alloc unsigned int[n2+1];
	unsigned int* r;

	p[0] = 0;
	for( unsigned int j = 1; j <= n2; ++j ) p[j] = p[j-1] + cost_ins;

	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	for( unsigned int i = 1; i <= n1; ++i ){
			q[0] = p[0] + cost_del;
			for( unsigned int j = 1; j <= n2; ++j ){
					unsigned int d_del = p[j] + cost_del;
					unsigned int d_ins = q[j-1] + cost_ins;
					//printf("%d %d",i,j);
					unsigned int d_sub = p[j-1] + ( s1[0][i-1] == s2[0][j-1] ? 0 : cost_sub );
					//q[j] = min( min( d_del, d_ins ), d_sub );
					q[j]=d_ins;
					if (d_del < q[j])q[j]=d_del;
					if (d_sub < q[j])q[j]=d_sub;	
			}
			r = p;
			p = q;
			q = r;
	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	}
	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	*answer = p[n2];
	free(p);
	free(q);
}

void levenSub(int *answer, char **s1, char **s2) {
	const unsigned int cost_del = 1;
	const unsigned int cost_ins = 1;
	const unsigned int cost_sub = 1;
	unsigned int n1 = strlen(s1[0]);// s1.length();
	unsigned int n2 = strlen(s2[0]);
	//printf("%d %d\n",n1,n2);
	unsigned int* p = malloc(sizeof(unsigned int)*(n2+1));
	unsigned int* q = malloc(sizeof(unsigned int)*(n2+1));//alloc unsigned int[n2+1];
	unsigned int* r;
	unsigned int min;

	p[0] = 0;
	for( unsigned int j = 1; j <= n2; ++j ) p[j] = 0;

	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	for( unsigned int i = 1; i <= n1; ++i ){
			q[0] = p[0] + cost_del;
			for( unsigned int j = 1; j <= n2; ++j ){
					unsigned int d_del = p[j] + cost_del;
					unsigned int d_ins = q[j-1] + cost_ins;
					//printf("%d %d",i,j);
					unsigned int d_sub = p[j-1] + ( s1[0][i-1] == s2[0][j-1] ? 0 : cost_sub );
					//q[j] = min( min( d_del, d_ins ), d_sub );
					q[j]=d_ins;
					if (d_del < q[j])q[j]=d_del;
					if (d_sub < q[j])q[j]=d_sub;	
			}
			r = p;
			p = q;
			q = r;
	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	}
	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	min=p[0];
	for( unsigned int j = 1; j <= n2; ++j ){
		if(p[j]<min)min=p[j];
	}
	*answer = min;
	free(p);
	free(q);
}

void levenHomo(int *answer, char **s1, char **s2, int *homoLimit) {
	unsigned int cost_del = 1;
	unsigned int cost_ins = 1;
	const unsigned int cost_sub = 1;
	unsigned int n1 = strlen(s1[0]);// s1.length();
	unsigned int n2 = strlen(s2[0]);
	//printf("%d %d\n",n1,n2);
	unsigned int* p = malloc(sizeof(unsigned int)*(n2+1));
	unsigned int* q = malloc(sizeof(unsigned int)*(n2+1));//alloc unsigned int[n2+1];
	unsigned int* r;

	p[0] = 0;
	for( unsigned int j = 1; j <= n2; ++j ) p[j] = p[j-1] + cost_ins;

	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	for( unsigned int i = 1; i <= n1; ++i ){
			cost_del=0;
			if(i>*homoLimit){
				for(unsigned int step = 1; step <= *homoLimit+1; step++){
					if(s1[0][i]!=s1[0][i-step]){
						cost_del = 1;
						break;	
					}
				}
			}else cost_del=1;
			q[0] = p[0] + cost_del;
			for( unsigned int j = 1; j <= n2; ++j ){
				cost_ins=0;
				if(j>*homoLimit){
					for(unsigned int step = 1; step <= *homoLimit; step++){
						if(s2[0][j]!=s2[0][j-step]){
							cost_ins = 1;
							break;	
						}
					}
				}else cost_ins=1;
				unsigned int d_del = p[j] + cost_del;
				unsigned int d_ins = q[j-1] + cost_ins;
				//printf("%d %d",i,j);
				unsigned int d_sub = p[j-1] + ( s1[0][i-1] == s2[0][j-1] ? 0 : cost_sub );
				//q[j] = min( min( d_del, d_ins ), d_sub );
				q[j]=d_ins;
				if (d_del < q[j])q[j]=d_del;
				if (d_sub < q[j])q[j]=d_sub;	
			}
			r = p;
			p = q;
			q = r;
	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	}
	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	*answer = p[n2];
	free(p);
	free(q);
}

void levenAll(int *answer, char **s1, char **s2,int *subString, int *homoLimit) {
	unsigned int cost_del = 1;
	unsigned int cost_ins = 1;
	const unsigned int cost_sub = 1;
	unsigned int n1 = strlen(s1[0]);// s1.length();
	unsigned int n2 = strlen(s2[0]);
	//printf("%d %d\n",n1,n2);
	unsigned int* p = malloc(sizeof(unsigned int)*(n2+1));
	unsigned int* q = malloc(sizeof(unsigned int)*(n2+1));//alloc unsigned int[n2+1];
	unsigned int* r;
	unsigned int min;

	p[0] = 0;
	if(*subString) for( unsigned int j = 1; j <= n2; ++j ) p[j] = 0;
	else for( unsigned int j = 1; j <= n2; ++j ) p[j] = p[j-1] + cost_ins;
	

	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	for( unsigned int i = 1; i <= n1; ++i ){
		cost_del=0;
		if(*homoLimit&&i>*homoLimit){
			for(unsigned int step = 1; step <= *homoLimit+1; step++){
				if(s1[0][i]!=s1[0][i-step]){
					cost_del = 1;
					break;	
				}
			}
		}else cost_del=1;
		q[0] = p[0] + cost_del;
		for( unsigned int j = 1; j <= n2; ++j ){
			cost_ins=0;
			if(*homoLimit&&j>*homoLimit){
				for(unsigned int step = 1; step <= *homoLimit; step++){
					if(s2[0][j]!=s2[0][j-step]){
						cost_ins = 1;
						break;	
					}
				}
			}else cost_ins=1;
			unsigned int d_del = p[j] + cost_del;
			unsigned int d_ins = q[j-1] + cost_ins;
			//printf("%d %d",i,j);
			unsigned int d_sub = p[j-1] + ( s1[0][i-1] == s2[0][j-1] ? 0 : cost_sub );
			//q[j] = min( min( d_del, d_ins ), d_sub );
			q[j]=d_ins;
			if (d_del < q[j])q[j]=d_del;
			if (d_sub < q[j])q[j]=d_sub;	
		}
		r = p;
		p = q;
		q = r;
//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	}
	//for (int i = 1; i <=n2; ++i) printf(" %d ",p[i]);
	if(*subString){
		min=p[0];
		for( unsigned int j = 1; j <= n2; ++j ){
			if(p[j]<min)min=p[j];
		}
		*answer = min;
	}else *answer = p[n2];
	//printf(" %d ",*answer);
	free(p);
	free(q);
	
}



void levenMany(int *answer, char **strings, int *subString, int *homoLimit, int *numberStrings) {
//loop through strings comparing each to each don't overlap
//pass in answer num*num but fill 2 on each
//levenAll(int *answer, char **s1, char **s2,int *subString, int *homoLimit) 
	unsigned int dist;
	int* thisAnswer = malloc(sizeof(unsigned int));
	//char **s1;
	//char **s2;
	for(unsigned int i=0;i<*numberStrings;i++){
		for(unsigned int j=i;j<*numberStrings;j++){	
			if(i==j) answer[i+j * *numberStrings]=0;
			else{
				//s1=&strings[i];
				//s2=&strings[j];
				levenAll(thisAnswer,&strings[i],&strings[j],homoLimit,numberStrings);
				answer[i+j * *numberStrings]=*thisAnswer;
				answer[j+i * *numberStrings]=*thisAnswer;
			}
		}
	}
}

void levenOneToMany(int *answer, char **string, char **strings, int *subString, int *homoLimit, int *numberStrings) {
//loop through strings comparing each to each don't overlap
//pass in answer num*num but fill 2 on each
//levenAll(int *answer, char **s1, char **s2,int *subString, int *homoLimit) {
	for(unsigned int i=0;i<*numberStrings;i++){
	
	}
}

#define SWAP_CHAR( x, y ) {char c; c = x; x = y; y = c;}
void reverse(char t[])
{
  int i,j;
  for(i = 0, j = strlen(t)-1; i < j; i++, j--)
    SWAP_CHAR(t[i], t[j]);
}


#define INS_ID 1
#define DEL_ID 2
#define MATCH_ID 4

void levenAlign(int *answer, char **align, char **s1, char **s2,int *subString, int *homoLimit, int *gapLimit,int *debug) {
	unsigned int cost_del = 1;
	unsigned int cost_ins = 1;
	const unsigned int cost_sub = 1;
	unsigned int n1 = strlen(s1[0]);// s1.length();
	unsigned int n2 = strlen(s2[0]);
	//printf("% d \n",*gapLimit);
	unsigned int* p = malloc(sizeof(unsigned int)*(n2+1));
	unsigned int* q = malloc(sizeof(unsigned int)*(n2+1));//alloc unsigned int[n2+1];
	unsigned int* r;
	unsigned int min;
	unsigned int gaps =0;

	int **traceback = (int **)malloc((n1+1) * sizeof(int *));
	traceback[0] = (int *)malloc((n1+1) * (n2+1) * sizeof(int));
	for(unsigned int i = 1; i <= n1; i++)
		traceback[i] = traceback[0] + i * (n2+1);
	traceback[0][0]=0;

	p[0] = 0;
	for( unsigned int j = 1; j <= n2; ++j ){
		if(*subString) p[j] = 0;
		else p[j] = p[j-1] + cost_ins;
		traceback[0][j] = INS_ID;	
	}
	

	//for (int i = 1; i <=n2; ++i) printf(" %d ",q[i]);
	for( unsigned int i = 1; i <= n1; ++i ){
		cost_del=1;
		if(*homoLimit&&i>*homoLimit){
			cost_del=0;
			for(unsigned int step = 1; step <= *homoLimit+1; step++){
				if(s1[0][i]!=s1[0][i-step]){
					cost_del = 1;
					break;	
				}
			}
		}
		q[0] = p[0] + cost_del;
		traceback[i][0] = DEL_ID;	
		for( unsigned int j = 1; j <= n2; ++j ){
			cost_ins=1;
			if(*homoLimit&&j>*homoLimit){
				cost_ins=0;
				for(unsigned int step = 1; step <= *homoLimit; step++){
					if(s2[0][j]!=s2[0][j-step]){
						cost_ins = 1;
						break;	
					}
				}
			}
			if(*gapLimit && gaps>*gapLimit) cost_ins=0;	

			//printf(" %d-%d",gaps,cost_ins);
			unsigned int d_del = p[j] + cost_del;
			unsigned int d_ins = q[j-1] + cost_ins;
			//printf("%d %d",i,j);
			unsigned int d_sub = p[j-1] + ( s1[0][i-1] == s2[0][j-1] ? 0 : cost_sub );
			//q[j] = min( min( d_del, d_ins ), d_sub );
			q[j]=d_ins;
			traceback[i][j]=0;
			if (d_del < q[j]){
				//insert gap sequence 2 and move up trackback
				q[j]=d_del;
				traceback[i][j]+=DEL_ID;
			}
			if (d_sub < q[j]){
				//match and move diagonal trackback
				q[j]=d_sub;
				traceback[i][j]+=MATCH_ID;
			}
			if (d_ins == q[j]){
				//insert gap sequence 1 and move left trackback
				traceback[i][j]+=INS_ID;
				gaps++;
			}else gaps=0;
		}
		r = p;
		p = q;
		q = r;
		if(*debug){for (int i = 1; i <=n2; ++i) printf(" %d ",p[i]);printf("\n");}
	}
	if(*debug){for (int i = 0; i <=n1; i++){for (int j = 0; j<=n2;j++){ printf(" %d ",traceback[i][j]);} printf("\n");}}
	//for (int i = 1; i <=n2; ++i) printf(" %d ",p[i]);
	int cursorx=n1;
	int cursory=n2;
	if(*subString){
		min=p[0];
		for( unsigned int j = 1; j <= n2; ++j ){
			if(p[j]<min){
				min=p[j];
				cursory=j-1;
			}
		}
		*answer = min;
	}else *answer = p[n2];
	int counter=0;
	char convert[2];
	convert[1]='\0';
	char convert2[2];
	convert2[1]='\0';
	strcpy(align[0],"");
	strcpy(align[1],"");
	while(cursorx > 0 || (cursory >0 && !*subString)){
		if(*debug)printf(" %d %d \n",cursorx,cursory);
		//counter++;
		//if(counter>20)break;
		if(traceback[cursorx][cursory] & MATCH_ID){
			convert[0]=s1[0][cursorx-1];
			convert2[0]=s2[0][cursory-1];
			cursorx--;
			cursory--;
		}else if(traceback[cursorx][cursory] & DEL_ID){
			convert[0]=s1[0][cursorx-1];
			convert2[0]='-';	
			cursorx--;
		}else if(traceback[cursorx][cursory] & INS_ID){
			convert[0]='-';	
			convert2[0]=s2[0][cursory-1];
			cursory--;
		}else{
			//We've got a problem
			printf("Something went wrong");
			break;
		}
		strcat(align[0],convert);	
		strcat(align[1],convert2);	
		if(*debug){printf("|%s|%s|%s|%s|\n",align[0],align[1],convert,convert2);}
	}
	
	reverse(align[0]);
	reverse(align[1]);
	
	//printf(" %d ",*answer);
	free(p);
	free(q);
	
	free(traceback[0]);
	free(traceback);
}

