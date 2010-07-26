##File Name: levenshtein.R
##Creation Date: Jan 01, 2009
##Last Modified: Wed 10 Jun 2009 05:03:31 PM EDT
##Created By: scott
##Summary: Functions to calculate levenshtein distance between strings (uses C code from c/leven.c)

levenSO<-'~/scripts/R/c/leven.so'
loader<-try(dyn.load(levenSO),TRUE)


if(any(grep("Error",loader))) stop(simpleError('Error loading levenshtein c functions'))


multiMismatch<-function(patterns, subject,...){
	matches<-do.call(rbind,lapply(patterns,function(x)bestMismatch(x,subject,TRUE,...)))
	matches<-cbind(matches,1:nrow(matches))
	colnames(matches)<-c('mismatch','pos','number')
	return(matches[matches[,1]==min(matches[,1]),,drop=FALSE][1,])
}

bestMismatch <- function(pattern, subject, pos=FALSE, weights=rep(1,nchar(pattern))){
	if(nchar(pattern)>nchar(subject))stop(simpleError('Pattern longer than subject'))
	if(nchar(pattern)!=length(weights))stop(simpleError('Weights and pattern length do not match'))
	ans<-.C('bestMismatch',as.integer(c(-1,-1)),as.character(pattern),as.character(subject),as.integer(weights))[[1]]
	if(any(ans<0))stop(simpleError('Something did not work in mismatch'))
	if(pos)return(ans)
	else return(ans[[1]][1])
}

levenAll <- function(string1, string2,distance=FALSE,homoLimit=0,debug=FALSE,prepend=NULL,append=NULL) {
	if(is.null(prepend))prepend<-c(FALSE,FALSE)
	if(is.null(append))append<-c(FALSE,FALSE)
	if(length(append)!=2&length(prepend!=2))stop(simpleError('Specify prepend and append as a vector of two logicals'))
	if(any(!is.logical(c(prepend,append))))stop(simpleError('Specify prepend and append as a vector of two logicals'))
	if(any(prepend & rev(append))){
		runs<-list(c(prepend[1],FALSE,append[1],FALSE),c(FALSE,prepend[2],FALSE,append[2]))
	}else runs<-list(c(prepend,append))

	#make sure we're only dealing with 1 string each
	string1<-string1[1];string2<-string2[1]
	n1<-nchar(string1)
	n2<-nchar(string2)
	if(is.na(string1)|is.na(string2))stop(simpleError('NULL or NA string in levenAll'))
	####WORK HERE, CHECK APPEND PREPEND IF 1|1 2|2 |1 |2 1| 2| 1,2| |1,2 THEN RUN ONCE 1,2|1,2 1|2 2|1 1,2|1 1,2|2 1|1,2 2|1,2 RUN TWICE
	ans<-sapply(runs,function(x).C('levenAll',as.integer(1),as.character(string1),as.character(string2),as.integer(homoLimit),as.integer(x[1:2]),as.integer(x[3:4]),as.integer(debug))[[1]][1])
	#FIX UP DISTANCE CALCULATION
#	if(distance)ans[[1]]
#		ans2<-.C('levenAll',as.integer(1),as.character(string2),as.character(string1),as.integer(subString),as.integer(homoLimit),as.integer(prepend),as.integer(append),as.integer(debug))
#		if(distance){if(ans2[[1]][1]/n2<ans[[1]][1]/n1){ans<-ans2;string1<-string2;n1<-n2}} #replacing answer and string for later calculation
#		else{if(ans2[[1]][1]<ans[[1]][1])ans<-ans2}
#	}
	output<-min(ans)
	#if (distance){
	#	if(subString) return(ans[[1]][1]/n1)
	#	else return(ans[[1]][1]/max(n1,n2))
	#}
	return(output)
}

#levenStringsToStrings
#strings1: vector of strings to compare to strings2 (or as strings1 x strings1 if strings==NULL)
#strings2: vector of strings to compare to strings1 or NULL (either strings1 or strings2 should be 1 or should have equal lengths)
#oneToOne: compare strings1[1] to strings2[1], strings1[2] to strings2[2],... (both vectors must be same length)
#distance: if FALSE return edit distance, if TRUE return edit distance/min(nchar(seq1,seq2)
#homoLimit: deletions or insertions in homopolymers > homoLimit cost 0
#WORK HERE
#substring1
#substring2
#prepend
#append
levenStringsToStrings<-function(strings1,strings2=NULL,oneToOne=FALSE,distance=FALSE,homoLimit=0,vocal=0,debug=FALSE,prepend=NULL,append=NULL,substring1=FALSE,substring2=FALSE){
	if(any(!prepend %in% 1:2)|any(!append %in% 1:2))stop(simpleError('Specify prepend and append with 1 or 2'))
	multiToMulti<-FALSE
	if(is.null(strings2)){
		strings2<-strings1
		multiToMulti<-TRUE
	}
	nStrings1=length(strings1)
	nStrings2=length(strings2)
	if(oneToOne){
		if(nStrings1!=nStrings2)stop(simpleError('Length of strings1 and strings2 not equal for 1 to 1 dists'))
		nStrings2<-1
	}
	if(substring1){append<-c(append,1);prepend<-c(prepend,1)}
	if(substring2){append<-c(append,2);prepend<-c(prepend,2)}
	append<-1:2 %in% append
	prepend<-1:2 %in% prepend
	dist.mat<-matrix(NA,nrow=nStrings1,ncol=nStrings2)
	#doesn't seem to be any speed benefit to using lapply or doing the looping in C (strange)
	for(i in 1:ifelse(nStrings1>1,nStrings1,1)){
		if(vocal>0&i%%vocal==0)message("Working on string ",i)
		for(j in ifelse(multiToMulti,i,1):nStrings2){
			dist.mat[i,j]<-levenAll(strings1[i],strings2[ifelse(oneToOne,i,j)],distance=distance,homoLimit=homoLimit,debug=debug,append=append,prepend=prepend)
			if(multiToMulti)dist.mat[j,i]<-dist.mat[i,j]
		}
	}
	return(dist.mat[,])
}



#DEVELOPMENT FUNCTIONS###########################################
#source('levenshtein.R');levenAlign('ACTGACCC','TTTATATATAGGACTGATACCTATACCCAT',debug=TRUE,gapLimit=1,subString=TRUE)
#is acting funny at start without subString and with gapLimit
levenAlign <- function(string1, string2,distance=FALSE,subString=FALSE,homoLimit=0,gapLimit=0,debug=FALSE) {
	if(is.null(homoLimit))homoLimit<-0
	n1<-nchar(string1)
	n2<-nchar(string2)
	align<-rep(paste(string1,string2,sep=''),2)
	ans<-.C('levenAlign',as.integer(1),as.character(align),as.character(string1),as.character(string2),as.integer(subString),as.integer(homoLimit),as.integer(gapLimit),as.integer(debug))
	#print(ans);
	if (distance){
		if(subString) dist<-(ans[[1]][1]/n1)
		else dist<-(ans[[1]][1]/max(n1,n2))
	}else dist<-ans[[1]][1]
	return(list(ans[[2]],dist))
}


#LEGACY FUNCTIONS##################################################
levenDist<-function(strings,distance=FALSE,vocal=0,subString=FALSE,homoLimit=0,subBoth=FALSE){
	num<-length(strings)
	if(num<1)stop(simpleError('No strings input to levenDist'))
	dist.mat<-matrix(NA,ncol=num,nrow=num)
	for(i in 1:num){
		if(vocal>0&i%%vocal==0)message("Working on string ",i)
		for(j in i:num){
			dist.mat[i,j]<-levenAll(strings[i],strings[j],distance=distance,subString=subString,homoLimit=homoLimit,subBoth)			
			dist.mat[j,i]<-dist.mat[i,j]
		}
	}
	return(dist.mat)
}

levenOnetoMany<-function(string,strings,distance=FALSE,subString=FALSE,homoLimit=0,subBoth=FALSE,vocal=0,warn=TRUE){
	if(length(string)>1){warning("Only using first element in string");string<-string[1]}
	if(all(nchar(string)!=nchar(strings))&!subString&warn)warning("String lengths not equal")
	num<-length(strings)
	dist.mat<-rep(NA,num)
	for(i in 1:num){
			if(vocal>0&i%%vocal==0)message('Working on string ',i)
			dist.mat[i]<-levenAll(string,strings[i],distance,subString,homoLimit,subBoth)	
	}
	return(dist.mat)
}


levenManytoMany<-function(strings1,strings2,distance=FALSE,subString=FALSE,homoLimit=0,vocal=0,warn=TRUE){
	if(length(strings1)!=length(strings2))stop(simpleError('Length of string vectors differ'))
	if(any(nchar(strings1)!=nchar(strings2))&!subString&warn)warning("String lengths not equal")
	num<-length(strings1)
	dist.mat<-rep(NA,num)
	for(i in 1:num){
			if(vocal>0&&i%%vocal==0)message('Working on string ',i)
			dist.mat[i]<-levenAll(strings1[i],strings2[i],distance,subString,homoLimit)	
	}
	return(dist.mat)
}

levenManytoOne<-function(strings,string,distance=FALSE,subString=FALSE,homoLimit=0,vocal=0,warn=TRUE){
	if(length(string)>1){warning("Only using first element in string");string<-string[1]}
	if(all(nchar(string)!=nchar(strings))&!subString&warn)warning("String lengths not equal")
	num<-length(strings)
	dist.mat<-rep(NA,num)
	for(i in 1:num){
			if(vocal>0&&i%%vocal==0)message('Working on string ',i)
			dist.mat[i]<-levenAll(strings[i],string,distance,subString,homoLimit)	
	}
	return(dist.mat)
}

levenMany<- function(strings, distance=FALSE,substring=FALSE,homoLimit=0) {
	numStrings<-length(strings)
	ans<-.C('levenMany',as.integer(rep(1,length(strings)*numStrings)),as.character(strings),as.integer(substring),as.integer(homoLimit),as.integer(numStrings))
	answer<-matrix(ans[[1]],ncol=numStrings)
	if (distance) return(ans[[1]][1]/max(nchar(string1),nchar(string2)))
	return(ans[[1]][1])
}



levenR <- function(string1, string2, substring=FALSE) {
	s1 <- strsplit(paste(" ", string1, sep=""), NULL)[[1]]
	s2 <- strsplit(paste(" ", string2, sep=""), NULL)[[1]]
	
	l1 <- length(s1)
	l2 <- length(s2)
	
	d <- matrix(nrow = l1, ncol = l2)
 
	#if(substring) for(i in 1:l1) d[i,1] <- 0
	#else for(i in 1:l1) d[i,1] <- i-1
	for(i in 1:l1) d[i,1] <- i-1
	if(substring) for(i in 1:l2) d[1,i] <- 0
	else for(i in 1:l2) d[1,i] <- i-1


	for(i in 2:l1){
		#deletion cost
		if(i>4){
			if(s1[i]==s1[i-1]&s1[i]==s1[i-2]&s1[i]==s1[i-3]) delcost<-0
			else delcost<-1
		}
		else delcost<-1

		for(j in 2:l2){
			#substitution cost
			if(s1[i]==s2[j]) subcost<-0
			else subcost<-1

			#insertion cost
			if(j>4){
				if(s2[j]==s2[j-1]&s2[j]==s2[j-2]&s2[j]==s2[j-3]) inscost<-0
				else inscost<-1
			}else inscost<-1

			d[i,j] <- min( 
				(d[i-1,j]+delcost), 
				(d[i,j-1]+inscost), 
				(d[i-1,j-1]+subcost)
			)
		}
	}
	rownames(d)<-s1
	colnames(d)<-s2
	if (substring) answer<-min(d[l1,])
	else answer<-d[l1,l2]
	return(list(d,answer))
}


leven <- function(string1, string2, case=TRUE, map=NULL) {
	########
	# levenshtein algorithm in R
	# Author  : Hans-Joerg Bibiko
	# Date    : 29/06/2006
	# Contact : bibiko@eva.mpg.de
	########
	# string1, string2 := strings to compare
	# case = TRUE := case sensitivity; case = FALSE := case insensitivity
	# map := character vector of c(regexp1, replacement1, regexp2, replacement2, ...)
	#   example:
	#      map <- c("[aeiou]","V","[^aeiou]","C") := replaces all vowels with V and all others with C
	#      levenshtein("Bank","Bond", map=map)   =>  0
	########
	if(!is.null(map)) {
		m <- matrix(map, ncol=2, byrow=TRUE)
		s <- c(ifelse(case, string1, tolower(string1)), ifelse(case, string2, tolower(string2)))
		for(i in 1:dim(m)[1]) s <- gsub(m[i,1], m[i,2], s)
		string1 <- s[1]
		string2 <- s[2]
	}
 
	if(ifelse(case, string1, tolower(string1)) == ifelse(case, string2, tolower(string2))) return(0)
 
	s1 <- strsplit(paste(" ", ifelse(case, string1, tolower(string1)), sep=""), NULL)[[1]]
	s2 <- strsplit(paste(" ", ifelse(case, string2, tolower(string2)), sep=""), NULL)[[1]]
	
	l1 <- length(s1)
	l2 <- length(s2)
	
	d <- matrix(nrow = l1, ncol = l2)
 
	for(i in 1:l1) d[i,1] <- i-1
	for(i in 1:l2) d[1,i] <- i-1
	for(i in 2:l1) for(j in 2:l2) d[i,j] <- min((d[i-1,j]+1) , (d[i,j-1]+1) , (d[i-1,j-1]+ifelse(s1[i] == s2[j], 0, 1)))
	d[l1,l2]
}
