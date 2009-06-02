levenSO<-'~/scripts/R/c/leven.so'
leven2SO<-'~/scripts/R/c/leven2.so'
loader<-try(dyn.load(levenSO),TRUE)
loader2<-try(dyn.load(leven2SO),TRUE)
if (any(grep("Error",loader))){
	message('Using R leven function')
	leven <- function(string1, string2, case=TRUE, map=NULL) {
		
		########
		#
		# levenshtein algorithm in R
		#
		# Author  : Hans-Joerg Bibiko
		# Date    : 29/06/2006
		#
		# Contact : bibiko@eva.mpg.de
		#
		########
		#
		# string1, string2 := strings to compare
		#
		# case = TRUE := case sensitivity; case = FALSE := case insensitivity
		#
		# map := character vector of c(regexp1, replacement1, regexp2, replacement2, ...)
		#
		#   example:
		#      map <- c("[aeiou]","V","[^aeiou]","C") := replaces all vowels with V and all others with C
		#
		#      levenshtein("Bank","Bond", map=map)   =>  0
		#
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

}else{
	leven <- function(string1, string2,distance=FALSE,subString=FALSE) {
		if(subString) ans<-.C('levenSub',as.integer(1),as.character(string1),as.character(string2))
		else ans<-.C('leven',as.integer(1),as.character(string1),as.character(string2))	
		if (distance) return(ans[[1]][1]/max(nchar(string1),nchar(string2)))
		return(ans[[1]][1])
	}
}
levenAll <- function(string1, string2,distance=FALSE,subString=FALSE,homoLimit=0,subBoth=FALSE) {
	if(is.null(homoLimit))homoLimit<-0
	n1<-nchar(string1)
	n2<-nchar(string2)
	ans<-.C('levenAll',as.integer(1),as.character(string1),as.character(string2),as.integer(subString),as.integer(homoLimit))
	if(subString&subBoth&ans[[1]][1]>0){
		ans2<-.C('levenAll',as.integer(1),as.character(string2),as.character(string1),as.integer(subString),as.integer(homoLimit))
		if(distance){if(ans2[[1]][1]/n2<ans[[1]][1]/n1){ans<-ans2;string1<-string2;n1<-n2}} #replacing answer and string for later calculation
		else{if(ans2[[1]][1]<ans[[1]][1])ans<-ans2}
	}
	if (distance){
		if(subString) return(ans[[1]][1]/n1)
		else return(ans[[1]][1]/max(n1,n2))
	}
	return(ans[[1]][1])
}
#levenSub<-function(string1, string2,distance=FALSE) {
#	ans<-.C('levenSub',as.integer(1),as.character(string1),as.character(string2))	
#	if (distance) return(ans[[1]][1]/(nchar(string1)))
#	return(ans[[1]][1])
#}

#leven2 <- function(string1, string2,distance=FALSE) {
#	ans<-.C('leven2',as.integer(1),as.character(string1),as.character(string2))	
#	if (distance) return(ans[[1]][1]/max(nchar(string1),nchar(string2)))
#	return(ans[[1]][1])
#}

#levenHomo <- function(string1, string2,distance=FALSE,homoLimit) {
#	ans<-.C('levenHomo',as.integer(1),as.character(string1),as.character(string2),as.integer(homoLimit))	
#	if (distance) return(ans[[1]][1]/max(nchar(string1),nchar(string2)))
#	return(ans[[1]][1])
#}

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

levenOnetoMany<-function(string,strings,distance=FALSE,subString=FALSE,homoLimit=0,subBoth=FALSE,vocal=0){
	if(length(string)>1){warning("Only using first element in string");string<-string[1]}
	if(all(nchar(string)!=nchar(strings))&!subString)warning("String lengths not equal")
	num<-length(strings)
	dist.mat<-rep(NA,num)
	for(i in 1:num){
			if(vocal>0&i%%vocal==0)message('Working on string ',i)
			dist.mat[i]<-levenAll(string,strings[i],distance,subString,homoLimit,subBoth)	
	}
	return(dist.mat)
}


levenManytoMany<-function(strings1,strings2,distance=FALSE,subString=FALSE,homoLimit=0,vocal=0){
	if(length(strings1)!=length(strings2))stop(simpleError('Length of string vectors differ'))
	if(any(nchar(strings1)!=nchar(strings2))&!subString)warning("String lengths not equal")
	num<-length(strings1)
	dist.mat<-rep(NA,num)
	for(i in 1:num){
			if(vocal>0&&i%%vocal==0)message('Working on string ',i)
			dist.mat[i]<-levenAll(strings1[i],strings2[i],distance,subString,homoLimit)	
	}
	return(dist.mat)
}

levenManytoOne<-function(strings,string,distance=FALSE,subString=FALSE,homoLimit=0,vocal=0){
	if(length(string)>1){warning("Only using first element in string");string<-string[1]}
	if(all(nchar(string)!=nchar(strings))&!subString)warning("String lengths not equal")
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
