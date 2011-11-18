##File Name: levenshtein.R
##Creation Date: Jan 01, 2009
##Last Modified: Wed 10 Jun 2009 05:03:31 PM EDT
##Created By: scott
##Summary: Functions to calculate levenshtein distance between strings (uses C code from c/leven.c)

levenSO<-'~/scripts/R/c/leven.so'
loader<-try(dyn.load(levenSO),TRUE)


if(any(grep("Error",loader))) stop(simpleError('Error loading levenshtein c functions'))


multiMismatch<-function(patterns, subject,drop=TRUE,...){
	matches<-do.call(rbind,mapply(function(x,y){matches<-bestMismatch(x,subject,TRUE,...);ifelse(is.null(dim(matches)),list(c(matches,y)),list(cbind(matches,y)))[[1]]},patterns,1:length(patterns),SIMPLIFY=FALSE))
	colnames(matches)<-c('mismatch','pos','number')
	rownames(matches)<-NULL
	if(drop) return(matches[matches[,1]==min(matches[,1]),,drop=FALSE][1,])
	else return(matches)
}

mismatchTable<-function(patterns,subjects,...){
	out<-do.call(rbind,lapply(patterns,function(x)sapply(subjects,function(y)bestMismatch(x,y,FALSE,...),USE.NAMES=FALSE)))
	return(out)
}
bestMismatch <- function(pattern, subject, pos=findAll, weights=rep(1,nchar(pattern)),findAll=FALSE,dummyChar='|'){
	if(nchar(pattern)>nchar(subject))stop(simpleError('Pattern longer than subject'))
	if(nchar(pattern)!=length(weights))stop(simpleError('Weights and pattern length do not match'))
	ans<-.C('bestMismatch',as.integer(c(-1,-1)),as.character(pattern),as.character(subject),as.integer(weights))[[1]]
	if(any(ans<0))stop(simpleError('Something did not work in mismatch'))
	stillLooking<-TRUE;nextMatch<-out<-ans
	if(findAll){
		while(stillLooking&findAll){
			substring(subject,nextMatch[2],nextMatch[2]+nchar(pattern)-1)<-paste(rep(dummyChar,nchar(pattern)),collapse='')
			nextMatch=bestMismatch(pattern,subject,TRUE,weights)
			if(nextMatch[1]==ans[1])out<-rbind(out,nextMatch)
			else stillLooking<-FALSE
		}
		return(out)
	}else{
		if(pos)return(ans)
		else return(ans[[1]][1])
	}
}


alignStringsToString<-function(strings,ref,homoLimit=0,prepend=NULL,append=NULL,substring1=FALSE,substring2=FALSE){
	if(substring1){append<-c(append,1);prepend<-c(prepend,1)}
	if(substring2){append<-c(append,2);prepend<-c(prepend,2)}
	append<-1:2 %in% append
	prepend<-1:2 %in% prepend
	aligns<-lapply(strings,function(x)levenAll(x,ref,homoLimit=homoLimit,prepend=prepend,append=append,align=TRUE)[[2]])
	refs<-c(ref,sapply(aligns,'[[',2))
	aligns<-c(ref,sapply(aligns,'[[',1))
	out<-rep('',length(strings)+1)
	#rather inefficient
	while(any(nchar(refs)>0)){
		selector<-substring(refs,1,1)=='-'
		if(any(selector)){
			out[!selector]<-sprintf('%s%s',out[!selector],'-')
			refs[selector]<-substring(refs[selector],2)
			out[selector]<-sprintf('%s%s',out[selector],substring(aligns[selector],1,1))
			aligns[selector]<-substring(aligns[selector],2)
		}else{
			out<-sprintf('%s%s',out,substring(aligns,1,1))
			refs<-substring(refs,2)
			aligns<-substring(aligns,2)
		}
	}
	#refGaps<-lapply(gregexpr('-',refs,fixed=TRUE),function(x)x[x!=-1])
	#uniqueGaps<-sort(unique(unlist(refGaps)))
	#finalOut<-vector('list',length(strings)+1)
	#for(i in uniqueGaps){
		#selector<-sapply(refGaps,function(x)i %in% refGaps)#could optimize here
	#}
	return(list('ref'=out[1],'align'=out[-1]))
}

levenAll <- function(string1, string2,distance=FALSE,homoLimit=0,debug=FALSE,prepend=NULL,append=NULL,align=FALSE) {
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
	
	#don't forget to deal with runs
	if(align){
		isAlign<-1
		outStrings<-rep(paste(rep('Z',2*(n1+n2)),collapse=''),2)
	}else{
		isAlign<-0
		outStrings<-c('','')
	}
	ans<-lapply(runs,function(x).C('levenAll',as.integer(1),as.character(string1),as.character(string2),as.integer(homoLimit),as.integer(x[1:2]),as.integer(x[3:4]),as.integer(debug),as.integer(isAlign),as.character(outStrings))[c(1,9)])
	#FIX UP DISTANCE CALCULATION
#	if(distance)ans[[1]]
#		ans2<-.C('levenAll',as.integer(1),as.character(string2),as.character(string1),as.integer(subString),as.integer(homoLimit),as.integer(prepend),as.integer(append),as.integer(debug))
#		if(distance){if(ans2[[1]][1]/n2<ans[[1]][1]/n1){ans<-ans2;string1<-string2;n1<-n2}} #replacing answer and string for later calculation
#		else{if(ans2[[1]][1]<ans[[1]][1])ans<-ans2}
#	}
	selector<-which.min(sapply(ans,'[[',1))
	if(!align)output<-ans[[selector]][[1]]
	else output<-ans[[selector]]
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
		if(vocal[1]>0&i%%vocal[1]==0)message("Working on string ",i)
		for(j in ifelse(multiToMulti,i,1):nStrings2){
			if(!is.na(vocal[2])&j%%vocal[2]==0)message("Working on string ",i,' - ',j)
			dist.mat[i,j]<-levenAll(strings1[i],strings2[ifelse(oneToOne,i,j)],distance=distance,homoLimit=homoLimit,debug=debug,append=append,prepend=prepend)
			if(multiToMulti)dist.mat[j,i]<-dist.mat[i,j]
		}
	}
	return(dist.mat[,])
}


