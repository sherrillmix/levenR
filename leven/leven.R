#' Simple Levenshtein alignments
#'
#' \tabular{ll}{
#' Package: \tab leven\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1\cr
#' Date: \tab 2011-11-16\cr
#' License: \tab CC BY-NC 3.0\cr
#' LazyLoad: \tab yes\cr
#' }
#' 
#' @name leven-package
#' @aliases levenstein
#' @docType package
#' @title Simple Levenshtein alignments
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#' @references
#' \url{http://en.wikipedia.org/wiki/Levenshtein_distance}
#' @keywords alignment, string, distance, levenshtein
#' @seealso \code{Biostrings}
#' @examples
#' library(leven)
#' leven('ACAC','ACACA')
NULL

#' Reverse a string
#'
#' Reverse a string e.g. ABCDE becomes EDCBA
#' @param strings Vector of strings to be revered
#' @return Vector of reversed strings
reverseString<-function(strings){
	output<-sapply(strsplit(strings,''),function(x)paste(rev(x),collapse=''))	
	return(output)
}

#' Compliment DNA
#'
#' Compliment a string of DNA e.g. A<->T,G<->C. Doesn't attempt to deal with ambigous bases.
#' @param dnas Vector of DNA sequences
#' @return Vector of complimented DNA strings
complimentDna<-function(dnas,brackets=TRUE,ambigs=TRUE){
	finds<-'TGAC'
	replaces<-'ACTG'
	return(chartr(finds,replaces,dnas))
}

#' Reverse compliment dna
#'
#' Reverse compliment a string of DNA
#' @param dnas Vector of sequences
#' @return Vector of reverse complimented DNA strings
revComp<-function(dnas){
	return(complimentDna(reverseString(dnas)))
}

#' Find fewest mismatches for multiple patterns
#'
#' Find fewest mismatches (and position if desired) between pattern and string without allowing gaps. Could be programmed better but quick enough for now
#' @param patterns A vector of string patterns to search for in subject
#' @param subject A single subject string to search in
#' @param ... Additional arguments for bestMismatch
#' @param drop If TRUE return only return the first best match
#' @return Three column matrix with columns minMismatch, position, pattern number 
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#" @export
multiMismatch<-function(patterns, subject,drop=TRUE,...){
	matches<-do.call(rbind,mapply(function(x,y){matches<-bestMismatch(x,subject,TRUE,...);ifelse(is.null(dim(matches)),list(c(matches,y)),list(cbind(matches,y)))[[1]]},patterns,1:length(patterns),SIMPLIFY=FALSE))
	colnames(matches)<-c('mismatch','pos','number')
	rownames(matches)<-NULL
	if(drop) return(matches[matches[,1]==min(matches[,1]),,drop=FALSE][1,])
	else return(matches)
}

#' Find fewest mismatches
#'
#' Find fewest mismatches (and position if desired) between pattern and string without allowing gaps. Could be programmed better but quick enough for now.
#'
#' @param pattern String pattern to search for in subject
#' @param subject Subject string to search in
#' @param pos If TRUE return position in addition to distance
#' @param weights Integer weights for each position of pattern
#' @param findAll If TRUE return all best matches
#' @param dummyChar A character that does not appear in either string (only necessary if findAll==TRUE [and could probably program better])
#' @return Single value giving the least mismatch if pos=FALSE, c(minMismatches,position) if pos=TRUE, if findAll==TRUE matrix with columns minMismatches, position and rows for each occurrence
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#" @export
bestMismatch <- function(pattern, subject, pos=findAll, weights=rep(1,nchar(pattern)),findAll=FALSE,dummyChar='|'){
	if(nchar(pattern)>nchar(subject))stop(simpleError('Pattern longer than subject'))
	if(nchar(pattern)!=length(weights))stop(simpleError('Weights and pattern length do not match'))
	ans<-.C('bestMismatch',as.integer(c(-1,-1)),as.character(pattern),as.character(subject),as.integer(weights),PACKAGE='leven')[[1]]
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

#' Combine gapped alignments
#'
#' Combine several pairs of reference, query alignments
#'
#' @param refs The gapped reference sequences (one for each query)
#' @param aligns The gapped query sequences
#' @return A vector of the aligned strings with gaps added
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#" @export
combineAligns<-function(refs,aligns){
	out<-rep('',length(aligns))

	#deal with starts by stepping backward from ref
	lastGap<-regexpr('[^-*]',refs)-1
	maxGap<-max(lastGap)
	dummy<-paste(rep('-',maxGap),collapse='')
	out<-sprintf('%s%s%s',out,substring(dummy,1,maxGap-lastGap),substring(aligns,1,lastGap))
	aligns<-substring(aligns,lastGap+1)
	refs<-substring(refs,lastGap+1)

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
	return(out)
}


#' Align many strings to a reference
#'
#' Align strings to a reference string based on Levenshtein distance with ends-free or homopolymer alignments possible
#'
#' @param strings The strings to be aligned
#' @param ref The reference to align to
#' @param homoLimit deletions or insertions in homopolymers > homoLimit cost 0
#' @param substring1 ends-free for strings1
#' @param substring2 ends-free for strings2
#' @param prepend 1 or 2 for ends-free on front of string 1 or 2
#' @param append 1 or 2 for ends-free on back of string 1 or 2
#' @param Find and return the best alignment of identity or reverse compliment of queries
#' @return A list with first entry for the gapped reference and second entry all the aligned strings
#' @references \url{http://en.wikipedia.org/wiki/Levenshtein_distance}
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#" @export
levenAlign<-function(strings,ref,homoLimit=0,prepend=NULL,append=NULL,substring1=FALSE,substring2=FALSE,trimOuterGaps=FALSE,revComp=FALSE){
	if(length(ref)>1)stop(simpleError('Only a single reference supported'))
	if(substring1){append<-c(append,1);prepend<-c(prepend,1)}
	if(substring2){append<-c(append,2);prepend<-c(prepend,2)}
	append<-1:2 %in% append
	prepend<-1:2 %in% prepend
	aligns<-lapply(strings,function(x)levenAll(x,ref,homoLimit=homoLimit,prepend=prepend,append=append,align=TRUE,revComp=revComp)[[2]])
	refs<-c(ref,sapply(aligns,'[[',2))
	aligns<-c(ref,sapply(aligns,'[[',1))
	out<-combineAligns(refs,aligns)
	if(trimOuterGaps){
		refRange<-range(gregexpr('[^-]',out[1])[[1]])
		nonGaps<-gregexpr('[^-]',out[-1])
		alignRange<-c(min(sapply(nonGaps,min)),max(sapply(nonGaps,max)))
		trimRange<-c(max(refRange[1],alignRange[1]),min(refRange[2],alignRange[2]))
		out<-substring(out,trimRange[1],trimRange[2])
	}

	#refGaps<-lapply(gregexpr('-',refs,fixed=TRUE),function(x)x[x!=-1])
	#uniqueGaps<-sort(unique(unlist(refGaps)))
	#finalOut<-vector('list',length(strings)+1)
	#for(i in uniqueGaps){
		#selector<-sapply(refGaps,function(x)i %in% refGaps)#could optimize here
	#}
	out<-list('ref'=out[1],'align'=out[-1])
	if(trimOuterGaps)attr(out,'start')<-trimRange[1]
	return(out)
}

#' Wrapper for .C Levenshtein function
#'
#' Align and calculate distance betweeen two strings based on Levenshtein distance with ends-free or homopolymer alignments possible (probably not called directly)
#'
#' @param string1 A single string 
#' @param string2 Another single string
#' @param homoLimit deletions or insertions in homopolymers > homoLimit cost 0
#' @param debug Output various debugging information
#' @param prepend 1 or 2 for ends-free on front of string 1 or 2
#' @param append 1 or 2 for ends-free on back of string 1 or 2
#' @param align Should an alignment be calculated
#' @param revComp Should the reverse compliment of string1 also be tested and the best of identity or reverse compliment returned?
#' @param nThreads If threads >1, run nThreads threads in parallel 
#' @return Either a distance or if align==TRUE a list with distance and alignment
#' @references \url{http://en.wikipedia.org/wiki/Levenshtein_distance}
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#" @export
levenAll <- function(string1, string2,homoLimit=0,debug=FALSE,prepend=NULL,append=NULL,revComp=FALSE,align=FALSE,nThreads=1) {
	if(align&nThreads>1)stop(simpleError('Parallel and align not compatible yet'))
	if(is.null(prepend))prepend<-c(FALSE,FALSE)
	if(is.null(append))append<-c(FALSE,FALSE)
	if(length(append)!=2&length(prepend!=2))stop(simpleError('Specify prepend and append as a vector of two logicals'))
	if(any(!is.logical(c(prepend,append))))stop(simpleError('Specify prepend and append as a vector of two logicals'))
	if(any(prepend & rev(append)))runs<-list(c(prepend[1],FALSE,append[1],FALSE),c(FALSE,prepend[2],FALSE,append[2]))
	else runs<-list(c(prepend,append))
	if(revComp)runs<-c(lapply(runs,function(x)c(x,FALSE)),lapply(runs,function(x)c(x,TRUE)))
	else runs<-lapply(runs,function(x)c(x,FALSE))

	if(any(is.na(string1))|any(is.na(string2)))stop(simpleError('NULL or NA string in levenAll'))

	if(nThreads==1){
		#make sure we're only dealing with 1 string each
		string1<-string1[1];string2<-string2[1]
		out<--99
		#don't forget to deal with runs
		if(align){
			n1<-nchar(string1)
			n2<-nchar(string2)
			isAlign<-1
			outStrings<-rep(paste(rep('Z',2*(n1+n2)),collapse=''),2)
		}else{
			isAlign<-0
			outStrings<-c('','')
		}
		ans<-lapply(runs,function(x).C('levenAll',as.integer(out),as.character(if(x[5]) revComp(string1) else string1),as.character(string2),as.integer(homoLimit),as.integer(x[1:2]),as.integer(x[3:4]),as.integer(debug),as.integer(isAlign),as.character(outStrings),PACKAGE='leven')[c(1,9)])
		selector<-which.min(sapply(ans,'[[',1))
		if(!align)output<-ans[[selector]][[1]]
		else output<-ans[[selector]]
	}else{
		out<-rep(-99,length(string1)*length(string2))
		nStrings<-c(length(string1),length(string2))
		ans<-do.call(rbind,lapply(runs,function(x).C('parallelLeven',as.integer(out),as.character(string1),as.character(string2),as.integer(nStrings),as.integer(homoLimit),as.integer(x[1:2]),as.integer(x[3:4]),as.integer(debug),as.integer(nThreads),PACKAGE='leven')[[1]]))
		ans<-apply(ans,2,min)
		output<-matrix(ans,nrow=nStrings[1],ncol=nStrings[2])
	}
		
	return(output)
}

#' Find Levenshtein distance between many strings
#' @param strings1 vector of strings to compare to strings2 (or as strings1 x strings1 if strings2==NULL)
#' @param strings2 vector of strings to compare to strings1 or NULL (either strings1 or strings2 should be 1 or should have equal lengths)
#' @param oneToOne compare strings1[1] to strings2[1], strings1[2] to strings2[2],... (both vectors must be same length)
#' @param homoLimit deletions or insertions in homopolymers > homoLimit cost 0
#' @param vocal Integer to display a status message every vocal iterations
#' @param debug If TRUE display various debugging information
#' @param prepend 1 or 2 for ends-free on front of string 1 or 2
#' @param append 1 or 2 for ends-free on back of string 1 or 2
#' @param substring1 ends-free for strings1
#' @param substring2 ends-free for strings2
#' @param nThreads Run using threads?
#' @return Levenshtein distance matrix
#' @references \url{http://en.wikipedia.org/wiki/Levenshtein_distance}
#' @author Scott Sherrill-Mix \email{R@@sherrillmix.com}
#" @export
leven<-function(strings1,strings2=NULL,oneToOne=FALSE,homoLimit=0,vocal=0,debug=FALSE,prepend=NULL,append=NULL,substring1=FALSE,substring2=FALSE,nThreads=1){
	if(nThreads>1&oneToOne)warnings('oneToOne not currently compatible with parallel processing')
	if(any(!prepend %in% 1:2)|any(!append %in% 1:2))stop(simpleError('Specify prepend and append with 1 or 2'))
	if(is.null(strings2)){
		strings2<-strings1
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

	#doesn't seem to be any speed benefit to using lapply or doing the looping in C (strange)
	if(oneToOne){
		distMat<-rep(NA,nStrings1)
		for(i in 1:nStrings1){
			if(vocal[1]>0&i%%vocal[1]==0)message("Working on string ",i)
			distMat[i]<-levenAll(strings1[i],strings2[i],homoLimit=homoLimit,debug=debug,append=append,prepend=prepend)
		}
	}else{
		distMat<-matrix(NA,nrow=nStrings1,ncol=nStrings2)
		if(nThreads==1){
			for(i in 1:nStrings1){
				for(j in 1:nStrings2){
					distMat[i,j]<-levenAll(strings1[i],strings2[j],homoLimit=homoLimit,debug=debug,append=append,prepend=prepend)
				}
			}
		}else{
			distMat<-levenAll(strings1,strings2,homoLimit=homoLimit,debug=debug,append=append,prepend=prepend,nThreads=nThreads)
		}
	}
	return(distMat)
}


.onLoad<-function(lib, pkg){
	library.dynam("leven", pkg, lib)
}


.onUnload<-function(path){
	library.dynam.unload("leven", path)
}


