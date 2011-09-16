fileName:=$(shell basename $$(pwd))
version:=$(shell grep "Version: " DESCRIPTION|sed "s/Version: //")
all: ${fileName}_${version}.tar.gz

${fileName}_${version}.tar.gz: ${fileName} ${fileName}.R
	@echo Compiling package ${fileName}_${version}  
	R CMD check ${fileName}
	R CMD build ${fileName}

${fileName}:
	@echo Running Roxygen
	R --vanilla -e "package.skeleton('${fileName}',code_files='${fileName}.R',force=TRUE)"
	cp DESCRIPTION ${fileName}/DESCRIPTION
	echo Date: `date`>>${fileName}/DESCRIPTION
	R --vanilla -e "library(roxygen);roxygenize('${fileName}',roxygen.dir='${fileName}',copy.package=FALSE)"
	@echo Copying C
	mkdir ${fileName}/src
	cp *.c ${fileName}/src
	rmdir ${fileName}/inst/doc ${fileName}/inst

clean:
	rm -r ${fileName} ${fileName}_${version}.tar.gz


