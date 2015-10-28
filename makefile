VERSION:=$(shell grep Version: DESCRIPTION|sed 's/Version: //')
NAME:=$(shell grep Package: DESCRIPTION|sed 's/Package: //')
PACKAGEFILE:=$(NAME)_$(VERSION).tar.gz
CURRENTDIR:=$(shell basename `pwd`)


all: ../$(PACKAGEFILE)  README.md

.PHONY: all install

localInstall:
	R -e 'devtools::install()'

install:
	R -e 'devtools::install_github("sherrillmix/$(NAME)")'

README.md: README.Rmd
	R -e 'knitr::opts_chunk$$set(fig.path="README_files/");knitr::knit("README.Rmd")'

man: R/*.R
	R -e 'roxygen2::roxygenize()'
	touch man

../$(PACKAGEFILE): man R/*.R DESCRIPTION tests/*.R src/*.c
	sed -i "s/^Date:.*$$/Date: `date +%Y-%m-%d`/" DESCRIPTION
	cd ..;\
	R CMD build $(CURRENTDIR);\
	R CMD check $(PACKAGEFILE)
