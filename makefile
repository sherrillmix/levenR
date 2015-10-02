VERSION:=$(shell grep Version: DESCRIPTION|sed 's/Version: //')
NAME:=$(shell grep Package: DESCRIPTION|sed 's/Package: //')
PACKAGEFILE:=../$(NAME)_$(VERSION).tar.gz

all: $(PACKAGEFILE) README.md

.PHONY: all install

install:
	R -e 'devtools::install_github("sherrillmix/$(NAME)")'

man: R/*.R
	R -e 'devtools::document()'
	touch man

$(PACKAGEFILE): man R/*.R DESCRIPTION tests/*.R src/*.c
	sed "s/^Date:.*$/Date: `date +%Y-%m-%d`/" DESCRIPTION>tmp
	cp tmp DESCRIPTION
	R -e 'devtools::check();devtools::build()'
