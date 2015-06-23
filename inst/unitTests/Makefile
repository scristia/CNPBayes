## Note that the tests are run on the scripts in the installed directly.
## Need to run 'make' to test any changes to the files in inst/unitTests
TOP=../..
PKG=${shell cd ${TOP};pwd}
SUITE=doRUnit.R
#R=${R_HOME}/bin/R

all: inst test

inst: # Install package
	cd ${TOP}/..;\
	R CMD INSTALL ${PKG}


##  test all unit tests
test:
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/3.2-bioc-devel/library R-3.2 --vanilla --slave < ${SUITE}

##  test all unit tests
test_jc:
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/.Rlibs R --vanilla --slave < ${SUITE}

## install CNPBayes and then test unit tests
devel:
	R_LIBS_USER=~/Library/R/3.2-bioc-devel/library R-3.2 CMD INSTALL --no-lock ${PKG}
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R_LIBS_USER=~/Library/R/3.2-bioc-devel/library R-3.2 --vanilla --slave < ${SUITE}

## RS: I'll be testing with the following setup
install:
	R_LIBS_USER=~/Library/R/3.2-bioc-devel/library R-3.2 CMD INSTALL --no-lock ${PKG}
