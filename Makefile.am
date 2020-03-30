AUTOMAKE_OPTIONS = dejagnu
SUBDIRS = src
dist_doc_DATA = README.md
EXTRA_DIST = autogen.sh src test
RUNTESTFLAGS = --tool mpl MPL=`pwd`/src/mpl --srcdir $$srcdir/test

save:
	gitpush

cov:
	make check
	make -C src cov
