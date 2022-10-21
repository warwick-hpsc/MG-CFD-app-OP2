#all: 
#	CC -fPIC -c message.cpp
#	CC -fPIC -c simpic.cpp
#	ar rcs libsimpic.a simpic.o message.o
#	CC coupler.cpp -L. -lsimpic
#debug: 
#	CC -g simpic.cpp message.cpp

lib:
	CC -I../src -I../src_op -fPIC -c message.cpp
	CC -I../src -I../src_op -fPIC -c simpic.cpp
	ar rcs libsimpic.a simpic.o message.o

clean:
	rm -f a.out *.o *.dat *.a

# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    = 
SPHINXBUILD   = sphinx-build
SOURCEDIR     = doc
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile all debug clean

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
