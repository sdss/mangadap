# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXAPI     = sphinx-apidoc
SPHINXAPIOPT  = members,private-members,undoc-members,show-inheritance
SPHINXOPTS    = -aE -w ./sphinx_warnings.out
SPHINXBUILD   = sphinx-build
BUILDDIR      = _build
STATICDIR     = _static
DOCTREE       = $(BUILDDIR)/doctrees
LOCALFILES    = $(BUILDDIR)/* api/*.rst sphinx_warnings.out

# Internal variables.
ALLSPHINXOPTS   = -d $(DOCTREE) $(SPHINXOPTS) .

clean:
	rm -rf $(LOCALFILES)

apirst:
	SPHINX_APIDOC_OPTIONS=$(SPHINXAPIOPT) $(SPHINXAPI) --separate -o ./api ../mangadap ../mangadap/tests/* ../mangadap/scripts/* 
	python3 ./scripts/write_bitmask_tables.py
	python3 ./scripts/write_parset_tables.py
	python3 ./scripts/write_datatables.py
	python3 ./scripts/write_defaultplan.py
	python3 ./scripts/write_script_help.py

html:
	make apirst
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

