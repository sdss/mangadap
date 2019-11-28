
# docs

This directory primarily contains the
[sphinx](http://www.sphinx-doc.org/en/master/) source files for the
documentation hosted on [ReadTheDocs](https://readthedocs.org/).  We
have also included some of the source scripts and files for the
published DAP papers.

----

To build the docs locally, make sure that you have sphinx installed, use
it to build the html, and then open the main page:

```
pip install sphinx
make clean ; make html
open _build/html/index.html
```

