
# docs

## Papers

Current Overleaf submodules linked to:

| Directory             | Overleaf |
| --------------------- | ------------- |
| papers/Overview/ms    | https://www.overleaf.com/16412681zkjtcrcyzrwc |

## Sphinx

We use [Sphinx](http://www.sphinx-doc.org/) and
[https://readthedocs.org/](Read the Docs) to maintain the documentation
for the DAP.  The documentation is rather bare-bones at the moment.

The API documentation is built using
[sphinx-apidoc](https://www.sphinx-doc.org/en/master/man/sphinx-apidoc.html)
as follows:

```
export SPHINX_APIDOC_OPTIONS='members,private-members,undoc-members,show-inheritance'
sphinx-apidoc -o ./docs/sphinx/src python/mangadap --separate -H mangadap -A "SDSS-IV/MaNGA Pipeline Group"
```

