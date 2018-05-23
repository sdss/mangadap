
# docs

## Papers

Current Overleaf submodules linked to:

| Directory             | Overleaf |
| --------------------- | ------------- |
| papers/Overview       | https://www.overleaf.com/16412681zkjtcrcyzrwc |
| papers/EmLines        | https://www.overleaf.com/14318910rcfbhbmxnvcy |

### Connecting this repo to an Overleaf project

When preparing new documentation, you can connect this repo to an
Overleaf repo as follows.  (Inspired by the description
[here](https://abyvinod.github.io/gitsubmodules.html).  Beware that
behavior may vary with different git versions.  I'm using version
2.17.0.)

The Overleaf share button will provide the repo location.  For the
Overview paper, the result is:

https://git.overleaf.com/16412681zkjtcrcyzrwc

To add this document as a submodule:

```
cd ~/Work/MaNGA/dap/repo/mangadap
git submodule add https://git.overleaf.com/16412681zkjtcrcyzrwc docs/papers/Overview
cd docs/papers/Overview
git checkout master
```

You can then treat the docs/papers/Overview directory as you would any
other git repository.  **Be sure to pull/push often to make sure that
you're in sync with anyone making edits directly via the Overleaf
browser interface.**

The status of the main repository should also track changes in the
submodule.  I'll update this with details of any complications I run
into.

## Sphinx

The sphinx documentation is built as follows:

 - Add the members to documentation via the environmental keyword

    ```
    export SPHINX_APIDOC_OPTIONS='members,private-members,undoc-members,show-inheritance'
    ```

 - Automatically write the src directory

    ```
    sphinx-apidoc -o ./sphinx/src ../python/mangadap --full --separate -H mangadap -A "SDSS-IV/MaNGA Pipeline Group" -V 2.2.2 -R 2.2.2
    ```

 - Edit conf.py extensions to:

    ```
    extensions = [
        'sphinx.ext.autodoc',
        'sphinx.ext.todo',
        'sphinx.ext.viewcode',
        'sphinx.ext.napoleon',
        'sphinx.ext.mathjax',
        'sphinx.ext.todo',
        'matplotlib.sphinxext.only_directives',
        'matplotlib.sphinxext.plot_directive',
    ]
    html_theme = 'sphinx_rtd_theme'
    ```

 - The "plot_directive" extension allows you to include code that will
   construct and save plots that will be included in the documentation.
   Here's an example of what you can include in the docstring:

    ```
    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        from scipy.interpolate import InterpolatedUnivariateSpline
        x = np.linspace(-3, 3, 50)
        y = np.exp(-x**2) + 0.1 * np.random.randn(50)
        spl = InterpolatedUnivariateSpline(x, y, k=1)
        plt.plot(x, y, 'ro', ms=5)
        xs = np.linspace(-3, 3, 1000)
        plt.plot(xs, spl(xs), 'g', lw=3, alpha=0.7)
        plt.show()
    ```

 - (Re)Build the documentation

    ```
    rm -fR ./sphinx/plot_directive ; sphinx-build -aEb html -w ./sphinx/warnings.out -T ./sphinx/src ./sphinx
    ```

 - Scorched-earth approach:  If you want to completely rebuild the docs,
   you can use the `from_scratch.bash` script, which basically the
   commands above, but with a group of `rm` statements to remove the
   existing sphinx output.
   
   ```
   source from_scratch.bash
   ```
   
   **Warning**: Any files *not* produced by sphinx will be lost; i.e.,
   you'll lose any by-hand changes you've made to any files.  Git should
   be able to then sort out the chaos.  The script could be edited to
   keep any human generated, e.g., `*.rst` files.


