"""
General package utilities.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

def all_subclasses(cls):
    """
    Collect all the subclasses of the provided class.

    The search follows the inheritance to the highest-level class.  Intermediate
    base classes are included in the returned set, but not the base class itself.

    Thanks to:
    https://stackoverflow.com/questions/3862310/how-to-find-all-the-subclasses-of-a-class-given-its-name

    Args:
        cls (object):
            The base class

    Returns:
        :obj:`set`: The unique set of derived classes, including any
        intermediate base classes in the inheritance thread.
    """
    return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in all_subclasses(c)])


def load_object(module, obj=None):
    """
    Load an abstracted module and object.

    Thanks to: https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path?rq=1

    Args:
        module (:obj:`str`):
            The name of a global python module, the root name of a local file
            with the object to import, or the full module + object type.  If
            ``obj`` is None, this *must* be the latter.
        obj (:obj:`str`, optional):
            The name of the object to import.  If None, ``module`` must be the
            full module + object type name.

    Return:
        :obj:`type`: The imported object.

    Raises:
        ImportError:
            Raised if unable to import ``module``.
    """
    import importlib
    from pathlib import Path

    if obj is None:
        _module = '.'.join(module.split('.')[:-1])
        obj = module.split('.')[-1]
    else:
        _module = module

    try:
        CubeModule = importlib.import_module(_module)
    except (ModuleNotFoundError, ImportError) as e:
        p = Path(module + '.py').resolve()
        if not p.exists():
            raise ImportError(f'Unable to load module {_module}!') from e
        spec = importlib.util.spec_from_file_location(_module, str(p))
        CubeModule = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(CubeModule)

    return getattr(CubeModule, obj)

