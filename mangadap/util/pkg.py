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


def load_object(module, object):
    """
    Load an abstracted module and object.

    Thanks to: https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path?rq=1

    Args:
        module (:obj:`str`):
            The name of a global python module, or the root name of a local file
            with the object to import.
        object (:obj:`str`):
            The name of the object to import

    Return:
        :obj:`type`: The imported object.

    Raises:
        ImportError:
            Raised if unable to import ``module``.
    """
    import importlib
    from pathlib import Path

    try:
        CubeModule = importlib.import_module(cube_module)
    except (ModuleNotFoundError, ImportError) as e:
        p = Path(cube_module + '.py').resolve()
        if not p.exists():
            raise ImportError(f'Unable to load module {cube_module}!') from e
        spec = importlib.util.spec_from_file_location(cube_module, str(p))
        CubeModule = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(CubeModule)

    return getattr(CubeModule, cube_object)

