"""Subcommand entry points.

Each task is imported only when it runs, so ``nanoflux prepare`` and
``nanoflux score`` work on machines without the inference dependencies
(liblinear), and ``nanoflux infer`` without the alignment tools.
"""


def infer(args):
    from .infer import main

    return main(args)


def prepare(args):
    from .prepare import main

    return main(args)


def score(args):
    from .score import main

    return main(args)


__all__ = ["infer", "prepare", "score"]
