from argparse import ArgumentError, ArgumentTypeError
from typing import Callable


def min_max_float(min_: float, max_: float) -> Callable:
    def _func(arg) -> float:
        try:
            f = float(arg)
        except ValueError:
            raise ArgumentTypeError("Must be of type float")

        if f <= min_ or f >= max_:
            msg = f"value needs to be between {min_} and {max_}"
            raise ArgumentTypeError(msg)

        return f

    return _func
