import time
from functools import wraps


def timeit(func):
    @wraps(func)
    def wrapper_timer(*args, **kwargs):
        tic = time.perf_counter()
        value = func(*args, **kwargs)
        toc = time.perf_counter()
        elapsed_time = toc - tic
        print(
            f"Elapsed time for function '{func.__name__}': {elapsed_time:0.4f} seconds"
        )
        return value

    return wrapper_timer
