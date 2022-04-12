import time

def record_timestamp(func):
    """Record timestamp before and after call of `func`."""
    def deco(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()

        return result, end_time - start_time

    return deco
