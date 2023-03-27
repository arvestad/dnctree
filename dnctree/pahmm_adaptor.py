import importlib.util
def pahmm_available() -> bool:
    """ 
    Check if the pahmm-library is available.
    """
    return bool(importlib.util.find_spec("pahmm"))

