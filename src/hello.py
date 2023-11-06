def print_name(name: str) -> None:
    """Function that prints the name given in the argument.
    If the given string is empty, it raises an error.
    
    Parameters
    ----------
    name: str
        A string indicating the name to print.
        
    Returns
    -------
    None
    """
    if len(name) < 1:
        raise ValueError("Name given is empty, please enter a valid name.")
    
    print(f"The name is {name}")