from typing import List

def assemble_path(path: List[str]) -> str:
    """
    Assemble a path from a list of strings.
    :param path: List of strings to assemble.
    :type path: List[str]
    :return: Assembled contig.
    """

    if len(path) == 0:
        return ""
    
    contig = path[0]
    if not isinstance(contig, str):
        raise TypeError("Path must be a list of strings.")
    
    for step in path[1:]:
        if not isinstance(step, str):
            raise TypeError("Path must be a list of strings")
        contig += step[-1]

    return contig