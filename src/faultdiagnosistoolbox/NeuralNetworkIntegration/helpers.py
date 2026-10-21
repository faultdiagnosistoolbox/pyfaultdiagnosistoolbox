from dataclasses import dataclass

# Constants used when generating a equations causility model
RELATED_VALUE = 2
STATE_VALUE = 3


@dataclass
class SsmodelElement:
    """A structured description of the variable dependencies of a state"""

    xout: list[str]
    xin: list[str]
    zin: list[str]
    elem_type: str  # 'state', 'der', 'out'
    xout_int: list[str]


@dataclass
class GenConfigParams:
    """Input parameters used to map signal names and equation variables to a variable"""

    y_vars: list[str]
    x_vars: list[str]
    signals: list[str] = None  # the names of the variables given in y_var


class FlowDict(dict):
    """Dictionary with specified string function used for writing in config file"""

    def __str__(self):
        return dict_to_str(self, 2)


def dict_to_str(d, indent=0):
    """Returns a string containing a nested dictionary with brackets, colons and indentation"""
    s = "{"
    if not d:
        return s + " }"
    elif not isinstance(d, dict):
        s += "\n" + (" " * indent) + str(d)
    else:
        for idx, (key, value) in enumerate(d.items()):
            s += "\n" + (" " * indent) + str(key) + ": "
            if isinstance(value, dict):
                s += dict_to_str(value, indent + 2)
            elif isinstance(value, list):
                s += list_to_str(value)
            else:
                s += str(value)
            s += ", " if idx < (len(d) - 1) else ""

    return s + "\n" + (" " * (indent - 2)) + "}"


def list_to_str(lst: list):
    """Returns a string containing a nested list with square brackets and commas"""
    s = "[ "
    for idx, elem in enumerate(lst):
        if isinstance(elem, list):
            s += list_to_str(elem)
        else:
            s += str(elem)

        s += ", " if idx < (len(lst) - 1) else " "

    return s + "]"


class QuotedString(str):
    """String class with specified string function used when writing in config file"""

    def __str__(self):
        return ("'") + self + "'"
