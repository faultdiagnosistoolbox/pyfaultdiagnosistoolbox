from dataclasses import dataclass

RELATED_VALUE = 2
STATE_VALUE = 3


@dataclass
class SsmodelElement:
    xout: list[str]
    xin: list[str]
    zin: list[str]
    elem_type: str  # 'state', 'der', 'out'
    xout_int: list[str]


class FlowDict(dict):
    pass


def flow_dict_representer(dumper, data):
    """
    Formating in .yml for dictionaries
    """
    return dumper.represent_mapping("tag:yaml.org,2002:map", data, flow_style=True)


class QuotedString(str):
    pass


def quoted_representer(dumper, data):
    """
    Formating in .yml for quoted strings
    """
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"')
