

class BaseOutput:
    """
    Base Class that takes the outfile as str and perfors various parsing
    operations
    """

    def __init__(self, outfile: str):
        self.result = self.ResultFromOutput(outfile)
        self.parsers = []
        self.parse()

    def ResultFromOutput(self, outfile: str):
        pass

    def parse(self):
        for parser in self.parsers:
            parser(self.result)


class BaseParser:
    """
    Base Class that takes the whole result section, find the relevant section
    and furthermore searches for
        1) a key in the section
        2) a relative result index, as by str.split()
        3) (optional) a relative unit index
    """

    def __init__(self, result):
        self.section = self.get_section(result)
        self.location_dict = {}
        return self.location_group()

    def get_section(self, result):
        return result

    def locate(self, value: tuple):
        pass
        # TODO

    def location_group(self):
        result_dict = dict()
        for key, value in self.location_dict.items():
            result_dict[key] = self.locate(value)
        return result_dict


class NumUnit:
    """
    Class to wrap a typical calculated value with its unit
    """

    def __init__(self, number, unit=None):
        self.number = float(number)
        self.unit = unit

    def __repr__(self):
        if self.unit:
            return f"{self.number} {self.unit}"
        else:
            return str(self.number)

    def __float__(self):
        return self.number
