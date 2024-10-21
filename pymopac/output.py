import os
from .helpers import optional_imports


class BaseOutput:
    """
    Base Class that takes the outfile as str and perfors various parsing
    operations.
    Subclasses are to be mainly defined via the self.parsers list.
    """

    def __init__(self, outfile: str):
        self.result = self.ResultFromOutput(outfile)
        self.parsers = []

    def ResultFromOutput(self, outfile: str):
        return outfile
        # TODO

    def parse_all(self):
        for parser in self.parsers:
            parser.parse(self, self.result)


class BaseParser:
    """
    Base Class that takes the whole result section, find the relevant section
    and furthermore searches for
        1) a key in the section
        2) a relative result index, as by str.split()
        3) (optional) a relative unit index
    Subclasses are to be mainly defined via specifying self.location_dict and,
    optionally, the get_section methods.
    """

    def __init__(self, result):
        self.section = self.get_section(result)
        self.location_dict = dict()

    def get_section(self, result):
        return result.split()
        # TODO

    def locate(self, search_tuple: tuple):
        start, end = self.find_sublist(search_tuple[0])
        if end is None:
            return None
            # raise Exception("sublist not found")

        if isinstance(search_tuple[1], int):
            number = self.section[end + search_tuple[1]-1]
        elif isinstance(search_tuple[1], list):
            number = " ".join([self.section[end + x-1]
                              for x in search_tuple[1]])
        else:
            raise KeyError("wrong type of search key passed")

        unit = None
        if len(search_tuple) == 3:
            if isinstance(search_tuple[2], int):
                unit = self.section[end + search_tuple[2]-1]
            elif isinstance(search_tuple[2], list):
                unit = " ".join([self.section[end + x-1]
                                 for x in search_tuple[2]])
            else:
                raise KeyError("wrong type of search key passed")

        return NumUnit(number, unit)

    def location_group(self):
        result_dict = dict()
        for key, value in self.location_dict.items():
            result_dict[key] = self.locate(value)
        return result_dict

    def find_sublist(self, search):
        if isinstance(search, str):
            search = search.split()
        search_length = len(search)
        try:
            i = self.section.index(search[0])
            while i < len(self.section):
                if self.section[i:i+search_length] == search:
                    return i, i+search_length
                i = self.section[i+1:].index(search[0])+i+1
        except Exception:
            return None, None
        return None, None

    def set_result(self, outputclass, search_tuple):
        search_result = self.locate(search_tuple)
        if isinstance(search_result, NumUnit):
            key = search_tuple[0].strip("=").strip().replace(" ", "_")
            outputclass[key] = search_result

    def parse(self, outputclass, result):
        pass


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


class MopacOutput(BaseOutput):
    def __init__(self, outfile: str, stdout=None, stderr=None):
        super().__init__(outfile)
        self.outfile = outfile
        try:
            lines = outfile.split("\n")[:2]
            self.header = lines[0].strip()
            self.comment = lines[1].strip()
        except:
            pass
        self.stdout = stdout
        self.stderr = stderr
        self.parsers = [XyzParser(self.result), StandardParser(self.result)]
        self.parse_all()

    def keys(self):
        return self.__dict__.keys()

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value


class XyzParser(BaseParser):
    def parse(self, outputclass, result):
        _, end = self.find_sublist("CARTESIAN COORDINATES")
        result = result.split()[end:]
        i = 0
        line_j = 1
        xyz = ""
        while isinstance(line_j, int):
            if str(line_j) == result[i]:
                line_j += 1
                xyz += " ".join(result[i+1:i+5])+"\n"
                i += 5
            else:
                xyz = f"{line_j-1}\n\n" + xyz
                line_j = None
        if xyz != "0\n\n":
            outputclass.xyz = xyz[:-1]


class StandardParser(BaseParser):
    def parse(self, outputclass, result):
        result = result.split()

        self.set_result(outputclass,
                        ("FINAL HEAT OF FORMATION =", 4, 5))
        self.set_result(outputclass,
                        ("COSMO AREA              =", 1, [2, 3]))
        self.set_result(outputclass,
                        ("COSMO VOLUME            =", 1, [2, 3]))
        self.set_result(outputclass,
                        ("GRADIENT NORM           =", 3, [4, 5]))
        self.set_result(outputclass,
                        ("IONIZATION POTENTIAL    =", 1, 2))

        homo = self.locate(("HOMO LUMO ENERGIES (EV) =", 1))
        if isinstance(homo, NumUnit):
            homo.unit = "EV"
            outputclass["HOMO"] = homo
        lumo = self.locate(("HOMO LUMO ENERGIES (EV) =", 1))
        if isinstance(lumo, NumUnit):
            lumo.unit = "EV"
            outputclass["LUMO"] = lumo

        self.set_result(outputclass,
                        ("NO. OF FILLED LEVELS    =", 1))
        self.set_result(outputclass,
                        ("MOLECULAR WEIGHT        =", 1))
        _, pointgroup_i = self.find_sublist("POINT GROUP:")
        if pointgroup_i:
            outputclass["POINT_GROUP"] = self.section[pointgroup_i]
