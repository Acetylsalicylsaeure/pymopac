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
        if start is None or end is None:
            return None

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
    def __init__(self, outfile):
        super().__init__(outfile)
        self.parsers = [XyzParser(self.result), StandardParser(self.result)]
        self.parse_all()


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
        outputclass.xyz = xyz[:-1]


class StandardParser(BaseParser):
    pass


class depMopacOutput():
    """
    reads the MOPAC .out file at the given out_path and parses datapoints

    outputs are available raw under self.outfile or parsed as a dictionary
    mostly in the format key: (value, unit) under self.dic

    some dictionary functionality is directly available for the class, i.e.
    keys can directly queried via self[key]

    standalone runs possible, but calling via MopacInput().run() recommended.
    """

    def __init__(self, out_path: str, stdout=None, stderr=None):
        self.stdout = stdout
        self.stderr = stderr
        if not os.path.isfile(out_path):
            raise Exception("output file not found")

        with open(out_path, "r") as file:
            self.outfile = file.read()

        self.result = ResultFromOutput(self.outfile)
        self.dic = ParseResult(self.result)
        self.mol = ExtractMol(self.result)

    def keys(self):
        return self.dic.keys()

    def values(self):
        return self.dic.values()

    def items(self):
        return self.dic.items()

    def __getitem__(self, key):
        return self.dic[key]


def ResultFromOutput(outfile: str):
    break_str = " -------------------------------------------------------------------------------"
    out_l = outfile.split("\n")
    if break_str in out_l:
        break_index = out_l.index(break_str)+1
    else:
        break_index = 0
    result_l = out_l[break_index:]
    return result_l


def ParseResult(result: str) -> dict:
    """
    takes the result section of the MOPAC outfile, returns targeted parts
    mostly in the format key: (value, unit)
    """
    dic = dict()
    dic["header"] = result[0][1:]
    dic["comment"] = result[1][1:]

    for line in result:
        parsed = ParseLine(line)
        if parsed is not None:
            key, value, unit, store = parsed
            if unit is not None:
                unit = " ".join(unit)
                try:
                    unit = float(unit)
                except Exception as e:
                    pass
                dic[key] = (value, unit)
            else:
                dic[key] = value

    return dic


def ParseLine(line):
    targets = {"FINAL HEAT OF FORMATION =": (-2, None),
               "COSMO AREA              =": (-3, None),
               "COSMO VOLUME            =": (-3, None),
               "GRADIENT NORM           =": (-3, None),
               "IONIZATION POTENTIAL    =": (-2, None),
               "HOMO LUMO ENERGIES (EV) =": (-2, None),
               "NO. OF FILLED LEVELS    =": -1,
               "MOLECULAR WEIGHT        =": 3
               }
    for key in targets.keys():
        if key in line:
            spli = line.split()
            target_key = targets[key]
            if isinstance(target_key, tuple):
                unit_i = target_key[1]
                return (key.strip("=").strip(), float(spli[target_key[0]]),
                        spli[target_key[0]+1:unit_i], None)
            elif isinstance(target_key, int):
                return (key.strip("=").strip(), float(spli[target_key]), None,
                        None)
    return None


def ExtractMol(result: str):
    start_i = None
    for i, n in enumerate(result):
        if "CARTESIAN COORDINATES" in n:
            start_i = i
    if start_i is None:
        return Exception("cartesian coordinates not found")
    start_i += 2
    end_i = result[start_i:].index("")
    end_i += start_i
    XYZRaw = result[start_i:end_i]
    XYZBlock = str(len(XYZRaw)) + "\n\n"

    for line in XYZRaw:
        XYZBlock += " ".join(line.split()[1:]) + "\n"

    optional_imports("from rdkit import Chem", globals())
    mol = Chem.MolFromXYZBlock(XYZBlock)

    return mol
