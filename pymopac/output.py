import warnings


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

    def parseAll(self):
        for parser in self.parsers:
            parser.parse(self)


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
        """
        preformating for the section which gets parsed
        """
        return self.split_with_newlines(result)
        # TODO

    def locate(self, search_tuple: tuple):
        """
        finds search term and returns NumUnit dataclass
        """
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
        """
        method to make it easier to search sublists resulting from str.split()
        returns start and end index
        """
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
        """
        shorthand function that directly sets the search result on the output class
        """
        search_result = self.locate(search_tuple)
        if isinstance(search_result, NumUnit):
            key = search_tuple[0].strip("=").strip().replace(" ", "_")
            outputclass[key] = search_result

    def split_with_newlines(self, text):
        return [word for line in text.splitlines() for word in (line.split() or [''])]

    def parse(self, outputclass):
        warnings.warn("parser not implemented")


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


class ValUnit:
    """
    Class to wrap an arbitrary value with a unit
    """

    def __init__(self, value, unit=None):
        try:
            self.value = float(value.replace("D", "e"))
        except:
            self.value = value
        self.unit = unit

    def __repr__(self):
        if self.unit:
            return f"{self.value} {self.unit}"
        else:
            return str(self.value)


class MopacOutput(BaseOutput):
    """
    Main Output class, calls parsers and structures outputs

    Custom parsers can be written, inherinting from input.BaseParser.
    Every Output class has a list at self.parser which, custom parsers can
    simply be added. If a parser has been added after Output initialization,
    parsing can be redone using Output.parseAll()
    By convention, custom parsers are able to set attributes on the Ouput class
    by having the Output passed as an argument at parse time.

    a .aux file is passed, as is standard when creating the output via the
    MopacInput.run() method. Using this, all properties can be parsed in an
    unsupervised manner. Results of this can be found under self.auxDict
    """

    def __init__(self, outfile: str, stdout=None, stderr=None, aux: str = None):
        super().__init__(outfile)
        self.outfile = outfile
        if aux:
            self.aux = aux
        try:
            lines = outfile.split("\n")[:2]
            self.header = lines[0].strip()
            self.comment = lines[1].strip()
        except:
            self.header = "Not found"
            self.comment = "Not found"
        self.stdout = stdout
        self.stderr = stderr
        self.parsers = [XyzParser(self.result), StandardParser(self.result)]
        if hasattr(self, "aux"):
            self.parsers.append(AuxParser(self.aux))
        self.parseAll()

    def keys(self):
        return self.__dict__.keys()

    def __getitem__(self, key):
        return self.__dict__[key]

    def __setitem__(self, key, value):
        self.__dict__[key] = value

    def toMol(self):
        """
        returns an rdkit.Chem mol object
        """
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        if hasattr(self, "xyz"):
            mol = Chem.MolFromXYZBlock(self.xyz)
            try:
                rdDetermineBonds.DetermineBonds(mol)
            except:
                warnings.warn("failed to infer bond order")
            return mol


class XyzParser(BaseParser):
    def parse(self, outputclass):
        _, end = self.find_sublist("CARTESIAN COORDINATES")
        if isinstance(end, int):
            end += 1
        result = self.section[end:]
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
    def parse(self, outputclass):
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
        _, formula_i = self.find_sublist("Empirical Formula:")
        if formula_i:
            formula_e = self.section[formula_i:].index("=")
            outputclass["Empirical_Formula"] = " ".join(
                self.section[formula_i:formula_i+formula_e])
            if self.section[formula_i+formula_e+2] == "atoms":
                outputclass["Atom_Count"] = self.section[formula_i+formula_e+1]
                if "xyz" in outputclass.keys() and outputclass["xyz"].split("\n")[0] != outputclass["Atom_Count"]:
                    warnings.warn(
                        "Atom count of molecular formula incongruent with xyz block, proceed with caution")


class AuxParser(BaseParser):
    def get_section(self, result):
        return result.split("\n")

    def parse(self, outputclass):
        try:
            assert self.section[0] == " START OF MOPAC PROGRAM"
            assert self.section[1] == " START OF MOPAC FILE"
            assert self.section[-3] == " END OF MOPAC FILE"
            assert self.section[-2] == " END OF MOPAC PROGRAM"
        except Exception as e:
            warnings.warn("AUX assertions failed, proceed with caution")

        self.section = self.section[2:-3]
        main_dic = dict()
        section_dic = dict()
        section_header = ""
        in_header = False

        title = ""
        unit = None
        line_memory = []

        for line in self.section:
            if " ########" in line:
                in_header = not in_header
                if len(section_dic) > 0:
                    main_dic[section_header] = section_dic
                if in_header:
                    section_header = ""
                continue
            if in_header:
                section_header += line.strip().strip("#").strip()
            else:
                if "=" in line:
                    subline = line[:line.index("=")]
                    if section_header != "":
                        if title != "" and len(line_memory) > 0:
                            section_dic[title] = ValUnit(
                                "\n".join(line_memory), unit)
                        if ":" in line:
                            dpoint_index = line.index(":")
                            title = subline[:dpoint_index].strip()
                            unit = subline[dpoint_index+1:].strip("=")
                        else:
                            title = subline.strip().strip("=")
                            unit = None
                        line_memory = []
                        eq_index = line.index("=")+1
                        if len(line) > eq_index:
                            line_memory.append(line[eq_index:].strip())
                else:
                    line_memory.append(line.strip())

        outputclass.auxDict = main_dic
