import pytest
from pymopac.output import BaseParser, MopacOutput
from rdkit import Chem


def test_find_sublist():
    input = "a b c d e f g"
    target = "c d e"
    parser = BaseParser(input)
    start, end = parser.find_sublist(target)
    print(start, end)
    assert input.split()[start:end] == target.split()


def test_locate():
    result = """
 -------------------------------------------------------------------------------
 PM7
 #



     GEOMETRY OPTIMISED USING EIGENVECTOR FOLLOWING (EF).     
     SCF FIELD WAS ACHIEVED                                   


                              PM7 CALCULATION
                                                       MOPAC v22.1.1 Linux
                                                       Wed Sep 25 20:41:39 2024




          FINAL HEAT OF FORMATION =        -18.19805 KCAL/MOL =     -76.14064 KJ/MOL


          COSMO AREA              =         75.89 SQUARE ANGSTROMS
          COSMO VOLUME            =         57.43 CUBIC ANGSTROMS

          GRADIENT NORM           =          0.93583          =       0.33086 PER ATOM
          IONIZATION POTENTIAL    =         11.965113 EV
          HOMO LUMO ENERGIES (EV) =        -11.965  4.509
          NO. OF FILLED LEVELS    =          7
          MOLECULAR WEIGHT        =         30.0694         POINT GROUP:  D3d 

          MOLECULAR DIMENSIONS (Angstroms)

            Atom       Atom       Distance
            H     6    H     5     3.09342
            H     7    H     4     2.90185:
            H     8    H     3     2.80353
        """
    parser = BaseParser(result)
    target = ("IONIZATION POTENTIAL    =", 1, 2)
    out = parser.locate(target)
    assert out.number == 11.965113


def test_multiitem_locate():
    result = """
 -------------------------------------------------------------------------------
 PM7
 #



     GEOMETRY OPTIMISED USING EIGENVECTOR FOLLOWING (EF).     
     SCF FIELD WAS ACHIEVED                                   


                              PM7 CALCULATION
                                                       MOPAC v22.1.1 Linux
                                                       Wed Sep 25 20:41:39 2024




          FINAL HEAT OF FORMATION =        -18.19805 KCAL/MOL =     -76.14064 KJ/MOL


          COSMO AREA              =         75.89 SQUARE ANGSTROMS
          COSMO VOLUME            =         57.43 CUBIC ANGSTROMS

          GRADIENT NORM           =          0.93583          =       0.33086 PER ATOM
          IONIZATION POTENTIAL    =         11.965113 EV
          HOMO LUMO ENERGIES (EV) =        -11.965  4.509
          NO. OF FILLED LEVELS    =          7
          MOLECULAR WEIGHT        =         30.0694         POINT GROUP:  D3d 

          MOLECULAR DIMENSIONS (Angstroms)

            Atom       Atom       Distance
            H     6    H     5     3.09342
            H     7    H     4     2.90185:
            H     8    H     3     2.80353
        """

    parser = BaseParser(result)
    target2 = ("COSMO VOLUME =", 1, [2, 3])
    out2 = parser.locate(target2)
    assert out2.unit == "CUBIC ANGSTROMS"


def test_xyz():
    result = """

     5       H          0.57936688  *   0.37495992  *  -0.95004855  *
     6       H          2.90323228  *  -0.18866112  *   1.01240523  *
     7       H          2.90289981  *   1.08409546  *  -0.21565979  *
     8       H          2.90296080  *  -0.61571630  *  -0.70313293  *

                             CARTESIAN COORDINATES

   1    C        0.979177621     0.093385431     0.030953656
   2    C        2.503310891     0.093095057     0.030798060
   3    H        0.579354935     0.802186213     0.765270335
   4    H        0.579660431    -0.897703934     0.277635668
   5    H        0.579366881     0.374959918    -0.950048547
   6    H        2.903232280    -0.188661124     1.012405235
   7    H        2.902899813     1.084095458    -0.215659787
   8    H        2.902960799    -0.615716303    -0.703132935


           Empirical Formula: C2 H6  =     8 atoms



      MOLECULAR POINT GROUP   :   D3d 


        """
    output = MopacOutput(result)
    print("printing xyz", output.xyz)
    assert isinstance(output.xyz, str)
    assert len(output.xyz.split("\n")) == 10
    mol = Chem.MolFromXYZBlock(output.xyz)
    assert mol.GetNumAtoms() == int(output.xyz.split("\n")[0])
