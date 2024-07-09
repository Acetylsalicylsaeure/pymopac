import pymopac


def test_k2co3():
    infile = pymopac.MopacInput("C(=O)([O-])[O-].[K+].[K+]")
    outfile = infile.run()
    assert outfile.stderr == b""
