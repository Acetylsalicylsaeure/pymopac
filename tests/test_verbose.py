from pymopac import MopacInput


def test_silent(capsys):
    out = MopacInput("CCC", preopt=True, addHs=True, verbose=False).run()

    captured = capsys.readouterr()
    assert captured.out == ""
    assert "ended normally" in out.outfile.lower()


def test_verbose(capsys):
    MopacInput("CCC", preopt=True, addHs=True, verbose=True).run()

    captured = capsys.readouterr()
    print(captured.out)
    assert "created path at" in captured.out
