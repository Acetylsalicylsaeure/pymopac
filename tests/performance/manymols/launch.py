#!/bin/python
import pymopac
import sys


def test(smi, count):
    infile = pymopac.MopacInput(smi)
    with open("./infile.mop", "w") as file:
        file.write(infile.inpfile())
    for i in range(count):
        infile.run()


if __name__ == "__main__":
    test(sys.argv[1], int(sys.argv[2]))
