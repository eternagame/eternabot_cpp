import os
import subprocess

file_path = os.path.realpath(__file__)
spl = file_path.split("/")
BASE_DIR = "/".join(spl[:-2])

class SequenceDesignerData(object):
    def __init__(self, sequence, score):
        self.sequence, self.score = sequence, score

class SequenceDesigner(object):
    def __init__(self):
        self.__exe_path = BASE_DIR + "/cmake/build/eternabot"
        if not os.path.isfile(self.__exe_path):
            raise IOError("cannot find eternabot exe")

    def __call_cpp_eternabot(self, sequence, structure, steps, solutions):
        cmd = self.__exe_path + " -seq {seq} -ss \"{ss}\" -steps {s} -n {sols}".format(
            seq=sequence, ss=structure, s=steps, sols=solutions)
        output = subprocess.check_output(cmd, shell=True)
        return output

    def __parse_cpp_output(self, output):
        lines = output.decode().split("\n")
        sols = []

        lines.pop()
        for l in lines:
            spl = l.split()
            sols.append(SequenceDesignerData(spl[1],float(spl[0])))
        return sols

    def design(self, structure, sequence, steps=100, solutions=1):
        output = self.__call_cpp_eternabot(sequence, structure, steps, solutions)
        sols = self.__parse_cpp_output(output)
        return sols


def main():
    sd = SequenceDesigner()
    sols = sd.design("((((....))))", "NNNNUUCGNNNN")

if __name__ == "__main__":
    main()
