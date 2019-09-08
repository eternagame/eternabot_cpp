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
        subprocess.call(cmd, shell=True)

    def __parse_cpp_output(self, file_name):
        f = open(file_name)
        lines = f.readlines()
        f.close()

        lines.pop(0)
        sols = []

        for l in lines:
            spl = l.split(",")
            sols.append(SequenceDesignerData(spl[2],float(spl[1])))
        return sols

    def design(self, sequence, structure, steps=100, solutions=1):
        self.__call_cpp_eternabot(sequence, structure, steps, solutions)
        sols = self.__parse_cpp_output("eternabot.scores")
        return sols


def main():
    sd = SequenceDesigner()
    sols = sd.design("NNNNUUCGNNNN", "((((....))))")
    
if __name__ == "__main__":
    main()
