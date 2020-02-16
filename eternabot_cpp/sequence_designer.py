import os
import subprocess

file_path = os.path.realpath(__file__)
spl = file_path.split("/")
BASE_DIR = "/".join(spl[:-2])

class SequenceDesignerData(object):
    def __init__(self, sequence, score, target_structure, folded_structure):
        self.sequence, self.score = sequence, score
        self.target_structure = target_structure
        self.folded_struture = folded_structure

class SequenceDesigner(object):
    def __init__(self):
        self.__exe_path = BASE_DIR + "/cmake/build/eternabot"
        if not os.path.isfile(self.__exe_path):
            raise IOError("cannot find eternabot exe")

    def __call_cpp_eternabot(self, sequence, structure, steps, solutions):
        cmd = [self.__exe_path, '-seq', sequence, '-ss', structure, '-steps', str(steps), '-n', str(solutions)]
        output = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return output.stdout

    def __parse_cpp_output(self, output):
        lines = output.decode().split("\n")
        sols = []

        lines.pop()
        for l in lines:
            spl = l.split()
            sols.append(SequenceDesignerData(spl[2],float(spl[0]), spl[3], spl[4]))
        return sols

    def design(self, structure, sequence, solutions=1, steps=100):
        output = self.__call_cpp_eternabot(sequence, structure, steps, solutions)
        sols = self.__parse_cpp_output(output)
        return sols


def main():
    sd = SequenceDesigner()
    sols = sd.design("((((....))))", "NNNNUUCGNNNN")
    for s in sols:
        print(s.sequence, s.score, s.target_structure, s.folded_struture)

if __name__ == "__main__":
    main()
