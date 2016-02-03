


class MotifConnection(object):
    def __init__(self, i, j, name_i, name_j):
        self.i, self.j = i, j
        self.name_i, self.name_j = name_i, name_j
        if self.name_i is None:
            self.name_i = ""
        if self.name_j is None:
            self.name_j

    def copy(self):
        return MotifConnection(self.i, self.j, self.name_i, self.name_j)

    def to_str(self):
        s = str(self.i) + "," + str(self.j) + "," + self.name_i + "," + self.name_j
        return s

def str_to_motif_connection(s):
    spl = s.split(",")
    return MotifConnection(int(spl[0]), int(spl[1]), spl[2], spl[3])
