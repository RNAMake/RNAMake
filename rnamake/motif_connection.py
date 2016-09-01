import exceptions


class MotifConnection(object):
    def __init__(self, i, j, name_i, name_j):
        self.i, self.j = i, j
        self.name_i, self.name_j = name_i, name_j
        if self.name_i is None:
            self.name_i = ""
        if self.name_j is None:
            self.name_j = ""

    def copy(self):
        return MotifConnection(self.i, self.j, self.name_i, self.name_j)

    def to_str(self):
        s = str(self.i) + "," + str(self.j) + "," + self.name_i + "," + self.name_j
        return s


class MotifConnections(object):
    def __init__(self):
        self.connections = []

    def __len__(self):
        return len(self.connections)

    def __iter__(self):
        return iter(self.connections)

    def add_connection(self, i, j, name_i, name_j):
        c = MotifConnection(i, j, name_i, name_j)
        self.connections.append(c)

    def remove_connections_to(self, i):
        pos = 0
        found = 0
        while pos < len(self.connections):
            if i == self.connections[pos].i:
                self.connections.pop(pos)
                found = 1
                pos -= 1
            elif i == self.connections[pos].j:
                self.connections.pop(pos)
                found = 1
                pos -= 1
            pos += 1

        if found == 0:
            raise exceptions.MotifConnectionException(
                "tried to remove connection to node index: " + str(i) + " but "
                "none were present to remove")

    def in_connection(self, index, name):
        for c in self.connections:
            if c.i == index and c.name_i == name:
                return 1
            if c.j == index and c.name_j == name:
                return 1
        return 0


def str_to_motif_connection(s):
    spl = s.split(",")
    return MotifConnection(int(spl[0]), int(spl[1]), spl[2], spl[3])
