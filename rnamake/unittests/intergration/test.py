class Edge(object):
    pass

class Vertex(object):
    def __init__(self, item):
        self.item = item

class Graph(object):
    _edge = Edge
    _vertex = Vertex

    def __init__(self):
        self.vertices = []
        self.edges = []

    def add_vertex(self, item):
        self.vertices.append( self._vertex(item) )

class FancyVertex(Vertex):
    def __init__(self, item):
        super(self.__class__, self).__init__(item)
        self.fancy = 1

class FancyGraph(Graph):
    _vertex = FancyVertex


gf = FancyGraph()
gf.add_vertex(1)
print gf.vertices[0].fancy

g = Graph()
g.add_vertex(1)
print g.vertices[0]


