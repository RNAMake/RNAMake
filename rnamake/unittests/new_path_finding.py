from rnamake import motif_factory, util, basic_io
from rnamake import transformations as t
import numpy as np
import random

def make_path(origin, vec, ex, ey, ez):
    path = []
    current = origin
    print ex, ey, ez
    direction = util.normalize(vec)
    for i in range(10):
        path.append(current)
        next = current + direction*5
        current = next

    r = t.euler_matrix(ex, ey, ez)
    direction = np.dot(direction, r.T[:3, :3])

    for i in range(10):
        path.append(current)
        next = current + direction*5
        current = next

    return path

def build_complex_path():
    start = motif_factory.factory.ref_motif.ends[0]
    start.flip()

    full_path = make_path(start.d(), start.r()[2], random.uniform(0, 6.2),
                          random.uniform(0, 6.2), random.uniform(0, 6.2))

    for i in range(5):
        s_direction = full_path[-1] - full_path[-2]
        s_pos = full_path[-1]
        path = make_path(s_pos, s_direction, random.uniform(0, 6.2),
                          random.uniform(0, 6.2), random.uniform(0, 6.2))
        full_path.extend(path)

    basic_io.points_to_pdb("full_path.pdb", full_path)
    f = open("full_path.str", "w")
    f.write(basic_io.points_to_str(full_path))
    f.close()



if __name__ == '__main__':
    build_complex_path()
    exit()

    start = motif_factory.factory.ref_motif.ends[0]
    start.flip()
    cx,cy,cz = -6.2, -6.2, -6.2
    count = 0

    f = open("paths.str", "w")
    while cx <= 6.2:
        cy = -6.2
        while cy <= 6.2:
            cz = -6.2
            while cz <= 6.2:
                path = make_path(start.d(), start.r()[2], cx, cy, cz)
                f.write(str(cx) + "|" + str(cy) + "|" + str(cz) + "|" + basic_io.points_to_str(path) + "\n")
                if count > 0 and count < 100:
                    basic_io.points_to_pdb("test."+str(count)+".pdb", path)
                    fn = open("test_"+str(count)+"path.str", "w")
                    fn.write(basic_io.points_to_str(path))
                    fn.close()
                count += 1
                cz += 0.4
            cy += 0.4
        cx += 0.4
    f.close()
