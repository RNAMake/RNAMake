from rnamake import transformations as t
import numpy as np
import math

def rmsd(p1, p2):
    return np.linalg.norm(p1 - p2) / np.sqrt(len(p1))

def center(points):
    length = points.shape[0]
    sum_x = np.sum(points[:, 0])
    sum_y = np.sum(points[:, 1])
    sum_z = np.sum(points[:, 2])

    return np.array([sum_x/length, sum_y/length, sum_z/length])

R = t.random_rotation_matrix(np.random.uniform(size=[3]))[:3,:3]
R2 = t.random_rotation_matrix(np.random.uniform(size=[3]))[:3,:3]

trans = np.random.uniform(-10,10,[3])

s = np.random.uniform(-1,1,[1,3])


#print R.dot(R2).T
#print R2.T.dot(R.T)

print R.dot(s.T)
exit()

end_points = np.dot(s, R.T) + trans
c = center(s)

act_rmsd =  rmsd(s, end_points)
print rmsd(np.array([np.dot(c, R.T) + trans]), np.array([c]))
print act_rmsd

