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
trans = np.random.uniform(-10,10,[3])

start_points = np.random.uniform(-1,1,[10,3])
end_points = np.dot(start_points, R.T) + trans
c = center(start_points)

act_rmsd =  rmsd(start_points, end_points)
print rmsd(np.array([np.dot(c, R.T) + trans]), np.array([c]))
print act_rmsd

