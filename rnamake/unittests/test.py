from rnamake import transformations as t
from rnamake import util, basic_io
import rnamake.resource_manager as rm

import numpy as np
import numpy.linalg as lin

import math

def rmsd(p1, p2):
    return np.linalg.norm(p1 - p2) / np.sqrt(len(p1))

def center(points):
    length = len(points)
    sum_x = np.sum(points[:, 0])
    sum_y = np.sum(points[:, 1])
    sum_z = np.sum(points[:, 2])

    return np.array([sum_x/length, sum_y/length, sum_z/length])

def get_moment_tensor(coords):

    sum = np.zeros([3,3])
    #get center of coordinates
    center = util.center_points(coords)
    for c in coords:
        #difference from center
        diff = c - center
        c_T = diff[np.newaxis].T
        c = diff[np.newaxis]
        sum +=  c_T.dot(c)
    sum /= len(coords)
    eig, vec = lin.eig(sum)
    return eig, vec, center, len(coords)


def part_1_and_2():
    #R**2 * r**2 = rT*RT*R*r
    s_rotated = R*s
    #part 1
    print np.dot(s_rotated.T, s_rotated)
    #print s.T*R.T*R*s
    #part 2
    #print s.T*s

def get_coords(res):
    coords = []
    for i, a in enumerate( res.atoms):
        if i < 12:
            continue
        #print a.name,
        coords.append(a.coords)
    #print
    return coords

def convert_to_principal(d, r):
    d_principal = s_principal_in_conventional_frame.dot(r.T) + d
    r_principal = util.unitarize(r_principal_in_conventional_frame.dot(r.T))
    return d_principal, r_principal

motif = rm.manager.get_motif(name='HELIX.IDEAL')

#orientation of basepair 1
motif.ends[0].flip()
r1 = motif.ends[0].r()
r1 = util.unitarize(r1)
#orientation of basepair 2
motif.ends[1].flip()
r2 = motif.ends[1].r()
r2 = util.unitarize(r2)

d1 = motif.ends[0].d()
d2 = motif.ends[1].d()

#rotaton between r1 and r2
R = util.unitarize(r1.T.dot(r2))
#distance between center of two basepairs
#vector between betweeen both centers
#diff = motif.ends[0].base_d() - motif.ends[1].base_d()
diff = d1 - d2
trans = lin.norm(diff)

#array of coordinates of atoms in bases
coords = np.array(get_coords(motif.ends[0].res1) + get_coords(motif.ends[0].res2))
coords_2 = np.array(get_coords(motif.ends[1].res2) + get_coords(motif.ends[1].res1))

e_eig, e_vec, e_cent, e_length  = get_moment_tensor(coords)


r_principal_in_conventional_frame = e_vec.T.dot(r1)
s_principal_in_conventional_frame = (e_cent - d1).dot(r1)

transformed_coords = [ np.dot((c - d1),r1).dot(r2.T) + d2 for c in coords]
e_eig2, e_vec2, e_cent2, e_length2 = get_moment_tensor(transformed_coords)
print e_vec2
print e_cent2
print convert_to_principal(d2, r2)[1]
print e_eig
print e_eig2

#transformed_coords = [np.dot(c, R.T) + diff for c in coords]
print rmsd(coords, coords_2)
print rmsd(coords_2, transformed_coords)
print rmsd(coords, transformed_coords)
print e_cent, e_vec
print convert_to_principal(d1, r1)
d1_principal, r1_principal = convert_to_principal(d1, r1)
d2_principal, r2_principal = convert_to_principal(d2, r2)

diff = d1_principal - d2_principal
R_principal = util.unitarize(r2_principal.T.dot(r1_principal))
R_principal = e_vec.T.dot(R_principal).dot(e_vec)
print R_principal
print math.sqrt(diff.T.dot(diff) + 2*(np.sum(e_eig) - np.trace(np.dot(R_principal, np.diag(e_eig)))))


#basic_io.points_to_pdb("start.pdb", coords)
#basic_io.points_to_pdb("transform.pdb", transformed_coords)
#basic_io.points_to_pdb("final.pdb", coords_2)


exit()

centered_coords = [ c - e_cent for c in coords]
transformed_centered_coords = [np.dot(c, R.T) + diff for c in centered_coords ]

R_principal = e_vec.T.dot(R).dot(e_vec)
print math.sqrt(diff.T.dot(diff) + 2*(np.sum(e_eig) - np.trace(np.dot(R_principal, np.diag(e_eig)))))
print "rmsd", rmsd(np.array(centered_coords), np.array(transformed_centered_coords))
#first term
#second term is the contribution of the distance from the origin
print np.sum(e_eig)

s = []
for c in centered_coords:
    s.append(c)

sum_rotated = 0
sum_rotated_R = 0
sum = 0
for e in s:
    e_rotated = np.dot(R,e)
    sum_rotated += np.dot(e_rotated.T, e_rotated)
    sum_rotated_R +=  e.T.dot(R.T).dot(R).dot(e)
    sum += e.T.dot(e)
print sum_rotated/len(coords)
print sum_rotated_R/len(coords)
print sum/len(coords)



exit()
#c = center(s)
#s = np.matrix(s)

#print R.dot(R2).T
#print R2.T.dot(R.T)

s_rotated = R*s
print np.dot(s_rotated.T, s_rotated)
print np.trace(np.dot(s_rotated.T, s_rotated))
exit()

end_points = np.dot(s, R.T) + trans
c = center(s)

act_rmsd =  rmsd(s, end_points)
print rmsd(np.array([np.dot(c, R.T) + trans]), np.array([c]))
print act_rmsd

