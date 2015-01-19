import base
import option
import motif_type
import chain
import structure
import motif
import pose
import atom
import util
import math
import numpy as np

class ChainEndPairMap(object):
    def __init__(self, chain1, chain2):
        self.p5_chain, self.p3_chain = chain1, chain2

    def is_hairpin(self):
        if self.p5_chain == self.p3_chain:
            return 1
        else:
            return 0

    def chains(self):
        return [self.p5_chain, self.p3_chain]


class MotifTreeMerger(base.Base):
    def __init__(self, **options):
        self.seen_connections, self.chains, self.nodes = {}, [], []
        self.setup_options_and_constraints()

    def setup_options_and_constraints(self):
        options = { 'include_head'         : 0,
                    'chain_closure'        : 1,
                    }
        self.options = option.Options(options)
        self.constraints = {}

    def merge(self, mt, **options):
        self.options.dict_set(options)

        # TODO turn motif into pose
        if len(mt.nodes) == 2 and self.option('include_head') == 0:
            return mt.nodes[1].motif.copy()

        self.seen_constraints, self.chains, self.nodes  = {}, [], mt.nodes

        for i, n in enumerate(self.nodes):
            if i == 0 and self.option('include_head') == 0:
                continue
            self.chains.extend([c.subchain(0) for c in n.motif.chains()])

        start_node = self.nodes[0]
        if self.option('include_head') == 0:
            start_node = self.nodes[1]

        self._merge_chains_in_node(start_node)
        motif = self._build_pose()
        return motif

    def reset(self):
        self.chains = []
        self.nodes = []
        self.seen_connections = {}
    def _build_pose(self):
        new_structure = structure.Structure()
        new_structure.assembled = 1
        for c in self.chains:
            new_structure.chains.append(c.copy())
        new_structure.renumber()

        residues = new_structure.residues()
        basepairs = []

        uuids = {}
        for res in residues:
            uuids[res.uuid] = res

        new_pose = pose.Pose()
        for node in self.nodes:
            for bp in node.motif.basepairs:
                if bp.res1.uuid in uuids and bp.res2.uuid in uuids:
                    cbp = bp.copy()
                    cbp.res1 = uuids[bp.res1.uuid]
                    cbp.res2 = uuids[bp.res2.uuid]

                    if node.motif.mtype == motif_type.HELIX:
                        new_pose.designable[cbp.uuid] = 1

                    basepairs.append(cbp)

        new_pose.name = "assembled"
        new_pose.structure = new_structure
        new_pose.basepairs = basepairs
        new_pose.setup_basepair_ends()

        if self.option('chain_closure'):
            for i,c in enumerate(new_pose.chains()):
                close_chain(c)

        return new_pose

    def _merge_chains_in_node(self, node):
        for c in node.connections:
            if c in self.seen_connections:
                continue
            self.seen_connections[c] = 1
            partner = c.partner(node)
            if not self.option('include_head') and partner == self.nodes[0]:
                continue
            #print node.index, partner.index
            node_chains    = self._get_chains_from_connection(node, c)
            partner_chains = self._get_chains_from_connection(partner, c)
            # TODO maybe figure out which is a better basepair to remove
            if partner.motif.mtype == motif_type.HELIX:
                merged_chains = self._helix_merge(node_chains, partner_chains)
            else:
                merged_chains = self._non_helix_merge(node_chains, partner_chains)
            used_chains = node_chains.chains() + partner_chains.chains()
            new_chains = []
            for c in self.chains:
                if c not in used_chains:
                    new_chains.append(c)
            for c in merged_chains:
                if c is not None:
                    new_chains.append(c)
            self.chains = new_chains
            self._merge_chains_in_node(partner)

    def _helix_merge(self, nc, pc):
        merged_chain_1, merged_chain_2 = None, None
        if   nc.is_hairpin() and pc.is_hairpin():
            raise ValueError("cannot merge an hairpin with another hairpin")
        elif nc.is_hairpin():
            p3_chain = pc.p3_chain.subchain(0,-1)
            p5_chain = pc.p5_chain.subchain(1)
            merged_chain_1 = self._get_merged_hairpin(p3_chain, p5_chain,
                                                      nc.p5_chain)
        elif pc.is_hairpin():
            merged_chain_1 = self._get_merged_hairpin(nc.p5_chain, nc_p3_chain,
                                                      pc.p5_chain, 1, 1)
        else:
            merged_chain_1 = self._get_merged_chain(nc.p5_chain, pc.p3_chain,
                                                    1, 1)
            merged_chain_2 = self._get_merged_chain(nc.p3_chain, pc.p5_chain,
                                                    0, 1)
        return merged_chain_1,merged_chain_2

    def _non_helix_merge(self, nc, pc):
        merged_chain_1, merged_chain_2 = None, None
        p3_chain, p5_chain = nc.p3_chain.subchain(0,-1), nc.p5_chain.subchain(1)
        if   nc.is_hairpin() and pc.is_hairpin():
            raise ValueError("cannot merge an hairpin with another hairpin")
        elif nc.is_hairpin():
            merged_chain_1 = self._get_merged_hairpin(pc.p3_chain, pc.p5_chain,
                                                      p5_chain)
        elif pc.is_hairpin():
            merged_chain_1 = self._get_merged_hairpin(p5_chain, p3_chain,
                                                      pc.p5_chain, 1)
        else:
            merged_chain_1 = self._get_merged_chain(p5_chain, pc.p3_chain, 1)
            merged_chain_2 = self._get_merged_chain(p3_chain, pc.p5_chain)
        return merged_chain_1, merged_chain_2

    def _get_merged_chain(self, c1, c2, join_by_3prime=0, remove_overlap=0):
        """
        Merges two chains together that share a common resiude

        :param c1: chain 1
        :param c2: chain 2
        :param join_by_3_prime: joins in the 3prime direction instead of the
            standard (optional)
        :param remove_overlap: removes the overlap residue between the two
            chains (optional)
        :type c1: chain object

        Example:

        chain1            chain2
        5'_|_|_|_|_|_3' + 5'_|_|_|_|_|_3' =
        5'_|_|_|_|_|_|_|_|_|_|_ 3'

        Notice one residue is lost at 5' end of the chain2, this is the overlap
        residue which was used to align the chains during alignment, you can
        set remove_overlap to 0 to stop that from happning

        """

        merged_chain = chain.Chain()
        chain1_res, chain2_res = c1.residues, c2.residues
        if join_by_3prime:
            chain1_res, chain2_res = chain1_res[::-1], chain2_res[::-1]
        merged_chain.residues = list(chain1_res)
        if remove_overlap:
            chain2_res.pop(0)
        merged_chain.residues.extend(list(chain2_res))
        if join_by_3prime:
            merged_chain.residues = merged_chain.residues[::-1]
        return merged_chain

    def _get_merged_hairpin(self, c1, c2, hairpin, join_by_3prime=0,
                            remove_overlap=0):
        hairpin.to_pdb("hairpin.pdb")
        merged_chain = self._get_merged_chain(c1, hairpin, join_by_3prime,
                                              remove_overlap)
        merged_chain = self._get_merged_chain(merged_chain, c2, join_by_3prime,
                                              remove_overlap)
        return merged_chain

    def _get_chains_from_connection(self, node, c):
        end = c.motif_end(node)
        return self._find_chains_for_end(end)

    def _find_chains_for_end(self,end):
        # TODO do I need this complicated function, not sure anymore
        chain_info = []
        seen = {}
        ci_index = 0

        for i,c in enumerate(self.chains):
            for res in end.residues():
                if res not in seen:
                    seen[res] = []

                if   res == c.first():
                    chain_info.append([c,0,i])
                    seen[res].append(ci_index)
                    ci_index += 1
                elif res == c.last():
                    chain_info.append([c,1,i])
                    seen[res].append(ci_index)
                    ci_index += 1

        if len(chain_info) > 2 and len(seen.values()) == 2:
            seen_chain = {}
            for ci in chain_info:
                if ci[0] not in seen_chain:
                    seen_chain[ci[0]] = []
                seen_chain[ci[0]].append(ci)

            for k,v in seen_chain.iteritems():
                if len(v) == 2:
                    return v

        if len(chain_info) != 2:
            print len(chain_info),end.res1,end.res2,end.res1.uuid
            #for i,c in enumerate(self.chains):
            #    chain_to_pdb("chain."+str(i)+".pdb",c)" " " "" " " " " " " "" " " " " " " "" " "" " " " " " " "" " " " " " " "" " " " " " " "" " "" " " " " " " "" " " " " " " "" " " " " " " "" " "" " " " " " " "" " " " " " " "" " " " " " " "" " " " " " " "
            for ci in chain_info:
                print ci[0],ci[0].first(),ci[0].last()
            raise ValueError("Could not find chain for end")
        elif len(chain_info) != 2 and not error:
            return None

        chain_info.sort(key=lambda x:x[1])

        chain_map = ChainEndPairMap(chain_info[0][0], chain_info[1][0])

        return chain_map


def create_coord_system(atoms):
    assert(len(atoms) == 3)
    e1 = util.normalize(atoms[0].coords - atoms[1].coords)
    e3 = np.cross(e1, atoms[2].coords - atoms[1].coords)
    e3 = util.normalize(e3)
    e2 = np.cross(e3, e1)
    matrix = np.array([e1, e2, e3]).T
    return matrix


def virtual_atom(name, l, theta, phi, parent_atoms):
    # get radians
    theta = math.radians(theta)
    phi   = math.radians(phi)
    matrix = create_coord_system(parent_atoms)
    v = [l*math.cos(theta),
         l*math.sin(theta)*math.cos(phi),
         l*math.sin(theta)*math.sin(phi)]
    coords = np.dot(matrix, v) + parent_atoms[0].coords
    return atom.Atom(name, coords)


def get_projection(coord, current_pos, projection_axis):
    r = coord - current_pos
    return r - np.dot(r, projection_axis)*projection_axis


def rotation_matrix(angle, direction, point=None):
    """Return matrix to rotate about axis defined by point and direction.

    >>> R = rotation_matrix(math.pi/2, [0, 0, 1], [1, 0, 0])
    >>> np.allclose(numpy.dot(R, [0, 0, 0, 1]), [1, -1, 0, 1])
    True
    >>> angle = (random.random() - 0.5) * (2*math.pi)
    >>> direc = np.random.random(3) - 0.5
    >>> point = np.random.random(3) - 0.5
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(angle-2*math.pi, direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(-angle, -direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> I = np.identity(4, numpy.float64)
    >>> np.allclose(I, rotation_matrix(math.pi*2, direc))
    True
    >>> np.allclose(2, numpy.trace(rotation_matrix(math.pi/2,
    ...                                               direc, point)))
    True

    """
    sina = math.sin(angle)
    cosa = math.cos(angle)
    direction = util.normalize(direction[:3])
    # rotation matrix around unit vector
    R = np.diag([cosa, cosa, cosa])
    R += np.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += np.array([[ 0.0,         -direction[2],  direction[1]],
                      [ direction[2], 0.0,          -direction[0]],
                      [-direction[1], direction[0],  0.0]])
    M = np.identity(4)
    M[:3, :3] = R
    if point is not None:
        # rotation not around origin
        point = np.array(point[:3], dtype=numpy.float64, copy=False)
        M[:3, 3] = point - np.dot(R, point)
    return M


def close_torsion(which_dir, parent_atoms, daughter_atoms, match_atoms_1,
                  match_atoms_2):
    matrix = create_coord_system(parent_atoms)
    x = matrix.T[0]
    current_atom_xyz = parent_atoms[0].coords
    weighted_sine = 0
    weighted_cosine = 0
    for i in range(len(match_atoms_1)):
        rho1 = get_projection(match_atoms_1[i].coords, current_atom_xyz, x)
        rho2 = get_projection(match_atoms_2[i].coords, current_atom_xyz, x)
        current_sine = which_dir * np.dot(x, np.cross(rho1,rho2))
        current_cosine = np.dot(rho1,rho2)
        weighted_sine += current_sine
        weighted_cosine += current_cosine

    twist_torsion = math.atan2(weighted_sine, weighted_cosine)
    for daughter_atom in daughter_atoms:
        daughter_atom.coords = np.dot(rotation_matrix( twist_torsion , x )[:3,:3],
                                      ( daughter_atom.coords - current_atom_xyz )) + \
                               current_atom_xyz

# Copied from Rhiju's Rosetta code: protocols/farna/RNA_LoopCloser.cc
def close_chain(chain):
    for i in range(0,len(chain.residues)-1):
        res1 = chain.residues[i]
        res2 = chain.residues[i+1]

        if res1.connected_to(res2,cutoff=2.0):
            continue

        atoms = []
        res1_atoms = [res1.get_atom(name) for name in [ "C4'", "C3'", "O3'" ]]
        res2_atoms = [res2.get_atom(name) for name in "C5',O5',OP2,OP1".split(",")]
        atoms.extend(res1_atoms)
        atoms.extend(res2_atoms)
        fail=0
        for a in atoms:
            if a is None:
                fail=1
        if fail:
            continue

        ovl1 = virtual_atom("OVL1", 1.606497, 60.314519, 0.0,
                            [res1.get_atom("O3'"), res1.get_atom("C3'"),
                             res1.get_atom("C4'")])
        ovl2 = virtual_atom("OVL2", 1.593180, 71.059360, 0.0,
                            [ovl1, res1.get_atom("O3'"),
                             res1.get_atom("C3'")])
        ovu1 = virtual_atom("OVU1", 1.593103, 71.027062, 114.600417,
                            [res2.get_atom("P"), res2.get_atom("O5'"),
                            res2.get_atom("OP2")])

        res1.atoms.extend([ovl1, ovl2])
        res2.atoms.append(ovu1)

        match_atoms_1 = [res1.get_atom("O3'"), ovl1,               ovl2                ]
        match_atoms_2 = [ovu1,                 res2.get_atom("P"), res2.get_atom("O5'")]

        for i in range(100):
            close_torsion(+1,
                              [res1.get_atom("O3'"), res1.get_atom("C3'"),res1.get_atom("C4'")],
                              [ovl1,ovl2], match_atoms_1, match_atoms_2)
            close_torsion(+1,
                              [ovl1,res1.get_atom("O3'"),res1.get_atom("C3'")],
                              [ovl2], match_atoms_1, match_atoms_2)
            close_torsion(-1,
                              [res2.get_atom("P"),res2.get_atom("O5'"),res2.get_atom("C5'")],
                              [ovu1,res2.get_atom("OP1"),res2.get_atom("OP2")],
                              match_atoms_1, match_atoms_2)
            close_torsion(-1,
                              [res2.get_atom("O5'"),res2.get_atom("C5'"),res2.get_atom("C4'")],
                              [ovu1,res2.get_atom("OP1"),res2.get_atom("OP2"),res2.get_atom("P")],
                              match_atoms_1, match_atoms_2)
            close_torsion(-1,
                              [res2.get_atom("C5'"),res2.get_atom("C4'"),res2.get_atom("C3'")],
                              [ovu1,res2.get_atom("OP1"),res2.get_atom("OP2"),
                               res2.get_atom("P"),res2.get_atom("O5'")],
                              match_atoms_1, match_atoms_2)

        res1.atoms.pop()
        res1.atoms.pop()
        res2.atoms.pop()


