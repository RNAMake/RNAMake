from . import base
from . import option
from . import motif_type
from . import chain
from . import structure
from . import motif
from . import pose

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

        for node in self.nodes:
            for bp in node.motif.basepairs:
                if bp.res1.uuid in uuids and bp.res2.uuid in uuids:
                    cbp = bp.copy()
                    cbp.res1 = uuids[bp.res1.uuid]
                    cbp.res2 = uuids[bp.res2.uuid]

                    if node.motif.mtype == motif_type.HELIX:
                        cbp.designable = 1

                    basepairs.append(cbp)

        new_pose = pose.Pose()
        new_pose.name = "assembled"
        new_pose.structure = new_structure
        new_pose.basepairs = basepairs
        new_pose.setup_basepair_ends()

        # if self.chain_closure:
        #    for i,c in enumerate(merged_motif.chains):
        #        close_chain(c)
        #        chain_to_pdb("chain."+str(i)+".pdb",c)

        return new_pose

    def _merge_chains_in_node(self, node):
        for c in node.connections:
            if c in self.seen_connections:
                continue
            self.seen_connections[c] = 1
            partner = c.partner(node)
            if not self.option('include_head') and partner == self.nodes[0]:
                continue
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
            merged_chain_1 = self._get_merged_hairpin(p3_chain, p5_chain,
                                                      nc.p5_chain)
        elif pc.is_hairpin():
            merged_chain_1 = self._get_merged_hairpin(nc.p5_chain, nc.p3_chain,
                                                      p5_chain, 1)
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
            #    chain_to_pdb("chain."+str(i)+".pdb",c)
            for ci in chain_info:
                print ci[0],ci[0].first(),ci[0].last()
            raise ValueError("Could not find chain for end")
        elif len(chain_info) != 2 and not error:
            return None

        chain_info.sort(key=lambda x:x[1])

        chain_map = ChainEndPairMap(chain_info[0][0], chain_info[1][0])

        return chain_map


