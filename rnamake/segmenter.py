import heapq
import pose


class PriorityQueue(object):
    def __init__(self):
        self.elements = []

    def empty(self):
        return len(self.elements) == 0

    def put(self, item, priority):
        heapq.heappush(self.elements, (priority, item))

    def get(self):
        return heapq.heappop(self.elements)[1]


class Pair(object):

    __slots__ = ["res1", "res2", "dist"]

    def __init__(self, res1, res2, dist):
        self.res1 = res1
        self.res2 = res2
        self.dist = dist

    def contains(self, res):
        if res == self.res1 or res == self.res2:
            return 1
        else:
            return 0


class PairSearchNode(object):
    def __init__(self, pairs):
        self.pairs = pairs
        self.score = 0
        for p in self.pairs:
            self.score += p.dist

    def contains(self, res):
        for p in self.pairs:
            if p.contains(res):
                return 1
        return 0


class PairSearch(object):
    def __init__(self):
        self.queue = PriorityQueue()

    def _default_values(self, res, pairs, end_pairs):
        values = {}
        all_pairs = []
        all_pairs.extend(pairs)
        all_pairs.extend(end_pairs)
        for r in res:
            total = 0.0
            count = 0.0
            for p in all_pairs:
                if p.contains(r):
                    total += p.dist
                    count += 1
            values[r] = (total / count) / 2
        return values

    def search(self, res, pairs, end_pairs):
        self.values = self._default_values(res, pairs, end_pairs)
        all_pairs = []
        all_pairs.extend(pairs)
        all_pairs.extend(end_pairs)
        self.res = res

        for p in pairs:
            node = PairSearchNode([p])
            estimated_score = self._estimated_score(node)
            self.queue.put(node, node.score + estimated_score)

        solutions = []
        while not self.queue.empty():
            current = self.queue.get()
            can_add = []
            for p in all_pairs:
                if p in current.pairs:
                    continue
                if current.contains(p.res1) or current.contains(p.res2):
                    continue
                can_add.append(p)

            for p in can_add:
                ps = current.pairs[:]
                ps.append(p)
                node = PairSearchNode(ps)
                estimated_score = self._estimated_score(node)
                if estimated_score == 0:
                    solutions.append(node)
                    if len(solutions) > 10:
                        return solutions
                else:
                    self.queue.put(node, node.score + estimated_score)

        return solutions

    def _estimated_score(self, node):
        estimated_score = 0
        for r in self.res:
            if node.contains(r):
                continue
            estimated_score += self.values[r]
        return estimated_score


class Segments(object):
    def __init__(self, remaining=None, removed=None):
        self.remaining = remaining
        self.removed = removed


class Segmenter(object):
    def __init__(self):
        pass

    def _get_pairs(self, m, res):
        pairs = []
        end_pairs = []
        for c in m.chains():
            for i, res1 in enumerate(res):
                for j, res2 in enumerate(res):
                    if i >= j:
                        continue
                    sc = c.subchain(start_res=res1, end_res=res2)
                    if sc is None:
                        continue
                    dist = len(sc.residues)
                    pairs.append(Pair(res1, res2, dist))
                sc1 = c.subchain(start_res=res1, end_res=c.last())
                sc2 = c.subchain(start_res=c.first(), end_res=res1)
                if sc1 is None:
                    continue
                if res1 != c.last() and c.last() not in res:
                    end_pairs.append(Pair(res1, c.last(), len(sc1.residues)))
                if res1 != c.first() and c.first() not in res:
                    end_pairs.append(Pair(c.first(), res1, len(sc2.residues)))
        return pairs, end_pairs


    def _get_subchain(self, m, pair):
        for c in m.chains():
            sc = c.subchain(start_res=pair.res1, end_res=pair.res2)
            if sc is None:
                continue
            return sc
        raise ValueError("cannot create subchain")

    def _get_pose(self, org_pose, res, bps):
        p = pose.Pose()
        p.structure._build_chains(res)
        p.structure._cache_coords()
        p.basepairs = bps
        p._cache_basepair_frames()
        p.setup_basepair_ends()
        motifs = pose.Motifs()
        for m in org_pose.all_motifs():
            fail = 0
            for r in m.residues():
                if r not in res:
                    fail = 1
                    break
            if fail:
                continue
            motifs._assign_motif_by_type(m)
        p.motifs = motifs
        return p

    def _get_segments(self, m, res, bps, cutpoints):
        removed = self._get_pose(m, res[:], bps)
        other_res = cutpoints[:]
        other_bps = []
        for r in m.residues():
            if r not in res:
                other_res.append(r)
        for bp in m.basepairs:
            if bp.res1 in other_res and bp.res2 in other_res:
                other_bps.append(bp)
        remaining = self._get_pose(m, other_res, other_bps)
        segments = Segments(remaining, removed)
        return segments

    def apply(self, m, bps):
        res = []
        for bp in bps:
            res.extend(bp.residues())

        pairs, end_pairs = self._get_pairs(m, res)

        pair_search = PairSearch()
        solutions = pair_search.search(res, pairs, end_pairs)

        for s in solutions:
            subchains = []
            sub_res = []
            for p in s.pairs:
                sc = self._get_subchain(m, p)
                subchains.append(sc)
                sub_res.extend(sc.residues)
            basepairs = []
            missed_bps = 0
            for bp in m.basepairs:
                if   bp.res1 in res and bp.res2 in res:
                    basepairs.append(bp)
                elif bp.res1 in res and bp.res2 not in res:
                    missed_bps += 1
                elif bp.res1 not in res and bp.res2 in res:
                    missed_bps += 1
            if missed_bps > 2:
                continue
            return self._get_segments(m, sub_res, basepairs, res)







