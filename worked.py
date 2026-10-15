"""worked.py - find and solve a small 5-router example both ways; writes worked.json"""
import itertools, json, random, heapq
import networkx as nx

NODES = "ABCDE"
EDGES = [("A", "B"), ("B", "C"), ("C", "D"), ("D", "E"), ("A", "E"), ("B", "D")]
FIXED = [1, 1, 1, 1, 10, 5]


def build(ws):
    G = nx.Graph()
    for (a, b), w in zip(EDGES, ws):
        G.add_edge(a, b, w=w)
    return G


def unique_paths(G):
    for s in G:
        for t in G:
            if s < t:
                ps = list(nx.all_shortest_paths(G, s, t, weight="w"))
                if len(ps) != 1:
                    return False
    return True


def dv_rounds(G):
    """Synchronous DV, no split horizon (kept simple for teaching).
    Round 0 = neighbours only.  Returns list of tables per round."""
    tab = {u: {u: (0, u)} for u in G}
    for u in G:
        for v in G[u]:
            tab[u][v] = (G[u][v]["w"], v)
    hist = [{u: dict(t) for u, t in tab.items()}]
    while True:
        new = {u: dict(t) for u, t in tab.items()}
        for u in G:
            for v in G[u]:
                w = G[u][v]["w"]
                for d, (c, _) in tab[v].items():
                    if d == u:
                        continue
                    if d not in new[u] or c + w < new[u][d][0]:
                        new[u][d] = (c + w, v)
        if new == tab:
            return hist
        tab = new
        hist.append({u: dict(t) for u, t in tab.items()})


def dijkstra_steps(G, src):
    dist = {n: float("inf") for n in G}; dist[src] = 0
    prev = {}; done = set(); steps = []
    while len(done) < len(G):
        u = min((n for n in G if n not in done), key=lambda n: dist[n])
        done.add(u)
        for v in G[u]:
            if v not in done and dist[u] + G[u][v]["w"] < dist[v]:
                dist[v] = dist[u] + G[u][v]["w"]; prev[v] = u
        steps.append({"settled": u, "dist": {n: (None if dist[n] == float("inf") else dist[n])
                                              for n in sorted(G)}, "prev": dict(prev)})
    return steps, dist, prev


ws = FIXED
G = build(ws)
assert unique_paths(G)
hist = dv_rounds(G)
steps, dist, prev = dijkstra_steps(G, "A")
nhs = {t: nx.shortest_path(G, "A", t, weight="w")[1] for t in "BCDE"}
out = {"edges": [[a, b, w] for (a, b), w in zip(EDGES, ws)],
       "dv_rounds": len(hist) - 1,
       "dv_A": [{d: list(t) for d, t in h["A"].items()} for h in hist],
       "dv_E": [{d: list(t) for d, t in h["E"].items()} for h in hist],
       "dijkstra": steps, "truth_A": {k: v for k, v in dist.items()},
       "next_hops_A": nhs}
# self-check: DV result equals Dijkstra truth for every router
for u in G:
    truth = nx.single_source_dijkstra_path_length(G, u, weight="w")
    got = {d: c for d, (c, _) in hist[-1][u].items()}
    assert truth == got, (u, truth, got)
json.dump(out, open("worked.json", "w"), indent=1)
print(ws, "dv rounds:", len(hist) - 1, "A next hops:", nhs)
