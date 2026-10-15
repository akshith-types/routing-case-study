"""
routesim.py - Experiments for the case study "When Routing Protocols Fail".

Compares distance-vector (Bellman-Ford, RIP-style, infinity = 16) with
link-state (reliable flooding + Dijkstra, OSPF/ARPANET-SPF style).

Model: synchronous rounds ("ticks"); one tick = one message hop on every link.
Distance vector uses TRIGGERED updates (a router speaks only when its own
table changed), which is the favourable case for DV.  Every converged state is
checked against ground-truth shortest paths.

Run:  python3 -I routesim.py     (writes results.json and fig_*.png)
"""
import heapq, json, random, statistics
import networkx as nx
import matplotlib
import matplotlib.ticker
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INF = 16  # RIP-style infinity
random.seed(7)
OPS = {"spf": 0}          # instrumentation: Dijkstra edge relaxations


# --------------------------------------------------------------------------
# Distance vector (triggered updates)
# --------------------------------------------------------------------------
class DV:
    """mode in {'plain','split','poison'}"""
    def __init__(self, G, mode="plain"):
        self.mode = mode
        self.G = G.copy()
        self.dist = {u: {u: 0} for u in G}
        self.nh = {u: {u: u} for u in G}
        self.heard = {u: {} for u in G}     # last vector heard from each neighbour
        self.sent = {}                      # last vector sent on each directed link
        self.dirty = set(G)                 # routers that must advertise
        self.msgs = 0                       # vector messages
        self.entries = 0                    # total entries inside those messages
        self.cpu = 0                        # entries processed by receivers

    def advert(self, u, to):
        vec = {}
        for d, c in self.dist[u].items():
            if d != u and self.nh[u].get(d) == to:
                if self.mode == "split":
                    continue
                if self.mode == "poison":
                    vec[d] = INF
                    continue
            vec[d] = c
        return vec

    def recompute(self, u):
        new_d, new_nh = {u: 0}, {u: u}
        for v in self.G[u]:
            w = self.G[u][v]["w"]
            for d, c in self.heard[u].get(v, {}).items():
                self.cpu += 1
                if d == u:
                    continue
                cand = min(INF, c + w)
                best = new_d.get(d, INF + 1)
                if cand < best or (cand == best and cand < INF and self.nh[u].get(d) == v):
                    new_d[d], new_nh[d] = cand, v
        for d in [d for d, c in new_d.items() if c >= INF]:
            del new_d[d]; new_nh.pop(d, None)
        changed = (new_d != self.dist[u]) or (new_nh != self.nh[u])
        self.dist[u], self.nh[u] = new_d, new_nh
        return changed

    def step(self):
        senders, self.dirty = list(self.dirty), set()
        receivers = set()
        for u in senders:
            for v in self.G[u]:
                vec = self.advert(u, v)
                prev = self.sent.get((u, v), {})
                diff = sum(1 for d in set(vec) | set(prev) if vec.get(d) != prev.get(d))
                if diff == 0 and (u, v) in self.sent:
                    continue                     # nothing new for this neighbour
                self.sent[(u, v)] = vec
                self.heard[v][u] = vec
                self.msgs += 1; self.entries += max(diff, 1)   # only changed routes travel
                receivers.add(v)
        for v in receivers:
            if self.recompute(v):
                self.dirty.add(v)
        return bool(self.dirty)

    def link_down(self, a, b):
        self.G.remove_edge(a, b)
        self.heard[a].pop(b, None); self.heard[b].pop(a, None)
        self.sent.pop((a, b), None); self.sent.pop((b, a), None)
        for x in (a, b):
            if self.recompute(x):
                self.dirty.add(x)

    def run(self, max_ticks=500, watch=None):
        trace, t = [], 0
        while self.dirty and t < max_ticks:
            self.step(); t += 1
            if watch:
                trace.append(self.dist[watch[0]].get(watch[1], INF))
        return t, trace

    def correct(self):
        """True iff every table equals ground-truth shortest-path cost."""
        for u in self.G:
            truth = nx.single_source_dijkstra_path_length(self.G, u, weight="w")
            truth = {d: c for d, c in truth.items() if c < INF}
            if {d: c for d, c in self.dist[u].items()} != truth:
                return False
        return True


# --------------------------------------------------------------------------
# Link state
# --------------------------------------------------------------------------
def dijkstra(adj, src):
    dist, pq = {src: 0}, [(0, src)]
    while pq:
        d, u = heapq.heappop(pq)
        if d > dist.get(u, 1e18):
            continue
        for v, w in adj.get(u, {}).items():
            OPS["spf"] += 1
            nd = d + w
            if nd < dist.get(v, 1e18):
                dist[v] = nd
                heapq.heappush(pq, (nd, v))
    return dist


class LS:
    def __init__(self, G):
        self.G = G.copy()
        self.lsdb = {u: {} for u in G}      # lsdb[u][origin] = (seq, {nbr: cost})
        self.seq = {u: 0 for u in G}
        self.msgs = 0
        self.entries = 0

    def flood(self, origins, max_ticks=500):
        inbox = {u: [] for u in self.G}
        for o in origins:
            self.seq[o] += 1
            l = (self.seq[o], {v: self.G[o][v]["w"] for v in self.G[o]})
            self.lsdb[o][o] = l
            for v in self.G[o]:
                inbox[v].append((o, l, o))
                self.msgs += 1; self.entries += len(l[1])
        t = 0
        while any(inbox.values()) and t < max_ticks:
            t += 1
            nxt = {u: [] for u in self.G}
            for u, items in inbox.items():
                for (o, l, frm) in items:
                    cur = self.lsdb[u].get(o)
                    if cur is None or l[0] > cur[0]:
                        self.lsdb[u][o] = l
                        for v in self.G[u]:
                            if v != frm:
                                nxt[v].append((o, l, u))
                                self.msgs += 1; self.entries += len(l[1])
            inbox = nxt
        return t

    def initial(self):
        return self.flood(list(self.G))

    def spf_all(self):
        """Every router runs Dijkstra on its own LSDB; returns True iff all correct."""
        ok = True
        for u in self.G:
            adj = {o: dict(l[1]) for o, l in self.lsdb[u].items()}
            d = dijkstra(adj, u)
            truth = nx.single_source_dijkstra_path_length(self.G, u, weight="w")
            ok &= (d == truth)
        return ok


def make_graph(edges):
    G = nx.Graph()
    for a, b, w in edges:
        G.add_edge(a, b, w=w)
    return G


def random_connected(n, extra=1.0, seed=0):
    r = random.Random(seed)
    G = nx.Graph()
    for i in range(1, n):
        G.add_edge(i, r.randrange(i), w=1)
    nodes = list(range(n))
    for _ in range(int(extra * n)):
        a, b = r.sample(nodes, 2)
        if not G.has_edge(a, b):
            G.add_edge(a, b, w=1)
    return G


C = {"plain": "#c0392b", "split": "#e08a1e", "poison": "#2c7fb8",
     "ls": "#1b7f5c", "grey": "#666666"}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": .25, "figure.dpi": 160})
results = {}

# --------------------------------------------------------------------------
# E1  Count-to-infinity.  (a) chain A-B-C-D, link C-D fails.
#                         (b) triangle A,B,C with stub D on A, link A-D fails.
# --------------------------------------------------------------------------
def fail_and_watch(G, mode, edge, watch):
    dv = DV(G, mode); dv.run()
    assert dv.correct()
    m0, e0 = dv.msgs, dv.entries
    dv.link_down(*edge)
    t, tr = dv.run(watch=watch)
    return dv, t, tr, dv.msgs - m0, dv.entries - e0

chain = make_graph([("A", "B", 1), ("B", "C", 1), ("C", "D", 1)])
tri = make_graph([("A", "B", 1), ("A", "C", 1), ("B", "C", 1), ("A", "D", 1)])
e1 = {"chain": {}, "triangle": {}}
traces = {"chain": {}, "triangle": {}}
for name, G, edge, watch, start in (("chain", chain, ("C", "D"), ("B", "D"), 2),
                                     ("triangle", tri, ("A", "D"), ("B", "D"), 2)):
    for mode in ("plain", "split", "poison"):
        dv, t, tr, m, e = fail_and_watch(G, mode, edge, watch)
        e1[name][mode] = {"ticks": t, "msgs": m, "entries": e,
                          "peak_metric": max(tr) if tr else start,
                          "final_correct": dv.correct()}
        traces[name][mode] = [start] + tr
    ls = LS(G); ls.initial(); ls.G.remove_edge(*edge)
    m0 = ls.msgs
    t = ls.flood(list(edge))
    e1[name]["ls"] = {"ticks": t, "msgs": ls.msgs - m0, "final_correct": ls.spf_all()}
results["E1"] = e1

fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.2), sharey=True)
for ax, name, title in ((axs[0], "chain", "(a) Chain: C-D fails"),
                        (axs[1], "triangle", "(b) Triangle + stub: A-D fails")):
    for mode, lab in (("poison", "Poison reverse"), ("plain", "Plain DV"),
                      ("split", "Split horizon (dashed)")):
        tr = traces[name][mode]
        ax.step(range(len(tr)), tr, where="post", color=C[mode], label=lab,
                lw=2.4 if mode != "split" else 1.5,
                ls="-" if mode != "split" else (0, (4, 3)),
                zorder=3 if mode == "split" else 2)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    results.setdefault("E1_identical_traces", {})[name] = (
        traces[name]["split"] == traces[name]["poison"])
    ax.axhline(INF, color=C["grey"], ls=":", lw=1)
    ax.set_title(title, fontsize=9.5)
    ax.set_xlabel("Rounds after failure")
    ax.set_ylim(0, 18)
axs[0].set_ylabel("B's metric to D")
axs[0].text(.3, INF + .4, "infinity = 16", color=C["grey"], fontsize=8)
axs[1].legend(frameon=False, loc="center right", fontsize=8)
fig.tight_layout(); fig.savefig("fig_e1_count_to_infinity.png"); plt.close(fig)

# --------------------------------------------------------------------------
# E2  Re-convergence after ONE non-bridge link failure vs network size
# --------------------------------------------------------------------------
sizes = [10, 20, 40, 80, 160]
trials = 15
keys = ["dv_ticks", "dv_msgs", "dv_entries", "dv_cpu",
        "ls_ticks", "ls_msgs", "ls_entries", "ls_cpu"]
e2 = {n: {k: [] for k in keys} for n in sizes}
for n in sizes:
    for tr in range(trials):
        G = random_connected(n, 1.0, seed=1000 * n + tr)
        br = set(map(frozenset, nx.bridges(G)))
        cand = [e for e in G.edges if frozenset(e) not in br]
        a, b = random.Random(tr).choice(cand)
        # distance vector, poison reverse + triggered updates (DV's best case)
        dv = DV(G, "poison"); dv.run(); assert dv.correct()
        m0, e0, c0 = dv.msgs, dv.entries, dv.cpu
        dv.link_down(a, b)
        t, _ = dv.run()
        assert dv.correct() or t >= 500, "DV did not converge to truth"
        e2[n]["dv_ticks"].append(t); e2[n]["dv_msgs"].append(dv.msgs - m0)
        e2[n]["dv_entries"].append(dv.entries - e0)
        e2[n]["dv_cpu"].append((dv.cpu - c0) / n)
        # link state
        ls = LS(G); ls.initial()
        ls.G.remove_edge(a, b)
        m0, e0 = ls.msgs, ls.entries
        t = ls.flood([a, b])
        OPS["spf"] = 0
        assert ls.spf_all()
        e2[n]["ls_ticks"].append(t); e2[n]["ls_msgs"].append(ls.msgs - m0)
        e2[n]["ls_entries"].append(ls.entries - e0)
        e2[n]["ls_cpu"].append(OPS["spf"] / n)
summ2 = {n: {k: round(statistics.mean(v), 1) for k, v in e2[n].items()} for n in sizes}
results["E2"] = summ2

fig, axs = plt.subplots(1, 3, figsize=(9.0, 3.0))
x = sizes
pairs = (("ticks", "Ticks to re-converge", False),
         ("entries", "Control entries sent", True),
         ("cpu", "Route computations per router", True))
for ax, (k, lab, lg) in zip(axs, pairs):
    ax.plot(x, [summ2[n]["dv_" + k] for n in x], "o-", color=C["poison"], lw=2, label="Distance vector")
    ax.plot(x, [summ2[n]["ls_" + k] for n in x], "s-", color=C["ls"], lw=2, label="Link state")
    ax.set_xlabel("Routers"); ax.set_ylabel(lab)
    if lg:
        ax.set_yscale("log")
axs[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig("fig_e2_scaling.png"); plt.close(fig)

# --------------------------------------------------------------------------
# E3  Cold start + steady-state memory per router
# --------------------------------------------------------------------------
e3 = {}
for n in sizes:
    G = random_connected(n, 1.0, seed=n)
    dv = DV(G, "poison"); t_dv, _ = dv.run(); assert dv.correct()
    ls = LS(G); t_ls = ls.initial(); assert ls.spf_all()
    e3[n] = {"edges": G.number_of_edges(),
             "dv_ticks": t_dv, "ls_ticks": t_ls,
             "dv_msgs": dv.msgs, "ls_msgs": ls.msgs,
             "dv_entries": dv.entries, "ls_entries": ls.entries,
             "dv_table": round(statistics.mean(len(dv.dist[u]) for u in G), 1),
             "dv_state": round(statistics.mean(
                 len(dv.dist[u]) + sum(len(v) for v in dv.heard[u].values()) for u in G), 1),
             "ls_state": round(statistics.mean(
                 sum(len(l[1]) for l in ls.lsdb[u].values()) for u in G), 1)}
results["E3"] = e3

# --------------------------------------------------------------------------
# E4  Replay of the 27 Oct 1980 ARPANET incident (RFC 789)
#     6-bit sequence numbers; n is LATER than m iff
#     (n > m and n - m <= 32) or (n < m and m - n > 32).
#     RFC 789 fix: replace "<=" by "<".
# --------------------------------------------------------------------------
def later(n, m, fixed=False):
    lim = (lambda d: d < 32) if fixed else (lambda d: d <= 32)
    return (n > m and lim(n - m)) or (n < m and (m - n) > 32)


def storm(fixed, ticks=80, n=14):
    """Stylised replay: every router floods an update it judges LATER, and
    answers an update it judges OLDER by sending its own copy back (the
    retransmission/echo behaviour of reliable flooding).  Metric: update
    messages processed per tick - a proxy for the CPU/buffer load that
    exhausted the IMPs."""
    G = random_connected(n, 0.8, seed=5)
    last = {u: 40 for u in G}
    inbox = {u: [] for u in G}
    for v in G[0]:                            # flaky node 0 emits the three
        inbox[v] += [(s, 0) for s in (8, 40, 44)]   # corrupted copies
    series, total = [], 0
    for _ in range(ticks):
        nxt = {u: [] for u in G}
        load = 0
        for u, items in inbox.items():
            for (s, frm) in items:
                load += 1
                if later(s, last[u], fixed):
                    last[u] = s
                    for v in G[u]:
                        nxt[v].append((s, u))
                elif later(last[u], s, fixed):
                    nxt[frm].append((last[u], u))
        inbox = nxt; total += load; series.append(load)
    return series, total

s_old, tot_old = storm(False)
s_new, tot_new = storm(True)
results["E4"] = {
    "ticks": 80,
    "load_total_original": tot_old, "load_total_fixed": tot_new,
    "load_last10_original": sum(s_old[-10:]), "load_last10_fixed": sum(s_new[-10:]),
    "load_peak_original": max(s_old), "load_peak_fixed": max(s_new),
    "ticks_until_quiet_fixed": next((i for i, v in enumerate(s_new) if v == 0), None),
    "cycle_original": {"44>40": later(44, 40), "40>8": later(40, 8), "8>44": later(8, 44)},
    "cycle_fixed": {"44>40": later(44, 40, True), "40>8": later(40, 8, True), "8>44": later(8, 44, True)}}

import math
fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.3))
pos = {8: (0.5, 0.88), 40: (0.14, 0.2), 44: (0.86, 0.2)}
for ax, fixed, title in ((axs[0], False, "Original rule: n-m <= 32"),
                         (axs[1], True, "RFC 789 fix: n-m < 32")):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.05); ax.axis("off"); ax.set_title(title, fontsize=9.5)
    for (n_, m_, lab) in ((44, 40, "44 later than 40\n(diff 4)"),
                          (40, 8, "40 later than 8\n(diff 32)"),
                          (8, 44, "8 later than 44\n(wrap: 44-8 = 36 > 32)")):
        ok = later(n_, m_, fixed)
        (x0, y0), (x1, y1) = pos[m_], pos[n_]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", lw=2.2 if ok else 1.4,
                                    color=C["plain"] if ok else C["grey"],
                                    ls="-" if ok else (0, (3, 3)),
                                    shrinkA=22, shrinkB=22,
                                    connectionstyle="arc3,rad=-0.18"))
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        off = {(44, 40): (0, -0.13), (40, 8): (-0.22, 0.06), (8, 44): (0.22, 0.06)}[(n_, m_)]
        ax.text(mx + off[0], my + off[1], lab if ok else lab.split("\n")[0].replace("later than", "NOT later than") + "\n(diff 32, not < 32)",
                ha="center", va="center", fontsize=7.2,
                color=C["plain"] if ok else C["grey"])
    for s, (x, y) in pos.items():
        ax.add_patch(plt.Circle((x, y), 0.075, fc="white", ec="#333", lw=1.4))
        ax.text(x, y, str(s), ha="center", va="center", fontsize=11, weight="bold")
fig.tight_layout(); fig.savefig("fig_e4_arpanet_storm.png"); plt.close(fig)

json.dump(results, open("results.json", "w"), indent=2)
print(json.dumps(results, indent=2))
