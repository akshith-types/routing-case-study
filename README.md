# Routing case study: distance-vector vs link-state

Code and results for the case study **"Who Wins the Comparison?"** (Computer Networks, BCSE308L).

A failure-first comparison of distance-vector and link-state routing, using a Python simulator of both protocols. Every converged state is checked against independently computed shortest paths.

## Files

| File | What it does |
|---|---|
| `routesim.py` | Experiments E1-E4: count-to-infinity, re-convergence vs network size, cold start, and the 1980 ARPANET sequence-number rule. Writes `results.json` and the `fig_e*.png` charts. |
| `worked.py` | Solves the 5-router worked example (Dijkstra and distance-vector). Writes `worked.json`. |
| `fig_topology.py` | Draws Figure 1 (the 5-router network). Needs `worked.json`, so run `worked.py` first. |
| `results.json`, `worked.json` | Output of a reference run, so results can be seen without running anything. |
| `fig_*.png` | The figures used in the report. |

## Run

```
pip3 install networkx matplotlib
python3 worked.py
python3 routesim.py
python3 fig_topology.py
```

Random seeds are fixed, so results are repeatable.

## Reference

E. C. Rosen, "Vulnerabilities of Network Control Protocols: An Example", RFC 789, 1981.
