# Tiny-MCTS
##### <i>A smol version of tic-tac-toe MCTS with little-to-no optimization.</i>
Monte Carlo Tree Search is a versatile search algorithm for making decisions. 
Its applications range from game-playing to vehicle-piloting.
## Algorithm
MCTS is an iterative algorithm that repeatedly executes four subroutines.

```
1. Selection
2. Expansion
3. Simulation
4. Back-propagation
```

### Selection
Starting at the root node, MCTS navigates the search tree in memory until reaching a node on the fringe.
At each internal node, a child is chosen according to some "tree policy."
A common policy is UCT, which treats each node like a multi-armed-bandit, exploring less and exploiting more as nodes accrue visits.

### Expansion
### Simulation
### Back-propagation



