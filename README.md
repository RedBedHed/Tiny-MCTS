# Tiny-MCTS
##### <i>A smol version of tic-tac-toe MCTS with little-to-no optimization.</i>
Monte Carlo Tree Search is a versatile search algorithm for making decisions. 
Its applications range from game-playing to vehicle-piloting.
# Algorithm Description
MCTS is an iterative algorithm that executes four subroutines until time runs out.

```
1. Selection
2. Expansion
3. Simulation
4. Back-propagation
```

## Selection
Starting at the root node, MCTS navigates the search tree in memory until reaching a node on the fringe.
At each internal node, a child is chosen according to some "tree policy." The tree policy should take into 
account the probability of winning stored in each child. A common policy is UCT, which treats each node 
like a multi-armed-bandit, exploring less and exploiting more as nodes accrue visits.

```c++
Node x = root;
while(true) {
  if(terminal node x) {
    value, num = evaluate();
    break;
  }
  if(leaf node x) {
    value, num = expand(x);
    break;
  }
  x = select_child(x);
  do_action(x.a);
}
back_propagate(x, prob, num);
```

## Expansion
Once a leaf node has been selected, MCTS expands it into its children. A new node is added to the tree for
each legal action according to the leaf's state.

```c++
score, count = 0;
for_each(action a) {
  do_action(a);
  prob = simulate();
  score += prob
  count++;
  undo_action(a);
  x.add_node(a, prob);
}
return score, count;
```

## Simulation
MCTS runs a simulation beneath each new node. The simulation is often dependent on the domain. A generic
approach is to run a playout with random legal actions. In zero-sum two-player games, the result of the
playout is usually the probability of winning for either player.
- 0   means that a player loses.
- 0.5 means that a player draws.
- 1   means that a player wins.

```c++
while(true) {
  if(terminal state) {
    prob = evaluate();
    break;
  }
  a = random_action();
  do_action(a);
  stack.push(a);
}
while(not stack.empty()) {
  undo_action(stack.pop());
}
return prob;
```

## Back-propagation
Once probabilities have been determined, MCTS propagates them up the tree, following the unique path from 
the leaf to the root. In zero-sum two-player games, the winning probability is added to the value of each
node owned by the relevant alliance. The simulation count is added to each node on the path, regardless of 
alliance.

```c++
while(x != root) {
  update(x, value, num);
  undo_action(x.a);
  x = x.parent;
}
update(root, value, num);
```

