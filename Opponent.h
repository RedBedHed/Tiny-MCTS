//
// Created by evcmo on 10/27/2021.
//
#pragma once
#ifndef BITTACTOE_OPPONENT_H
#define BITTACTOE_OPPONENT_H
#include "Board.h"
#include <cmath>
#include <cfloat>
#include <stack>

namespace opponent {

    using bit::Board;
    using bit::BoardLength;
    using bit::X;
    using bit::O;
    using bit::Node;
    using bit::Alliance;
    using std::vector;
    using std::stack;

    template<bool>
    void select
    (Board*, Node*);
    template<bool>
    void expand
    (Board*, double&, double&, double&, Node*);
    void old_simulate
    (Board*, double&, double&, double&, Alliance);
    void new_simulate
    (Board*, double&, double&, double&, Alliance, Node*);
    void back_propagate
    (Board*, double&, double&, double&, Node*, Node*);
    Node* selectNode
    (Node*);
    void randMove
    (Board*, Alliance&);


    constexpr double UCB1
        (
        const double vi, /* The value of the node.                         */
        const double n,  /* The number of visits to the parent             */
        const double ni  /* The number of visits to the node               */
        )
    {
        return 
            (ni <= 0)? DBL_MAX :
            (vi / ni) +
            1.42 * std::sqrt
            (std::log(n) / ni);
    }


    template<bool INIT>
    inline void select
        (
        Board * const b,/* The board.                                   */
        Node  * const n /* The root node.                               */
        )
    {
        Node*  x     = n;
        double winX  = 0,
               winO  = 0,
               total = 0;

        /**
         * If initializing the
         * root node, expand
         * in the rollout step
         * and update the
         * root.
         */
        if constexpr (INIT)
        {
            expand<INIT>
            (
            b, winX, winO, total, x
            );
            x->n += total;
            x->v += n->v += 
                ((int)~n->a) * winX + 
                ((int) n->a) * winO;
            return;
        }

        /**
         * Navigate the MC
         * Tree.
         */
        for(;;)
        {

            /*
             * Check if the
             * current node
             * is a winning
             * node. If it
             * is, we can
             * nudge the
             * search in
             * its direction.
             */
            if(x->a == X) 
            {
                if(b->hasVictory<X>())
                { winX = total = 1.0; break; }
            }
            else 
                if(b->hasVictory<O>())
                { winO = total = 1.0; break; }
            if(b->isFull())
            { 
                winX = winO = 0.5; 
                total = 1.0; break; 
            }

            /*
             * If we reach a
             * node without
             * children...
             */
            if(x->x.empty())
            {
                /**
                 * Rollout lines
                 * of play according
                 * to a random Default
                 * Policy. Expand
                 * the current node
                 * if it is promising.
                 */
                expand<INIT>
                (b, winX, winO, total, x);
                break;
            }

            /**
             * (1) Selection.
             * Select a node
             * according to
             * the UC Tree
             * policy.
             */
            x = selectNode(x);

            /**
             * Make the move.
             */
            b->mark(x->a, x->move);
        }

        /**
         * Update all nodes on the
         * path back to the root.
         */
        back_propagate
        (b, winX, winO, total, n, x);
    }


    template<bool INIT>
    inline void expand
        (
        Board* const b, /* The board.                                   */
        double& winX,   /* The win count of x.                          */
        double& winO,   /* The win count of o.                          */
        double& total,  /* The total number of simulations.             */
        Node* const x   /* The leaf tree node selected by tree policy   */
        )
    {
        /*
         * If this is an
         * initialization
         * rollout, expand
         * the node.
         */
        if constexpr (!INIT) {

            /**
             * Don't expand
             * the current node
             * until we get 30
             * good samples.
             * Collect a sample
             * now via a single
             * rollout.
             */
            if (x->n < 30)
            {
                old_simulate
                (
                    b, winX, winO, 
                    total, ~x->a
                );
                return;
            }
        }

        /**
         * (2) Expansion.
         * We now have
         * evidence that this
         * node is worth
         * expanding. Expand
         * the node and rollout
         * random simulations
         * from each of its
         * children.
         */
        uint16_t
        bb = b->legalMoves();
        for (; bb; bb &= bb - 1)
        {
            Node* n = new Node
            (
                8 - bit::bitScanFwd(bb), 
                ~x->a, x
            );
            b->mark(n->a, n->move);
            new_simulate
            (   
                b, winX, winO, 
                total, x->a, n
            );
            b->mark(n->a, n->move);
            x->x.push_back(n);
        }
    }


    inline void old_simulate
        (
        Board* const bx,  /* The board.                                 */
        double& winX,     /* The win count of x.                        */
        double& winO,     /* The win count of o.                        */
        double& total,    /* The total number of simulations.           */
        const Alliance ax /* The starting alliance                      */
        )
    {
        Board    b = *bx;
        Alliance a = ax;
        double   d = 0.0;

        for(;; a = ~a)
        {
            if(a == X) 
            {
                if(b.hasVictory<X>())
                { 
                    winX = 1.0 - d;
                    total = 1.0;
                    break; 
                }
            } else
                if(b.hasVictory<O>())
                { 
                    winO = 1.0 - d; 
                    total = 1.0;
                    break; 
                }
            if(b.isFull())
            {
                winX = winO = 0.5;
                total = 1.0;
                break;
            }
            d += 0.01;
            randMove(&b, a);
        }
    }


    inline void new_simulate
        (
        Board* const bx,  /* The board.                                 */
        double& winX,     /* The win count of x.                        */
        double& winO,     /* The win count of o.                        */
        double& total,    /* The total number of simulations.           */
        const Alliance ax,/* The starting alliance                      */
        Node* const n
        )
    {
        /**
         * (3) Simulation.
         * Rollout a random
         * line of play to 
         * the end of the
         * game.
         */
        Board    b = *bx;
        Alliance a = ax;
        double   d = 0.0;

        for(a = ax ;; a = ~a)
        {
            if(n->a == X) 
            {
                if(b.hasVictory<X>())
                {
                    winX += n->v = 1.0 - d;
                    total += n->n = 1.0;
                    break;
                }
            } else 
                if(b.hasVictory<O>())
                {
                    winO += n->v = 1.0 - d;
                    total += n->n = 1.0;
                    break;
                }
            if(b.isFull())
            {
                winO += 0.5;
                winX += n->v = 0.5;
                total += n->n = 1.0;
                break;
            }
            d += 0.01;
            randMove(&b, a);
        }
    }


    inline void back_propagate
        (
        Board* const b, /* The board.                                   */
        double& winX,   /* The win count of x.                          */
        double& winO,   /* The win count of o.                          */
        double& total,  /* The total number of simulations.             */
        Node* const n,  /* The root node                                */
        Node* x         /* The current node                             */
        )
    {
        /**
         * (4) Back-propagation.
         * return to the root
         * and update all nodes
         * on the path.
         */
        while(x != n)
        {
            x->n += total;
            x->v += 
                ((int)~x->a) * winX + 
                ((int) x->a) * winO;
            b->mark(x->a, x->move);
            x = x->parent;
        }
        x->n += total;
        x->v += winO;
    }


    inline Node* selectNode
        (
        Node* const n  /* The parent node.                              */
        )
    {
        /**
         * Select a
         * descendant
         * of the node
         * n according to
         * the Upper-
         * Confidence bound
         * (applied to trees)
         * Policy.
         */
        Node* s = nullptr;
        double m = -1;
        for(Node* const x: n->x)
        {
            double h = UCB1
            (x->v, n->n, x->n);
            if(h > m)
            { m = h; s = x; }
        }
        return s;
    }


    inline void randMove
        (
        Board* const b, /* The board.                                   */
        Alliance& a     /* The starting alliance.                       */
        )
    {
        /*
         * Make a random move. 
         * This is inefficient. 
         * Do better.
         */
        int m;
        do m = rand() % 9;
        while
        (b->occupiedSquare(m));
        b->mark(a, m);
    }


    inline Node* child(Board* b, Node* n, int i) 
    {
        if(n->x.empty())
            select<true>(b, n);
        for(Node* x: n->x)
            if (i == x->move)
                return x;
        return n;
    } 


    inline int treeWalk(Node* n, int depth)
    {
        if(n->x.empty())
            return 1;
        int c = 0;
        for(Node* x: n->x)
        {
            for(int i = 0; i < depth; ++i)
            { std::cout << '\t'; }
            std::cout << x->move << ": " << x->v << '/' << x->n << '\n';
            c += treeWalk(x, depth + 1);
        }
        return c + 1;
    }


    inline void destroyTree(Node* n)
    {
        for(Node* const x: n->x)
        { destroyTree(x); delete x; }
    }


    inline int search(Board * const b, Node*& n, bool f) 
    {
        clock_t const time = clock();
        int i = f? 100000: 10000;
        do select<false>(b, n);
        while((clock() - time) < i);
        int l = treeWalk(n, 0);
        std::cout << "Node Count:" << l << '\n';
        return (n = selectNode(n))->move;
    }
}

#endif //BITTACTOE_OPPONENT_H
