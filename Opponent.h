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

    constexpr double UCB1
        (
        const int vi, /* The MC value of the node.                      */
        const int n,  /* The number of simulations run under the parent */
        const int ni  /* The number of simulations run under the node   */
        )
    {
        if(ni <= 0)
            return DBL_MAX;
        return
            ((double) vi /
             (double) ni) +
            1.42 * std::sqrt
            (
            (double)
            std::log(n) /
            (double) ni
            );
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

    inline void rollup
        (
        Board* const b, /* The board.                                   */
        Alliance& a,    /* The starting alliance.                       */
        stack<int>& s   /* The rollout stack.                           */
        )
    {
        /*
         * Clean up after a
         * single rollout.
         */
        while (!s.empty())
        {
            b->mark
            (a = ~a, s.top());
            s.pop();
        }
    }

    template<bool INIT>
    inline void rollout
        (
        Board* const b, /* The board.                                   */
        int& winX,      /* The win count of x.                          */
        int& winO,      /* The win count of o.                          */
        int& total,     /* The total number of simulations.             */
        Node* const x   /* The leaf tree node selected by tree policy   */
        )
    {
        Alliance    a;
        Node*       l;
        stack<int>  s;
        total = 0;
        /*
         * If this is an
         * initialization
         * rollout, expand
         * the node.
         */
        if constexpr (INIT)
            goto expand;
        /**
         * (3) Simulation.
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
            for(a = ~x->a ;; a = ~a)
            {
                if(b->hasVictory<X>())
                {
                    winX += 2;
                    break;
                }
                if(b->hasVictory<O>())
                {
                    winO += 2;
                    break;
                }
                if(b->isFull())
                {
                    ++winX;
                    ++winO;
                    break;
                }
                ++total;
                int m;
                do m = rand() % 9;
                while
                (b->occupiedSquare(m));
                b->mark(a, m);
                s.push(m);
            }
            rollup(b, a, s);
            return;
        }
        /**
         * (2) Expansion.
         * (3) Simulation.
         * We now have
         * evidence that this
         * node is worth
         * expanding. Expand
         * the node and rollout
         * random simulations
         * from each of its
         * children.
         */
        expand:
        uint16_t
        bb = b->legalMoves();
        for (; bb; bb &= bb - 1)
        {
            l = new Node();
            l->a = ~x->a;
            l->parent = x;
            l->move =
            8 - bit::bitScanFwd(bb);
            b->mark(l->a, l->move);
            for(a = x->a ;; a = ~a)
            {
                if(b->hasVictory<X>())
                {
                    winX += 2;
                    if(l->a == X)
                        l->v = 2;
                    break;
                }
                if(b->hasVictory<O>())
                {
                    winO += 2;
                    if(l->a == O)
                        l->v = 2;
                    break;
                }
                if(b->isFull())
                {
                    ++winO;
                    ++winX;
                    l->v = 1;
                    break;
                }
                total += l->n = 1;
                int m;
                do m = rand() % 9;
                while
                (b->occupiedSquare(m));
                b->mark(a, m);
                s.push(m);
            }
            rollup(b, a, s);
            b->mark(l->a, l->move);
            x->x.push_back(l);
        }
    }

    template<bool INIT>
    inline void simulate
        (
        Board * const b,/* The board.                                   */
        Node  * const n /* The root node.                               */
        )
    {
        Node* x   = n;
        int winX  = 0,
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
            rollout<INIT>
            (
            b, winX, winO, total, x
            );
            x->n += total;
            x->v += winO;
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
            if(b->hasVictory<X>())
            {
                winX = 2;
                total = 1;
                break;
            }
            if(b->hasVictory<O>())
            {
                winO = 2;
                total = 1;
                break;
            }
            if(b->isFull()) {
                ++winX;
                ++winO;
                total = 1;
                break;
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
                rollout<INIT>
                (
                b, winX, winO, total, x
                );
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
         * (4) Back-propagation.
         * return to the root
         * and update all nodes
         * on the path.
         */
        while(x != n)
        {
            x->n += total;
            x->v += x->a == X?
                winX: winO;
            b->mark(x->a, x->move);
            x = x->parent;
        }
        x->n += total;
        x->v += winO;
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

    inline int search(Board * const b) {
        clock_t const time = clock();
        Node n; n.a = O; n.d = 0;
        simulate<true>(b, &n);
        do simulate<false>(b, &n);
        while((clock() - time) < 100000);
        int i = treeWalk(&n, 0);
        std::cout << "Node Count:" << i << '\n';
        int move = selectNode(&n)->move;
        destroyTree(&n);
        return move;
    }
}

#endif //BITTACTOE_OPPONENT_H
