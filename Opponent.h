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

    inline double UCB1
        (
        const int vi, /* The Q value of the node.                       */
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
        Node* s;
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

    inline void rolldwn
        (
        Board* const b, /* The board.                                   */
        Alliance& a,    /* The starting alliance.                       */
        int& c,         /* The crosses count.                           */
        int& o,         /* The o's count.                               */
        stack<int>& s   /* The rollout stack.                           */
        )
    {
        /*
         * Roll out a single
         * Monte-Carlo simulation.
         */
        while
            (
        ! ((c = b->hasVictory<X>())
        || (o = b->hasVictory<O>())
        || b->isFull())
            )
        {
            int m;
            do m = rand() % 9;
            while
            (b->occupiedSquare(m));
            b->mark(a, m);
            s.push(m);
            a = ~a;
        }
    }

    inline void rollup
        (
        Board* const b, /* The board.                                   */
        Alliance& a,    /* The starting alliance.                       */
        int& c,         /* The crosses count.                           */
        int& o,         /* The o's count.                               */
        stack<int>& s   /* The rollout stack.                           */
        )
    {
        /*
         * Clean up after a
         * single Monte-Carlo
         * simulation.
         */
        while (!s.empty())
        {
            b->mark
            (a = ~a, s.top());
            s.pop();
        }
    }

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
        int         c, o;
        total = 0;
        c = o = 0;
        /**
         * Don't expand
         * the current node
         * until we get 30
         * good samples.
         * Collect a sample
         * now via a single
         * Monte-Carlo
         * simulation.
         */
        if(x->n < 30)
        {
            a = ~x->a;
            rolldwn(b, a, c, o, s);
            if(c) ++winX;
            if(o) ++winO;
            ++total;
            rollup(b, a, c, o, s);
            return;
        }
        /**
         * We now have evidence
         * that this node
         * is worth expanding.
         * Expand the node
         * and run Monte-Carlo
         * simulations from
         * each of its children.
         */
        for (int i = 0; i < 9; ++i)
        {
            if (b->occupiedSquare(i))
                continue;
            l = new Node();
            l->a = ~x->a;
            l->parent = x;
            l->move = i;
            b->mark(l->a, l->move);
            a = x->a;
            rolldwn(b, a, c, o, s);
            if(c) winX += c = 1;
            if(o) winO += o = 1;
            l->n = 1;
            l->v += l->a == X? c: o;
            ++total;
            rollup(b, a, c, o, s);
            b->mark(l->a, l->move);
            x->x.push_back(l);
        }
    }

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
                winX = 1;
            if(b->hasVictory<O>())
                winO = 1;
            total = b->isFull()?
                0: 1;
            if(winX || winO || !total)
            {
                /**
                 * Only set the total
                 * simulations
                 * to one if there is
                 * a victory. Ties
                 * don't count.
                 */
                total = winX | winO;
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
                 * Rollout
                 * Monte-Carlo
                 * simulation(s)
                 * according to
                 * a random
                 * Default Policy.
                 * Expand the
                 * current node
                 * if it is
                 * promising.
                 */
                rollout
                (
                b, winX, winO, total, x
                );
                break;
            }
            /**
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
         * Back-up to the root
         * and update all nodes
         * encountered.
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
        x->v += winX;
    }

    inline void treeWalk(Node* n, int depth)
    {
        for(Node* x: n->x)
        {
            for(int i = 0; i < depth; ++i)
            { std::cout << '\t'; }
            std::cout << x->move << ": " << x->v << '/' << x->n << '\n';
            treeWalk(x, depth + 1);
        }
    }

    inline void destroyTree(Node* n)
    {
        for(Node* const x: n->x)
        { destroyTree(x); delete x; }
    }

    inline int search(Board * const b) {
        clock_t const time = clock();
        Node n; n.a = O; n.d = 0; //n.x.clear();
        //if(!n.x.empty()) std::cout << "woah";
        do simulate(b, &n);
        while((clock() - time) < 10000);
        treeWalk(&n, 0);
        int move = selectNode(&n)->move;
        destroyTree(&n);
        return move;
    }
}

#endif //BITTACTOE_OPPONENT_H
