//
// Created by evcmo on 2/6/2022.
//

#ifndef UNTITLED3_NODE_H
#define UNTITLED3_NODE_H

#include <vector>

namespace bit {

    /**
     * The Alliances, enumerated.
     */
    enum Alliance : uint8_t
    { X, O };

    struct Node final {
        int move = -1;
        int v = 0;
        int n = 0;
        int d = 0;
        Alliance a = X;
        Node *parent = nullptr;
        std::vector<Node *> x;
    };
}
#endif //UNTITLED3_NODE_H
