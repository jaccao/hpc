/*
# Copyright (C) 2017
# Luiz Fernando Jaccao <luizfernandojaccao@gmail.com>
# This file is part of HPC Project: Parallel 0-1 Knapsack Problem - UdL/FACENS Sem Fronteira.
#
# HPC Project is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HPC Project is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
*/


#include <omp.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <sys/time.h>

using namespace std;

struct Item {
    int weight;
    int value;
};

class NodeIn;
class NodeOut;

enum State { ToProcess,
             Processed };

class Node {
protected:
    Item* item;
    Item* next;
    NodeIn* nodeIn;
    NodeOut* nodeOut;
    Node* parent;
    int leftSpace;
    int totalValue;
    int level;
    State state;

public:
    Node(Node* parent, Item* item, Item* next, int level, int leftSpace, int totalValue);
    Node* process(int maxTotalValue);
    Node* createBranchIn();
    Node* createBranchOut();
    Node* getNodeOut();
    Node* getNodeIn();
    Node* createInitialBranches();
    int getTotalValue();
    int getLeftSpace();
    Item* getItem();
    State getState();
    Node* getParent();
    int getMaxPossibleValue();
};

class NodeIn : public Node {
public:
    NodeIn(Node* parent, Item* item, int level, int left_space, int total_value)
        : Node(parent, item, item + 1, level, left_space - item->weight, total_value + item->value)
    {
    }
};

class NodeOut : public Node {
public:
    NodeOut(Node* parent, Item* item, int level, int left_space, int total_value)
        : Node(parent, item, item + 1, level, left_space, total_value)
    {
        this->state = State::ToProcess;
    }
};

class NodeRoot : public Node {
public:
    NodeRoot(Item* item, int level, int left_space)
        : Node(NULL, NULL, item, level + 1, left_space, 0)
    {
    }
};

Node* Node::createInitialBranches()
{
    Node* node = this;
    while (node->createBranchIn()) {
        node = node->nodeIn;
    }
    return node;
}
int Node::getTotalValue()
{
    return this->totalValue;
}
int Node::getLeftSpace()
{
    return this->leftSpace;
}
Item* Node::getItem()
{
    return this->item;
}
State Node::getState()
{
    return this->state;
}
Node* Node::getParent()
{
    return this->parent;
}
int Node::getMaxPossibleValue()
{
    for (int i = 1; i < level - 1; i++) {
        Item* item = (this->next + i);
        if (item->weight <= this->leftSpace) {
            float profit = (float)item->value / (float)item->weight;
            return this->totalValue + profit * this->leftSpace;
        }
    }
    return this->totalValue;
}

Node::Node(Node* parent, Item* item, Item* next, int level, int left_space, int total_value)
{
    this->parent = parent;
    this->item = item;
    this->level = level - 1;
    this->nodeOut = NULL;
    this->nodeIn = NULL;
    this->leftSpace = left_space;
    this->totalValue = total_value;
    this->next = next;
    this->state = State::ToProcess;
    if (!this->level) {
        this->next = NULL;
    }
}

Node* Node::createBranchIn()
{
    if (!this->nodeIn) {
        if (this->next) {
            if (this->next->weight <= this->leftSpace) {
                this->nodeIn = new NodeIn(this, this->next, this->level, this->leftSpace, this->totalValue);
            }
        }
    }
    return this->nodeIn;
}

Node* Node::createBranchOut()
{
    if (!this->nodeOut) {
        if (this->next) {
            this->nodeOut = new NodeOut(this, this->next, this->level, this->leftSpace, this->totalValue);
        }
    }
    return this->nodeOut;
}

Node* Node::getNodeOut()
{
    return this->nodeOut;
}

Node* Node::getNodeIn()
{
    return this->nodeIn;
}

Node* Node::process(int maxTotalValue)
{
    this->state = State::Processed;
    if (this->leftSpace == 0) {
        return NULL;
    }
    int max = this->getMaxPossibleValue();
    if (max <= maxTotalValue) {
        return NULL;
    }
    if (max <= this->totalValue) {
        return NULL;
    }
    if (this->nodeIn) {
        if (max <= this->nodeIn->totalValue) {
            return NULL;
        }
    }
    if (this->createBranchOut()) {
        return this->nodeOut->createInitialBranches();
    }
    return NULL;
}

bool is_before(Item a, Item b)
{
    float profitA = (float)a.value / a.weight;
    float profitB = (float)b.value / b.weight;
    if (profitA == profitB) {
        return a.weight > b.weight;
    }
    return profitA > profitB;
}

#define NT 4

int knapsack(int max_weight, Item* items, int number_of_items)
{
    sort(items, items + number_of_items, is_before);
    NodeRoot* root = new NodeRoot(items, number_of_items, max_weight);
    Node* bestNode = root->createInitialBranches();

    int curSizeQ = 1;
    Node** curQ = new Node*[curSizeQ];
    curQ[0] = bestNode;

    int maxThread = NT;
    Node*** nextQ = new Node**[maxThread];

    while (curSizeQ) {
        int nextSizeQ[maxThread];
        Node* nextBestNode[maxThread];
        for (int t = 0; t < maxThread; t++) {
            nextSizeQ[t] = 0;
            nextQ[t] = NULL;
            nextBestNode[t] = NULL;
        }
        int perThread = (curSizeQ * 2) / maxThread + (curSizeQ * 2) % maxThread;
#pragma omp parallel num_threads(NT)
        {
            int idThread = omp_get_thread_num();
            Node* tmpBestNode = bestNode;
            int maxTotalValue = tmpBestNode->getTotalValue();
            if (nextQ[idThread]) {
                delete nextQ[idThread];
            }
            nextQ[idThread] = new Node*[perThread];
            for (int i = idThread; i < curSizeQ; i += maxThread) {
                Node* processing = curQ[i];
                if (processing->getTotalValue() > maxTotalValue) {
                    tmpBestNode = processing;
                    maxTotalValue = processing->getTotalValue();
                }
                Node* nextChild = processing->process(maxTotalValue);
                if (nextChild) {
                    nextQ[idThread][nextSizeQ[idThread]] = nextChild;
                    nextSizeQ[idThread]++;
                }
                if (processing->getParent()) {
                    if (processing->getParent()->getState() == State::ToProcess) {
                        nextQ[idThread][nextSizeQ[idThread]] = processing->getParent();
                        nextSizeQ[idThread]++;
                    }
                }
            }
            nextBestNode[idThread] = tmpBestNode;
#pragma omp barrier
#pragma omp master
            {
                delete curQ;
                curSizeQ = 0;
                for (int t = 0; t < maxThread; t++) {
                    curSizeQ += nextSizeQ[t];
                    if(nextBestNode[t]->getTotalValue() > bestNode->getTotalValue()){
                        bestNode = nextBestNode[t];
                    }
                }
                curQ = new Node*[curSizeQ];
                int countQ = 0;
                for (int s = 0; s < maxThread; s++) {
                    for (int q = 0; q < nextSizeQ[s]; q++) {
                        curQ[countQ] = nextQ[s][q];
                        countQ++;
                    }
                }
            }
        }
    }
    return bestNode->getTotalValue();
}

int main(int argc, char** argv)
{
    (void)argc;
    double tpivot1 = 0, tpivot2 = 0, tpivot3 = 0;
    struct timeval tim;

    long int number_of_items;
    long int max_load;

    FILE* test_file;
    if (!(test_file = fopen(argv[1], "r"))) {
        printf("Error opening Value file: %s\n", argv[1]);
        return 1;
    }
    fscanf(test_file, "%ld %ld\n", &number_of_items, &max_load);
    Item* items = new Item[number_of_items];

    gettimeofday(&tim, NULL);
    tpivot1 = tim.tv_sec + (tim.tv_usec / 1000000.0);

    Item temp;
    for (int cont = 0; cont < number_of_items; cont++) {
        fscanf(test_file, "%d,%d\n", &temp.value, &temp.weight);
        items[cont] = temp;
    }

    gettimeofday(&tim, NULL);
    tpivot2 = (tim.tv_sec + (tim.tv_usec / 1000000.0));
    cout << max_load << ":" << number_of_items << ":" << knapsack(max_load, items, number_of_items);
    gettimeofday(&tim, NULL);
    tpivot3 = (tim.tv_sec + (tim.tv_usec / 1000000.0));
    cout << ":" << tpivot3 - tpivot2 << ":" << tpivot3 - tpivot1 << endl;

    return 0;
}
