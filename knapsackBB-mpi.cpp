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

#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <sys/time.h>

using namespace std;

#define DEBUG_INFO 0

struct Item {
    int weight;
    int value;
};

class NodeIn;
class NodeOut;

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
    int center;
    int type;
public:
    Node(Node* parent, Item* item, Item* next, int level, int leftSpace, int totalValue, int center);
    ~Node();
    Node* process(int maxTotalValue);
    Node* createBranchIn();
    Node* createBranchOut();
    Node* getNodeOut();
    Node* getNodeIn();
    Node* createInitialBranches();
    int getTotalValue();
    int getLeftSpace();
    Item* getItem();
    Node* getParent();
    int getMaxPossibleValue();
    void destroyBranchOut();
    void destroyBranchIn();
    int getLevel();
    int getCenter();
    int getType();
};

class NodeIn : public Node {
public:
    NodeIn(Node* parent, Item* item, int level, int left_space, int total_value, int center)
        : Node(parent, item, item + 1, level, left_space - item->weight, total_value + item->value, center -1)
    {
        this->type=1;
    }
};

class NodeOut : public Node {
public:
    NodeOut(Node* parent, Item* item, int level, int left_space, int total_value, int center)
        : Node(parent, item, item + 1, level, left_space, total_value, center + 1)
    {
        this->type=-1;
    }
};

class NodeRoot : public Node {
public:
    NodeRoot(Item* item, int level, int left_space)
        : Node(NULL, NULL, item, level + 1, left_space, 0, 0)
    {
        this->type=0;
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

int Node::getType()
{
    return this->type;
}

Item* Node::getItem()
{
    return this->item;
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

Node::Node(Node* parent, Item* item, Item* next, int level, int left_space, int total_value, int center)
{
    this->parent = parent;
    this->item = item;
    this->level = level - 1;
    this->nodeOut = NULL;
    this->nodeIn = NULL;
    this->leftSpace = left_space;
    this->totalValue = total_value;
    this->next = next;
    this->center = center;
    if (!this->level) {
        this->next = NULL;
    }
}

Node::~Node()
{
    this->destroyBranchIn();
    this->destroyBranchOut();
}

void Node::destroyBranchOut()
{
    if(this->nodeOut){
        delete this->nodeOut;
        this->nodeOut = NULL;
    }
}

void Node::destroyBranchIn()
{
    if(this->nodeIn){
        delete this->nodeIn;
        this->nodeIn = NULL;
    }
}

Node* Node::createBranchIn()
{
    if (!this->nodeIn) {
        if (this->next) {
            if (this->next->weight <= this->leftSpace) {
                this->nodeIn = new NodeIn(this, this->next, this->level, this->leftSpace, this->totalValue, this->center);
            }
        }
    }
    return this->nodeIn;
}

Node* Node::createBranchOut()
{
    if (!this->nodeOut) {
        if (this->next) {
            this->nodeOut = new NodeOut(this, this->next, this->level, this->leftSpace, this->totalValue, this->center);
        }
    } else {
        return NULL;
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

int Node::getLevel(){
    return this->level;
}

int Node::getCenter(){
    return this->center;
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

class ParallelProcess {
public:
    bool processing;
    ParallelProcess(){
        this->processing=false;
    }
};

class KnapsackResult{
public:
    int maxTotalValue;
    int leftSpace;
    int center;
    KnapsackResult(int maxTotalValue,int leftSpace,int center){
        this->maxTotalValue=maxTotalValue;
        this->leftSpace=leftSpace;
        this->center=center;
    }
    void operator =(KnapsackResult o){
        this->leftSpace=o.leftSpace;
        this->center=o.center;
        this->maxTotalValue=o.maxTotalValue;
    }
};

KnapsackResult knapsackNode(Node* starter,int maxTotalValue){
    Node* process = starter;
    int leftSpace = -1;
    int center = 0;
    while(process){
        Node* node = process->process(maxTotalValue);
        if(node){
            process = node;
        } else {
            if(process != starter){
                process->destroyBranchOut();
            }
            process = process->getParent();
        }
        if(process){
            if(maxTotalValue < process->getTotalValue()){
                maxTotalValue = process->getTotalValue();
            }
            if(process->getLeftSpace() == 0){
                leftSpace = process->getLeftSpace();
                center = process->getCenter();
            }
        }
        if(process == starter){
            process = NULL;
        }
    }
    return KnapsackResult(maxTotalValue,leftSpace,center);
}

KnapsackResult knapsackLevel(NodeRoot* root,int level,int maxTotalValue){
    Node* node=root;

    while (node->getLevel()!=level) {
        node=node->getNodeIn();
    }
    int leftSpace=-1;
    int center=0;
    KnapsackResult t = knapsackNode(node, maxTotalValue);
    if(t.maxTotalValue > maxTotalValue){
        maxTotalValue = t.maxTotalValue;
    }
    if(t.leftSpace==0){
        leftSpace=0;
        center=t.center;
    }
    node->destroyBranchOut();
    return KnapsackResult(maxTotalValue,leftSpace,center);
}

int knapsack(NodeRoot* root, ParallelProcess* parallelProcess)
{
    Node* currentNode = root->createInitialBranches();
    int maxTotalValue = currentNode->getTotalValue();
    int rightLevel=root->getLevel();
    float factor=((float)currentNode->getLevel() / (float)root->getLevel())*0.5033;
    int minRight = root->getLevel()*(factor*2.0-1.0);
    int currentNumberOfProcess=0;
    if(currentNode->getLeftSpace()==0){
        return maxTotalValue;
    }

    const int recSize=4;
    int recMessage[recSize];
    const int sendSize=3;
    int sendMessage[sendSize];

    MPI_Status status;

    bool stopping=false;

    // while (!stopping || currentNumberOfProcess) {
    while (!stopping) {
        MPI_Recv((void *)recMessage,recSize,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);

        int command=recMessage[0];
        if(DEBUG_INFO) cout << "0: receiving from " << status.MPI_SOURCE << " command " << command << endl;

        if(command==1){
            parallelProcess[status.MPI_SOURCE].processing=false;
            currentNumberOfProcess--;
            KnapsackResult t = KnapsackResult(recMessage[1],recMessage[2],recMessage[3]);
            if(DEBUG_INFO) cout << "0: processing " << t.maxTotalValue << " : " << t.leftSpace << " : " << t.center << endl;
            // KnapsackResult t=knapsackLevel(root,currentNode->getLevel(),maxTotalValue);
            if(t.maxTotalValue > maxTotalValue){
                maxTotalValue = t.maxTotalValue;
            }
            if(t.leftSpace == 0 && rightLevel == root->getLevel()){
                rightLevel = t.center;
                if(rightLevel<minRight){
                    rightLevel=minRight;
                }
            }
            if(currentNode->getCenter()>rightLevel){
                stopping=true;
            }
        }

        if(stopping){
            sendMessage[0]=0; // exit
            if(DEBUG_INFO)  cout << "0: sending to " << status.MPI_SOURCE << " command 0" << endl;
            MPI_Send((void*)sendMessage,sendSize,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
        } else {
            sendMessage[0]=1; // process
            sendMessage[1]=currentNode->getLevel();
            sendMessage[2]=maxTotalValue;
            parallelProcess[status.MPI_SOURCE].processing=true;
            currentNumberOfProcess++;
            if(DEBUG_INFO)  cout << "0: sending to " << status.MPI_SOURCE << " command 1 : " << sendMessage[1] << " : " << sendMessage[2] << endl;
            MPI_Send((void*)sendMessage,sendSize,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
            if(currentNode == root){
                stopping = true;
            }else{
                currentNode = currentNode->getParent();
            }
        }
    }
    return maxTotalValue;
}

void knapsackReceive(NodeRoot* root,int rank){
    root->createInitialBranches();
    const int recSize=3;
    int recMessage[recSize];
    const int sendSize=4;
    int sendMessage[sendSize];
    MPI_Status status;

    sendMessage[0]=0; // request data only

    MPI_Send((void*)sendMessage,sendSize,MPI_INT,0,0,MPI_COMM_WORLD);

    MPI_Recv((void *)recMessage,recSize,MPI_INT,0,0,MPI_COMM_WORLD,&status);
    int command = recMessage[0];
    if(DEBUG_INFO)  cout << rank <<": receiving from " << status.MPI_SOURCE << " command " << command << endl;
    while(command==1){
        int level=recMessage[1];
        int maxTotalValue=recMessage[2];
        if(DEBUG_INFO) cout << rank <<": processing : " << level << " : "<< maxTotalValue << endl;
        KnapsackResult t=knapsackLevel(root,level,maxTotalValue);
        if(DEBUG_INFO) cout << rank <<": processed : " << level << " : "<< maxTotalValue << endl;
        sendMessage[0]=1; // with result
        sendMessage[1]=t.maxTotalValue;
        sendMessage[2]=t.leftSpace;
        sendMessage[3]=t.center;
        if(DEBUG_INFO) cout << rank <<": sending to 0 command 1 : " << sendMessage[1] << " : "<< sendMessage[2] << " : "<< sendMessage[3] << endl;
        MPI_Send((void*)sendMessage,sendSize,MPI_INT,0,0,MPI_COMM_WORLD);
        MPI_Recv((void *)recMessage,recSize,MPI_INT,0,0,MPI_COMM_WORLD,&status);
        command = recMessage[0];
        if(DEBUG_INFO) cout << rank <<": receiving from 0 command " << command << endl;
    }
}

int main(int argc, char** argv)
{
    if(argc<2){
        return -1;
    }
    struct timeval tim;

    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size < 2) {
        fprintf(stderr,"Requires at least two processes.\n");
        exit(-1);
    }

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    long int number_of_items;
    long int max_load;

    FILE* test_file;
    if (!(test_file = fopen(argv[1], "r"))) {
        printf("Error opening Value file: %s\n", argv[1]);
        return 1;
    }
    fscanf(test_file, "%ld %ld\n", &number_of_items, &max_load);
    Item* items = new Item[number_of_items];

    Item temp;
    for (int cont = 0; cont < number_of_items; cont++) {
        fscanf(test_file, "%d,%d\n", &temp.value, &temp.weight);
        items[cont] = temp;
    }

    sort(items, items + number_of_items, is_before);

    NodeRoot* root = new NodeRoot(items, number_of_items, max_load);

    if(world_rank == 0){
        gettimeofday(&tim, NULL);
        double tpivot2 = (tim.tv_sec + (tim.tv_usec / 1000000.0));

        ParallelProcess parallelProcess[world_size];

        cout << max_load << ":" << number_of_items << ":" << knapsack(root,parallelProcess);
        gettimeofday(&tim, NULL);
        double tpivot3 = (tim.tv_sec + (tim.tv_usec / 1000000.0));
        cout << ":" << tpivot3 - tpivot2 << endl;

        // stoping the nodes processing

        const int recSize=4;
        int recMessage[recSize];
        const int sendSize=3;
        int sendMessage[sendSize];
        MPI_Status status;

        int countParallel=0;

        for(int p=1; p<world_size; p++){
            if(parallelProcess[p].processing){
                parallelProcess[p].processing=false;
                countParallel++;
            }
        }
        while(countParallel){
            sendMessage[0]=0;
            MPI_Recv((void *)recMessage,recSize,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
            if(DEBUG_INFO) cout << "0: stoping: " << status.MPI_SOURCE << endl;
            MPI_Send((void*)sendMessage,sendSize,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
            countParallel--;
        }

    }else{
        knapsackReceive(root,world_rank);
    }

    MPI_Finalize();

    return 0;
}
