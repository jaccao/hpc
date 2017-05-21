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
#include <cmath>
#include <climits>

using namespace std;

int DEBUG_INFO=0;

struct Item {
    int weight;
    int value;
};

class KnapsackResult{
public:
    int maxTotalValue;
    int leftSpace;
    int center;
    int level;
    int depth;
    int processLevel;
    static const int messageSize=7;
    int buffer[7];
    KnapsackResult(){
        this->maxTotalValue = 0;
        this->leftSpace = 0;
        this->center = INT_MIN;
        this->level = 0;
        this->depth = 0;
        this->processLevel = 0;
    }
    KnapsackResult(int maxTotalValue,int leftSpace,int center, int level, int depth, int processLevel){
        this->maxTotalValue=maxTotalValue;
        this->leftSpace=leftSpace;
        this->center=center;
        this->level = level;
        this->depth = depth;
        this->processLevel = processLevel;
    }
    KnapsackResult(const KnapsackResult &o){
        this->leftSpace=o.leftSpace;
        this->center=o.center;
        this->maxTotalValue=o.maxTotalValue;
        this->level=o.level;
        this->depth = o.depth;
        this->processLevel = o.processLevel;
    }
    void printRes(const char *str){
        cout << str << endl;
        cout << "Max Value: " << this->maxTotalValue << endl;
        cout << "Left Space: " << this->leftSpace << endl;
        cout << "Center: " << this->center << endl;
        cout << "Level: " << this->level << endl;
        cout << "Depth: " << this->depth << endl;
        cout << "Process level: " << this->processLevel << endl;
        cout << endl;
    }
    void operator =(KnapsackResult o){
        this->leftSpace=o.leftSpace;
        this->center=o.center;
        this->maxTotalValue=o.maxTotalValue;
        this->level=o.level;
        this->depth = o.depth;
        this->processLevel = o.processLevel;
    }
    int* toMessage(int command){
        this->buffer[0]=command;
        this->buffer[1]=this->maxTotalValue;
        this->buffer[2]=this->leftSpace;
        this->buffer[3]=this->center;
        this->buffer[4]=this->level;
        this->buffer[5]=this->depth;
        this->buffer[6]=this->processLevel;
        return this->buffer;
    }
    int fromMessage(){
        int* message=this->buffer;
        this->maxTotalValue=message[1];
        this->leftSpace=message[2];
        this->center=message[3];
        this->level=message[4];
        this->depth=message[5];
        this->processLevel=message[6];
        return message[0];
    }
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
    Node* process(KnapsackResult b);
    Node* createBranchIn();
    Node* createBranchOut();
    Node* getNodeOut();
    Node* getNodeIn();
    Node* createInitialBranches(KnapsackResult b);
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
    static int getMaxItems(int maxItems = -1);
    KnapsackResult getResult();
    void printNode(const char *str);
};

class NodeIn : public Node {
public:
    NodeIn(Node* parent, Item* item, int level, int left_space, int total_value, int center)
        : Node(parent, item, item + 1, level, left_space - item->weight, total_value + item->value, center)
    {
        this->type=1;
    }
};

class NodeOut : public Node {
public:
    NodeOut(Node* parent, Item* item, int level, int left_space, int total_value, int center)
        : Node(parent, item, item + 1, level, left_space, total_value, center - 1)
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
        Node::getMaxItems(level+1);
    }
};

Node* Node::createInitialBranches(KnapsackResult b)
{
    Node* node = this;
    while (node->createBranchIn()) {
        node = node->nodeIn;
        if(node->getLevel() < b.level){
            return NULL;
        }
        if(node->getCenter() < b.center){
            return NULL;
        }
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

int Node::getMaxItems(int maxItems)
{
    static int staticMaxItems=0;
    if(maxItems > -1){
        staticMaxItems = maxItems;
    }
    return staticMaxItems;
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

Node* Node::process(KnapsackResult b)
{
    if (this->leftSpace == 0) {
        return NULL;
    }
    int max = this->getMaxPossibleValue();
    if (max <= b.maxTotalValue) {
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
        return this->nodeOut->createInitialBranches(b);
    }
    return NULL;
}

void Node::printNode(const char *str){
    cout << str << endl;
    cout << "Level: " << this->getLevel() << endl;
    cout << "Leftspace: " << this->getLeftSpace() << endl;
    cout << "Center: " << this->getCenter() << endl;
    cout << "Value: " << this->getTotalValue() << endl;
    cout << "Max Value: " << this->getMaxPossibleValue() << endl;
    cout << endl;
}

int Node::getLevel(){
    return this->level;
}

KnapsackResult Node::getResult(){
    return KnapsackResult(this->getTotalValue(),this->getLeftSpace(),this->getCenter(),this->getLevel(),0,0);
}

int Node::getCenter(){
    return this->center;
}

bool is_before(Item a, Item b)
{
    float profitA = (float)a.value / a.weight;
    float profitB = (float)b.value / b.weight;
    if (profitA == profitB) {
        return a.weight < b.weight;
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

KnapsackResult knapsackNode(Node* starter,KnapsackResult b){
    Node* process = starter;
    KnapsackResult r = b;
    int depth=starter->getLevel();
    while(process){
        Node* node = process->process(r);
        if(node){
            process = node;
        } else {
            if(process != starter){
                process->destroyBranchOut();
                process = process->getParent();
            }else{
                process = NULL;
            }
        }
        if(process){
            if(depth > process->getLevel()){
                depth = process->getLevel();
            }
            if(r.maxTotalValue < process->getTotalValue()){
                r = process->getResult();
            }
        }
        r.depth++;
    }
    return r;
}

KnapsackResult knapsackLevel(NodeRoot* root,KnapsackResult b){
    Node* node=root;
    while (node->getLevel()!=b.processLevel) {
        node=node->getNodeIn();
    }
    KnapsackResult t = knapsackNode(node, b);
    node->destroyBranchOut();
    return t;
}


int knapsack(NodeRoot* root)
{
    KnapsackResult b;
    Node* currentNode = root->createInitialBranches(b);
    b.maxTotalValue = currentNode->getTotalValue();
    int maxRange=5;
    double maxDepth=root->getLevel()*maxRange*pow(10,-6);
    int rootMaxLevel = currentNode->getLevel()+maxRange;
    if(DEBUG_INFO) currentNode->printNode("INITIAL");
    bool stopping=false;
    while (!stopping) {
        b.processLevel=currentNode->getLevel();
        KnapsackResult t=knapsackLevel(root,b);
        if(t.maxTotalValue > b.maxTotalValue){
            b.maxTotalValue = t.maxTotalValue;
            b.leftSpace = t.leftSpace;
            rootMaxLevel = currentNode->getLevel()+maxRange;
            if(DEBUG_INFO) currentNode->printNode("ROOT BEST");
            if(DEBUG_INFO) t.printRes("new result");
            if(DEBUG_INFO) b.printRes("new best");
        }
        double c=(double)t.depth/pow(2,23);
        maxDepth-=c;
        if(currentNode->getLevel() > rootMaxLevel && maxDepth < 0.0){
            stopping = true;
        }else{
            currentNode = currentNode->getParent();
            if(currentNode){
                if(DEBUG_INFO) currentNode->printNode("ROOT PROCESSING");
            }else{
                stopping=true;
            }
        }
    }
    return b.maxTotalValue;
}

int knapsackSend(NodeRoot* root, ParallelProcess* parallelProcess)
{
    KnapsackResult b;
    Node* currentNode = root->createInitialBranches(b);
    b.maxTotalValue = currentNode->getTotalValue();
    int maxRange=5;
    double maxDepth=root->getLevel()*maxRange*pow(10,-6);
    int rootMaxLevel = currentNode->getLevel()+maxRange;
    if(DEBUG_INFO) currentNode->printNode("INITIAL");

    int currentNumberOfProcess=0;
    MPI_Status status;
    bool stopping=false;
    while (!stopping) {
        KnapsackResult t;
        MPI_Recv((void *)t.buffer,t.messageSize,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
        int command = t.fromMessage();
        if(DEBUG_INFO) cout << "0: receiving from " << status.MPI_SOURCE << " command " << command << endl;
        if(command == 1){
            parallelProcess[status.MPI_SOURCE].processing=false;
            currentNumberOfProcess--;
            if(DEBUG_INFO) cout << "0: processing " << t.maxTotalValue << " : " << t.leftSpace << " : " << t.center << endl;
            if(t.maxTotalValue > b.maxTotalValue){
                b.maxTotalValue = t.maxTotalValue;
                b.leftSpace = t.leftSpace;
                rootMaxLevel = currentNode->getLevel()+maxRange;
                if(DEBUG_INFO) currentNode->printNode("ROOT BEST");
                if(DEBUG_INFO) t.printRes("new result");
                if(DEBUG_INFO) b.printRes("new best");
            }
            double c=(double)t.depth/pow(2,23);
            maxDepth-=c;
            if(currentNode->getLevel() > rootMaxLevel && maxDepth < 0.0){
                stopping = true;
            }else{
                currentNode = currentNode->getParent();
                if(currentNode){
                    if(DEBUG_INFO) currentNode->printNode("ROOT PROCESSING");
                }else{
                    stopping=true;
                }
            }
        }
        if(stopping){
            if(DEBUG_INFO)  cout << "0: sending to " << status.MPI_SOURCE << " command 0" << endl;
            MPI_Send((void*)b.toMessage(0),b.messageSize,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
        } else {
            b.processLevel = currentNode->getLevel();
            parallelProcess[status.MPI_SOURCE].processing=true;
            currentNumberOfProcess++;
            if(DEBUG_INFO) b.printRes("SENDING");
            MPI_Send((void*)b.toMessage(1),b.messageSize,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
        }
    }
    return b.maxTotalValue;
}

void knapsackReceive(NodeRoot* root,int rank){
    KnapsackResult b;
    root->createInitialBranches(b);
    MPI_Status status;

    MPI_Send((void*)b.toMessage(0),b.messageSize,MPI_INT,0,0,MPI_COMM_WORLD);

    MPI_Recv((void *)b.buffer,b.messageSize,MPI_INT,0,0,MPI_COMM_WORLD,&status);
    int command = b.fromMessage();
    if(DEBUG_INFO)  cout << rank <<": receiving from " << status.MPI_SOURCE << " command " << command << endl;
    while(command == 1){
        if(DEBUG_INFO) b.printRes("PROCESSING");
        KnapsackResult t=knapsackLevel(root,b);
        if(DEBUG_INFO) t.printRes("PROCESSED");
        MPI_Send((void*)t.toMessage(1),t.messageSize,MPI_INT,0,0,MPI_COMM_WORLD);
        MPI_Recv((void *)b.buffer,b.messageSize,MPI_INT,0,0,MPI_COMM_WORLD,&status);
        command = b.fromMessage();
        if(DEBUG_INFO) cout << rank <<": receiving from 0 command " << command << endl;
    }
}

int mainMPI(NodeRoot* root, int max_load, int number_of_items)
{
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

    if(world_rank == 0){
        gettimeofday(&tim, NULL);
        double tpivot2 = (tim.tv_sec + (tim.tv_usec / 1000000.0));

        ParallelProcess parallelProcess[world_size];

        cout << max_load << ":" << number_of_items << ":" << knapsackSend(root,parallelProcess);
        gettimeofday(&tim, NULL);
        double tpivot3 = (tim.tv_sec + (tim.tv_usec / 1000000.0));
        cout << ":" << tpivot3 - tpivot2 << endl;

        // stoping the nodes processing
        MPI_Status status;
        int countParallel=0;
        for(int p=1; p<world_size; p++){
            if(parallelProcess[p].processing){
                parallelProcess[p].processing=false;
                countParallel++;
            }
        }
        KnapsackResult t;
        while(countParallel){
            MPI_Recv((void *)t.buffer,t.messageSize,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
            if(DEBUG_INFO) cout << "0: stoping: " << status.MPI_SOURCE << endl;
            MPI_Send((void*)t.toMessage(0),t.messageSize,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
            countParallel--;
        }
    }else{
        knapsackReceive(root,world_rank);
    }

    MPI_Finalize();

    return 0;
}

int mainSerial(NodeRoot* root, int max_load, int number_of_items)
{
    struct timeval tim;

    gettimeofday(&tim, NULL);
    double tpivot2 = (tim.tv_sec + (tim.tv_usec / 1000000.0));

    cout << max_load << ":" << number_of_items << ":" << knapsack(root);
    gettimeofday(&tim, NULL);
    double tpivot3 = (tim.tv_sec + (tim.tv_usec / 1000000.0));
    cout << ":" << tpivot3 - tpivot2 << endl;

    return 0;
}

int main(int argc, char** argv)
{
    if(argc<2){
        return -1;
    }

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

#ifdef SERIAL
    mainSerial(root,max_load,number_of_items);
#else
    mainMPI(root,max_load,number_of_items);
#endif

    return 0;
}
