//Name : Sahil Punchhi
//zID  : z5204256
//Course: COMP 9024 Assignment 4 - Graphs
//Date: May 01, 2019

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>

// A vertex is a 2D point
typedef struct Vertex {
    int x; // x-coordinate
    int y; // y-coordinate
} Vertex;

// each edge is a pair of vertices (end-points)
typedef struct Edge {
    Vertex *p1; // first end point
    Vertex *p2; // second end point
} Edge;

// A vertex node stores a vertex and other information, and you need to expand this type
typedef struct VertexNode {
    Vertex *v;
    struct VertexNode *next;
    double weight;
} VertexNode;

typedef struct GraphRep *Graph;

typedef struct GraphRep { // graph header
    VertexNode **vertices; // an array of vertices or a linked list of vertices
    int nV; // #vertices
    int nE; // #edges
} GraphRep;

// the time complexity of CreateEmptyGraph() is O(k) where k is a constant
Graph CreateEmptyGraph()
{

    Graph g = malloc(sizeof(GraphRep));
    assert(g != NULL);
    // initially number of vertices is 0
    g->nV = 0;
    // initially number of edges is 0
    g->nE = 0;
    g->vertices = malloc(100 * sizeof(VertexNode*));
    return g;
    
}

// create new node
VertexNode *createNewNode(Vertex *k) {
    VertexNode *new = malloc(sizeof(VertexNode));
    assert(new != NULL);
    new->v = k;
    new->next = NULL;
    new->weight = 0;
    return new;
}

// to find if node exists in linked list
bool inLL(VertexNode *N, Vertex *k) {
    if (N == NULL)
        return false;
    if (N->v->x == k->x && N->v->y == k->y)
        return true;
    
    return inLL(N->next, k);
}

// insert node in a linked list
VertexNode *insertLL(VertexNode *N, Vertex *k) {
    if (inLL(N, k))
        return N;
    
    // add new node at the beginning
    VertexNode *new = createNewNode(k);
    new->weight = sqrt(abs(  ((k->x - N->v->x) * (k->x - N->v->x))  +  ((k->y - N->v->y) * (k->y - N->v->y))  ));
    VertexNode *head = N;
    
    while(head->next != NULL){
        head = head->next;
    }
    head->next = new;
    return N;
}

// delete node in a linked list
VertexNode * deleteLL(VertexNode * N, Vertex *k) {
    if (N == NULL)
        return N;
    if (N->v->x == k->x && N->v->y == k->y)
        return N->next;
    
    N->next = deleteLL(N->next, k);
    return N;
    
}

// free a linked list
void freeLL(VertexNode * N) {
    if (N != NULL) {
        freeLL(N->next);
        free(N);
    }
}

// assign rank of each vertex, starting from 0
int RankVertex(Graph g, Vertex *k){
    int i;
    
    for(i = 0; i < g->nV; i++){
        if(g->vertices[i]->v->x == k->x && g->vertices[i]->v->y == k->y)
            return i;
    }
    return -5;
}

// the time complexity of InsertEdge() is O(n) where n is the number of vertices in the graph
int InsertEdge(Graph g, Edge *e)
{
    Vertex *Vertex_1=malloc(sizeof(Vertex)); Vertex_1=e->p1;
    Vertex *Vertex_2=malloc(sizeof(Vertex)); Vertex_2=e->p2;
    
    if(RankVertex(g,Vertex_1)==-5 && RankVertex(g,Vertex_2)==-5){
        
        // we insert both the vertices
        g->vertices[g->nV] = createNewNode(Vertex_1);
        g->vertices[g->nV]->weight = 0;
        g->nV = g->nV + 1;
        g->vertices[g->nV] = createNewNode(Vertex_2);
        g->vertices[g->nV]->weight = 0;
        g->nV = g->nV + 1;
        g->vertices[RankVertex(g,Vertex_1)]=insertLL(g->vertices[RankVertex(g,Vertex_1)],Vertex_2);
        g->vertices[RankVertex(g,Vertex_2)]=insertLL(g->vertices[RankVertex(g,Vertex_2)],Vertex_1);
        g->nE=g->nE+1; return 1;}
    
    else if(RankVertex(g,Vertex_1)!=-5 && RankVertex(g,Vertex_2)==-5){
        
        // we insert second vertex, since first vertex is already there
        g->vertices[g->nV] = createNewNode(Vertex_2);
        g->vertices[g->nV]->weight = 0;
        g->nV = g->nV + 1;
        g->vertices[RankVertex(g,Vertex_1)]=insertLL(g->vertices[RankVertex(g,Vertex_1)],Vertex_2);
        g->vertices[RankVertex(g,Vertex_2)]=insertLL(g->vertices[RankVertex(g,Vertex_2)],Vertex_1);
        g->nE=g->nE+1; return 1;}
    
    else if(RankVertex(g,Vertex_1)==-5 && RankVertex(g,Vertex_2)!=-5){
        
        // we insert first vertex, since second vertex is already there
        g->vertices[g->nV] = createNewNode(Vertex_1);
        g->vertices[g->nV]->weight = 0;
        g->nV = g->nV + 1;
        g->vertices[RankVertex(g,Vertex_1)]=insertLL(g->vertices[RankVertex(g,Vertex_1)],Vertex_2);
        g->vertices[RankVertex(g,Vertex_2)]=insertLL(g->vertices[RankVertex(g,Vertex_2)],Vertex_1);
        g->nE=g->nE+1; return 1;}
    
    else if(RankVertex(g,Vertex_1)!=-5 && RankVertex(g,Vertex_2)!=-5){
        if(!inLL(g->vertices[RankVertex(g,Vertex_1)],Vertex_2)){
            
            // both the vertices are already there, only new edge is inserted in linked list
            g->vertices[RankVertex(g,Vertex_1)]=insertLL(g->vertices[RankVertex(g,Vertex_1)],Vertex_2);
            g->vertices[RankVertex(g,Vertex_2)]=insertLL(g->vertices[RankVertex(g,Vertex_2)],Vertex_1);
            g->nE=g->nE+1; return 1;}
        
        else{
            // edge already exists
            return 0;}
    }
    return 0;
}

// the time complexity of DeleteEdge() is O(n) where n is the number of vertices in the graph
void DeleteEdge(Graph g, Edge *e)
{
    Vertex *Vertex_1=malloc(sizeof(Vertex)); Vertex_1=e->p1;
    Vertex *Vertex_2=malloc(sizeof(Vertex)); Vertex_2=e->p2;
    
    // we check if edge exists
    if(inLL(g->vertices[RankVertex(g,Vertex_1)],Vertex_2)){
        g->vertices[RankVertex(g,Vertex_1)]=deleteLL(g->vertices[RankVertex(g,Vertex_1)],Vertex_2);
        g->vertices[RankVertex(g,Vertex_2)]=deleteLL(g->vertices[RankVertex(g,Vertex_2)],Vertex_1);
        g->nE = g->nE - 1;
        
        if(g->vertices[RankVertex(g,Vertex_1)]->next == NULL){
            int p1 = RankVertex(g,Vertex_1);
            g->vertices[RankVertex(g,Vertex_1)] = deleteLL(g->vertices[RankVertex(g,Vertex_1)],Vertex_1);
            for(int m1 = p1; m1 < g->nV - 1; m1++){
                g->vertices[m1] = g->vertices[m1+1];}
            g->vertices[g->nV-1] = NULL;
            free(g->vertices[g->nV-1]);
            g->nV = g->nV -1;
        }
        
        if(g->vertices[RankVertex(g,Vertex_2)]->next == NULL){
            int p2 = RankVertex(g,Vertex_2);
            g->vertices[RankVertex(g,Vertex_2)] = deleteLL(g->vertices[RankVertex(g,Vertex_2)],Vertex_2);
            for(int m2 = p2; m2 < g->nV - 1; m2++){
                g->vertices[m2] = g->vertices[m2+1];}
            g->vertices[g->nV-1] = NULL;
            free(g->vertices[g->nV-1]);
            g->nV = g->nV-1;
        }
    }
}

// to check if 2 vertices form an edge
bool adjacent(Graph g, Vertex *Vertex_1, Vertex *Vertex_2) {
    assert(g != NULL);
    return inLL(g->vertices[RankVertex(g,Vertex_1)],Vertex_2);
}

#define MAX_NODES 200
int visited[MAX_NODES];
int visited2[MAX_NODES];

// Depth First Search to find connected vertices, recursion is used
bool dfsPathCheck(Graph g, int nV, Vertex* v, Vertex* dest) {
    
    for (int t = 0; t < nV; t++){
        if (adjacent(g, v, g->vertices[t]->v) && visited2[t] == -1) {
            visited2[t] = 0;
            if (g->vertices[t]->v->x == dest->x && g->vertices[t]->v->y == dest->y)
                return true;
            else if (dfsPathCheck(g, nV, g->vertices[t]->v, dest))
                return true;
        }
    }
    return false;
}

bool findPathDFS(Graph g, int nV, Vertex* src, Vertex* dest) {
    for (int p = 0; p < nV; p++)
        visited2[p] = -1;
    visited2[RankVertex(g,src)] = 0;
    return dfsPathCheck(g, nV, src, dest);
}

// the time complexity of ReachableVertices() is O(m+n) where n is the number of vertices in the graph and m is the number of edges in the graph
// Recursive DFS is used, for which time complexity is O(m+n)
void ReachableVertices(Graph g, Vertex *v)
{   int s= RankVertex(g,v);
    printf("The output of ReachableVertices() for (%d,%d) is ",v->x, v->y);
    int m = 0;
    for(int f=0; f < g->nV; f++){
        if(f != s){
            
            // if the path between 2 vertices exists, means they are connected
            if(findPathDFS(g, g->nV, v, g->vertices[f]->v)){
                m = m+1;
                if(m==1){
                    printf("(%d,%d)",g->vertices[f]->v->x, g->vertices[f]->v->y);}
                else if(m != 1){
                    printf(",(%d,%d)",g->vertices[f]->v->x, g->vertices[f]->v->y);}
            }
        }
    }
    printf(".\n\n");
}

// min heap node
typedef struct MinHeapNode
{
    Vertex* v;
    int distance;
} MinHeapNode;

// min heap
typedef struct MinHeap
{
    int size;   // number of heap nodes
    int volume; // volume of min heap
    int *location; // for decrease_Key()
    MinHeapNode **list;
} MinHeap;

// creating a new minimum Heap Node
MinHeapNode* newMinHeapNode(Vertex* v, int distance)
{
    MinHeapNode* minHeapNode = (MinHeapNode*) malloc(sizeof(MinHeapNode));
    minHeapNode->v = v;
    minHeapNode->distance = distance;
    return minHeapNode;
}

// creating minimum heap
MinHeap* createMinHeap(int volume)
{
    MinHeap* minHeap = (MinHeap*) malloc(sizeof(MinHeap));
    minHeap->location = (int *)malloc(volume * sizeof(int));
    minHeap->size = 0; minHeap->volume = volume;
    minHeap->list =
    (MinHeapNode**) malloc(volume * sizeof(MinHeapNode*));
    return minHeap;
}

void swap_node_minHeap( MinHeapNode** m, MinHeapNode** n)
{
    MinHeapNode* s = *m;
    *m = *n;
    *n = s;
}

void toHeapify(MinHeap* minHeap, int index, Graph g)
{
    int least, left, right;
    least = index; right = 2 * index + 2; left = 2 * index + 1;
    
    if (right < minHeap->size &&
        minHeap->list[right]->distance < minHeap->list[least]->distance )
        least = right;
    
    if (left < minHeap->size &&
        minHeap->list[left]->distance < minHeap->list[least]->distance )
        least = left;
    
    if (least != index)
    {
        MinHeapNode *least_node = minHeap->list[least];
        MinHeapNode *index_node = minHeap->list[index];
        minHeap->location[RankVertex(g,least_node->v)] = index;
        minHeap->location[RankVertex(g,index_node->v)] = least;
        
        // Swapping nodes of min heap
        swap_node_minHeap(&minHeap->list[least], &minHeap->list[index]);
        
        toHeapify(minHeap, least, g);
    }
}

// to find if heap is empty
int isEmpty(MinHeap* minHeap)
{
    return minHeap->size == 0;
}

// we get minimum node from the heap
MinHeapNode* get_Min(MinHeap* minHeap, Graph g)
{
    if (isEmpty(minHeap))
        return NULL;
    
    MinHeapNode* first_node = minHeap->list[0];
    
    // replacing first_node node with last node
    MinHeapNode* lastNode = minHeap->list[minHeap->size - 1]; minHeap->list[0] = lastNode;
    
    // updating location of last node
    minHeap->location[RankVertex(g,first_node->v)] = minHeap->size-1;
    minHeap->location[RankVertex(g,lastNode->v)] = 0;
    
    // heapifying the first_node
    minHeap->size = minHeap->size - 1; toHeapify(minHeap, 0, g);
    return first_node;
}

// to decrease distance value of vertex v
void decrease_Key(MinHeap* minHeap, Vertex* v, int distance, Graph g)
{
    int s = minHeap->location[RankVertex(g,v)];
    
    // updating distance value of the node
    minHeap->list[s]->distance = distance;
    
    // time complexity of this is O(Log n) where n is the number of vertices in graph
    while (s && minHeap->list[s]->distance < minHeap->list[(s - 1) / 2]->distance)
    {
        // Swapping node with its parent
        minHeap->location[RankVertex(g,minHeap->list[s]->v)] = (s-1)/2;
        minHeap->location[RankVertex(g,minHeap->list[(s-1)/2]->v)] = s;
        swap_node_minHeap(&minHeap->list[s], &minHeap->list[(s - 1) / 2]);
        
        // moving on to the parent index
        s = (s - 1) / 2;
    }
}

// we check if vertex is in heap
bool is_In_MinHeap(MinHeap *minHeap, Vertex* v, Graph g)
{
    if (minHeap->location[RankVertex(g,v)] < minHeap->size)
        return true;
    return false;
}

// the time complexity of ShortestPath() is O((m+n)*log n) where n is the number of vertices in the graph and m is the number of edges in the graph
void ShortestPath(Graph g, Vertex* u, Vertex *v)
{
    if(RankVertex(g,u) == -5 || RankVertex(g,v) == -5){
        return;
    }
    if(!findPathDFS(g, g->nV, u, v)){
        return;
    }
    
    // swapping source with destination to get desired output
    // distance does not change on swapping end points
    Vertex* src1 = v;
    Vertex* dest1 = u;
    int V = g->nV;
    int visited3[V];
    int distance[V];
    
    for (int p = 0; p < g->nV; p++)
        visited3[p] = -1;
    visited3[RankVertex(g,src1)] = RankVertex(g,src1);
    
    // creating minimum heap
    MinHeap* minHeap = createMinHeap(V);
    
    // initially setting distance to 1000000
    for (int k = 0; k < V; k++)
    {
        distance[k] = 1000000;
        minHeap->list[k] = newMinHeapNode(g->vertices[k]->v, distance[k]);
        minHeap->location[k] = k;
    }
    
    minHeap->list[RankVertex(g,src1)] = newMinHeapNode(src1, distance[RankVertex(g,src1)]);
    minHeap->location[RankVertex(g,src1)] = RankVertex(g,src1);
    distance[RankVertex(g,src1)] = 0;
    decrease_Key(minHeap, src1, distance[RankVertex(g,src1)],g);
    
    minHeap->size = V;
    
    // min heap has nodes for which shortest distance is not yet finalized
    while (!isEmpty(minHeap))
    {
        // we get vertex with minimum distance
        MinHeapNode* minHeapNode = get_Min(minHeap,g);
        int u = RankVertex(g,minHeapNode->v);
        
        // we traverse through all adjacent vertices of u and update their distance values
        VertexNode* connected_vertex = g->vertices[u]->next  ;
        while (connected_vertex != NULL)
        {
            Vertex* v = connected_vertex->v  ;
            
            // if distance to v from u is less than its distance calculated before and shortest distance for v is not finalized
            if (is_In_MinHeap(minHeap, v,g) && distance[u] != 1000000 &&
                connected_vertex->weight + distance[u] < distance[RankVertex(g,v)])
            {
                distance[RankVertex(g,v)] = distance[u] + connected_vertex->weight;
                // to record from which adjacent vertex we have finalized this distance
                visited3[RankVertex(g,v)] = u;
                
                // updating distance value in minimum heap
                decrease_Key(minHeap, v, distance[RankVertex(g,v)], g);
            }
            connected_vertex = connected_vertex->next;
        }
    }
    
    // we find ranks of source and destination vertex
    int y = RankVertex(g,src1);
    int z = RankVertex(g,dest1);
    
    printf("The shortest path between (%d,%d) and (%d,%d) is ",u->x,u->y,v->x,v->y);
    int m = z;
    while (m != y) {
        printf("(%d,%d),", g->vertices[m]->v->x,g->vertices[m]->v->y);
        m = visited3[m];
    }
    printf("(%d,%d).\n\n", src1->x,src1->y);
}

// the time complexity of FreeGraph() is O(n*m) where n is the number of vertices in the graph and m is the number of edges in the graph
void FreeGraph(Graph g)
{
    assert(g != NULL);
    int i;
    
    // freeing linked list
    for (i = 0; i < g->nV; i++)
        freeLL(g->vertices[i]);
    
    // freeing graph
    free(g);
}

// implementing queues for breadth first search
typedef struct QueueRep *queue;

typedef struct QueueRep {
    int   length;
    VertexNode *head;
    VertexNode *tail;
} QueueRep;

// setting up empty queue
queue newQueue() {
    queue Q = malloc(sizeof(QueueRep));
    Q->length = 0;
    Q->head = NULL;
    Q->tail = NULL;
    return Q;
}

// removing unwanted queue
void dropQueue(queue Q) {
    VertexNode *curr = Q->head;
    while (curr != NULL) {
        VertexNode *temp = curr->next;
        free(curr);
        curr = temp;
    }
    free(Q);
}

// to check if queue is empty
int QueueIsEmpty(queue Q) {
    return (Q->length == 0);
}

// insert an int at end of queue
void QueueEnqueue(queue Q, Vertex* v) {
    VertexNode *new = malloc(sizeof(VertexNode));
    assert(new != NULL);
    new->v = v;
    new->next = NULL;
    if (Q->tail != NULL) {
        Q->tail->next = new;
        Q->tail = new;
    } else {
        Q->head = new;
        Q->tail = new;
    }
    Q->length++;
}

// remove int from front of queue
Vertex* QueueDequeue(queue Q) {
    assert(Q->length > 0);
    VertexNode *p = Q->head;
    Q->head = Q->head->next;
    if (Q->head == NULL) {
        Q->tail = NULL;
    }
    Q->length--;
    Vertex *d = p->v;
    free(p);
    return d;
}

// the time complexity of ShowGraph() is O(m+n) where n is the number of vertices in the graph and m is the number of edges in the graph
// here BFS algorithm is used for which time complexity is O(m+n)
void ShowGraph(Graph g)
{   printf("The output of ShowGraph() is");
    int k=g->nE;
    Vertex* z;
    
    // initially all the vertices are marked as unvisited, -1
    for (int h = 0; h < g->nV; h++){
        visited[h]= -1;
    }
    
    // to start with, first vertex is marked as visited, 0
    // Breadth First Search algorithm
    visited[0] = 0;
    queue Q = newQueue();
    QueueEnqueue(Q, g->vertices[0]->v);
    while (! QueueIsEmpty(Q)) {
        z = QueueDequeue(Q);
        visited[RankVertex(g,z)] = 0;
        
        for (int t = 0; t < g->nV; t++){
            if (adjacent(g, z, g->vertices[t]->v) && visited[t] == -1) {
                printf(" (%d,%d),(%d,%d)",z->x,z->y,g->vertices[t]->v->x,g->vertices[t]->v->y);
                k = k-1;
                if (k == 0){
                    printf(".\n\n");
                    return;}
                else{
                    QueueEnqueue(Q, g->vertices[t]->v);}
            }
        }
        
        // if number of edges is not equal to 0 but queue is empty means there are more connected components
        // we pick first unvisited vertex from all the vertices
        if(k != 0 && QueueIsEmpty(Q)){
            for(int f=0; f<g->nV; f++){
                if(visited[f]== -1){
                    QueueEnqueue(Q, g->vertices[f]->v);
                    break;
                }
            }
        }
    }
    return;
}

int main() //sample main for testing
{
    Graph g1;
    Edge *e_ptr;
    Vertex *v1, *v2;
    
    // Create an empty graph g1;
    g1=CreateEmptyGraph();
    
    // Create first connected component
    // Insert edge (0,0)-(0,10)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=0;
    v2->x=0;
    v2->y=10;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Insert edge (0,0)-(5,6)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=0;
    v2->x=5;
    v2->y=6;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Insert edge (0, 10)-(10, 10)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=10;
    v2->x=10;
    v2->y=10;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Insert edge (0,10)-(5,6)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=10;
    v2->x=5;
    v2->y=6;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Insert edge (0,0)-(5,4)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=0;
    v2->x=5;
    v2->y=4;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Insert edge (5, 4)-(10, 4)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=5;
    v1->y=4;
    v2->x=10;
    v2->y=4;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Insert edge (5,6)-(10,6)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=5;
    v1->y=6;
    v2->x=10;
    v2->y=6;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Insert edge (10,10)-(10,6)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=10;
    v1->y=10;
    v2->x=10;
    v2->y=6;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Insert edge (10, 6)-(10, 4)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=10;
    v1->y=6;
    v2->x=10;
    v2->y=4;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Create second connected component
    // Insert edge (20,4)-(20,10)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=20;
    v1->y=4;
    v2->x=20;
    v2->y=10;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Insert edge (20,10)-(30,10)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=20;
    v1->y=10;
    v2->x=30;
    v2->y=10;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    // Insert edge (25,5)-(30,10)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=25;
    v1->y=5;
    v2->x=30;
    v2->y=10;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    if (InsertEdge(g1, e_ptr)==0) printf("edge exists\n");
    
    //Display graph g1
    ShowGraph(g1);
    
    // Find the shortest path between (0,0) and (10,6)
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=0;
    v2->x=10;
    v2->y=6;
    ShortestPath(g1, v1, v2);
    free(v1);
    free(v2);
    
    // Delete edge (0,0)-(5, 6)
    e_ptr = (Edge*) malloc(sizeof(Edge));
    assert(e_ptr != NULL);
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=0;
    v2->x=5;
    v2->y=6;
    e_ptr->p1=v1;
    e_ptr->p2=v2;
    DeleteEdge(g1, e_ptr);
    free(e_ptr);
    free(v1);
    free(v2);
    
    // Display graph g1
    ShowGraph(g1);
    
    // Find the shortest path between (0,0) and (10,6)
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=0;
    v2->x=10;
    v2->y=6;
    ShortestPath(g1, v1, v2);
    free(v1);
    free(v2);
    
    // Find the shortest path between (0,0) and (25,5)
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v2=(Vertex *) malloc(sizeof(Vertex));
    assert(v2 != NULL);
    v1->x=0;
    v1->y=0;
    v2->x=25;
    v2->y=5;
    ShortestPath(g1, v1, v2);
    free(v1);
    free(v2);
    
    // Find reachable vertices of (0,0)
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v1->x=0;
    v1->y=0;
    ReachableVertices(g1, v1);
    free(v1);
    
    // Find reachable vertices of (20,4)
    v1=(Vertex*) malloc(sizeof(Vertex));
    assert(v1 != NULL);
    v1->x=20;
    v1->y=4;
    ReachableVertices(g1, v1);
    free(v1);
    
    // Free graph g1
    FreeGraph(g1);
    
    return 0;
}


