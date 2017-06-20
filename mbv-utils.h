#include <limits>
#include <iostream>
#include <forward_list>
#include <utility>
#include <algorithm>
#include <chrono>

using namespace std;

struct UnionOperation{
	int parent;
	int child;
	bool increasedParentRank;

	UnionOperation(): parent(0), child(0), increasedParentRank(false){}
	UnionOperation(int p, int c, bool i): parent(p), child(c), increasedParentRank(i){}
};

class UnionFind{
private:
	int* parent;
	int* rank;
	int elementCount;
	int disjointSetCount;

	int findSet (int x) {
		int top = x;

		while(parent[top] != top) {
			top = parent[top];
		}

		return top;
	}

public:
	UnionFind() : parent(nullptr), elementCount(0), disjointSetCount(0) {	}

	UnionFind(int n) : elementCount(n), disjointSetCount(n) {
		parent = new int[elementCount];
		rank = new int[elementCount];
		for (int i = 0; i < elementCount; ++i)
		{
			parent[i] = i;
			rank[i] = 0;
		}
		
	}

	~UnionFind() {
		if(parent){
			delete [] parent;
		}
		if (rank)
		{
			delete [] rank;
		}
	}

	void printElements() {
		for (int i = 0; i < elementCount; ++i)
		{
			cout<<parent[i]<<",";
		}
		cout<<endl;
	}

	bool unionSet (int x, int y, UnionOperation &op) {
		int setX = findSet(x);
		int setY = findSet(y);
		if(setX == setY) {
			return false;
		}
		
		if(rank[setX] > parent[setY]){
			parent[setY] = setX;
			op = UnionOperation(setX,setY,false);
		}
		else{
			parent[setX] = setY;
			rank[setY]++;
			op = UnionOperation(setY,setX,true);
		}
		
		disjointSetCount--;
		return true;
	}

	void undoUnion(UnionOperation &op) {
		parent[op.child] = op.child;
		if(op.increasedParentRank){
			rank[op.parent]--;
		}
		disjointSetCount++;
		return;
	}

	int getElementCount() {
		return elementCount;
	}

	int getElement(int e) {
		return parent[e];
	}

	int getDisjointSetCount() {
		printElements();
		return disjointSetCount;
	}
};

class Graph {
public:
	//For every vertex, a forward_list of edges it's connected to
	forward_list<int>* vertices;
	//A list with all the edges
	pair<int, int>* edges;
	//A list with vertex degrees
	int* vertexDegrees;
	//Sizes
	int nVertices, mEdges;
	//Amount of edges already added
	int addedEdges;

	Graph() : vertices(nullptr), edges(nullptr), vertexDegrees(nullptr), nVertices(0), mEdges(0), addedEdges(0) {}

	Graph(int n, int m) : nVertices(n), mEdges(m), addedEdges(0) {
		vertices = new forward_list<int>[nVertices];
		vertexDegrees = new int[nVertices];
		edges = new pair<int, int>[mEdges];

		for (int i = 0; i < nVertices; ++i)
		{
			vertexDegrees[i] = 0;
		}
	} 

	void addEdge(int u, int v) {
		if(addedEdges == mEdges) {
			//cout<<"Already added all edges"<<endl;
		}
		edges[addedEdges] = pair<int, int>(u, v);
		vertices[u].push_front(addedEdges);
		vertexDegrees[u]++;
		vertices[v].push_front(addedEdges);
		vertexDegrees[v]++;
		addedEdges++;
	}

};

class sort_indices
{
   private:
     float* valueArr;
   public:
     sort_indices(float* arr) : valueArr(arr) {}
     bool operator()(int i, int j) const { return valueArr[i]<valueArr[j]; }
};

class MBVSolutionUndo {
private:

public:
	int branchVertexCount, activeEdgeCount;
	UnionFind* vertexConnected;
	int* edgeState;
	float* edgeWeight;
	int* edgeIndex;
	int* vertexDegrees;
	int* vertexProhibitedDegrees;
	UnionOperation* unionsList;
	int unionCounter;
	int* unionEdges;
	Graph *G;
    float* relaxedSolution;
    int currentRelaxedSolution;
    float constantRelaxedSubtractor;

	MBVSolutionUndo(Graph &sourceGraph) : 
	G(&sourceGraph), branchVertexCount(0), activeEdgeCount(0), unionCounter(0), currentRelaxedSolution(0), constantRelaxedSubtractor(0.0) {
		vertexDegrees = new int[G->nVertices];
		vertexProhibitedDegrees = new int[G->nVertices];
		vertexConnected = new UnionFind(G->nVertices);
		edgeState = new int[G->mEdges];
		edgeWeight = new float[G->mEdges];
		edgeIndex = new int[G->mEdges];
		unionsList = new UnionOperation[G->nVertices-1];
		unionEdges = new int[G->nVertices-1];
        relaxedSolution = new float[G->nVertices];

		int deg = 0;

		for (int v = 0; v < G->nVertices; ++v)
		{
			vertexDegrees[v] = 0;
			vertexProhibitedDegrees[v] = 0;
            relaxedSolution[v] = 0.0;
		}

		for (int e = 0; e < G->mEdges; ++e)
		{
			edgeState[e] = 0;
			edgeWeight[e] = 0.0;
			edgeIndex[e] = e;
		}

		//sort
		sortEdges();

		//cout<<"Sorted!"<<endl;
	}

	MBVSolutionUndo(int bvc) : G(nullptr), branchVertexCount(bvc), activeEdgeCount(0) {}

	~MBVSolutionUndo() {
		/**/
		//if (vertexConnected)
		//{
		//	delete vertexConnected;
		//}
		if (edgeState)
		{
			delete [] edgeState;
		}
		if (vertexDegrees)
		{
			delete [] vertexDegrees;
		}
		if (edgeWeight)
		{
			delete [] edgeWeight;
		}
		if (edgeIndex)
		{
			delete [] edgeIndex;
		}
		/**/
	}

	void sortEdges() {
		int deg = 0;
		int v = 0;
		//Calculate edge weights
        for (int e = 0; e < G->mEdges; ++e)
		{
			edgeIndex[e] = e;
            edgeWeight[e] = 0.0;

			//Prepare "heuristic" weights for sorting
			v = G->edges[e].first;
			deg = G->vertexDegrees[v] - vertexProhibitedDegrees[v];
			if(deg > 2) {
				edgeWeight[e] += 1.0/deg;
			}

			v = G->edges[e].second;
			deg = G->vertexDegrees[v] - vertexProhibitedDegrees[v];
			if(deg > 2) {
				edgeWeight[e] += 1.0/deg;
			}
		}
        //Calculate constant weight subtractor
        constantRelaxedSubtractor = 0.0;
        for (int v = 0; v < G->nVertices; ++v)
        {
            deg = G->vertexDegrees[v] - vertexProhibitedDegrees[v];
			if(deg > 2) {
				constantRelaxedSubtractor += 2.0/deg;
			} 
        }

		//sort
		sort(edgeIndex, edgeIndex+G->mEdges, sort_indices(edgeWeight));
	}

	int getEdgeState(int e) {
		return edgeState[e];
	}

	bool activateEdge(int e) {	
		if(edgeState[e] != 0) {
			return false;
		}	
		//Vertices connected by edge
		int u = G->edges[e].first;
		int v = G->edges[e].second;
		UnionOperation op;

		bool doesntMakeCycle = vertexConnected->unionSet(u, v, op);
		
		//If this forms a cycle end with false
		if(!doesntMakeCycle){
			return false;
		}

		//store last union
		unionsList[unionCounter] = op;
		unionEdges[unionCounter] = e;
		unionCounter++;

		//activate edge
		edgeState[e] = 1;
		activeEdgeCount++;
        relaxedSolution[currentRelaxedSolution + 1] = relaxedSolution[currentRelaxedSolution] + edgeWeight[e];
        currentRelaxedSolution++;

		//Update branch vertice count
		if(vertexDegrees[u] == 2) {
			branchVertexCount++;
		}
		if(vertexDegrees[v] == 2) {
			branchVertexCount++;
		}

		//Increase vertices' degrees
		vertexDegrees[u]++; vertexDegrees[v]++;

		/*
		cout<<"  activate [";
		for (int i = 0; i < G->mEdges; i++) {
			cout<<i<<": "<<edgeState[i]<<", ";
		}
		cout<<"]"<<endl;
		*/

		return true;
	}

	void undoActivateEdge() {
		unionCounter--;

		vertexConnected->undoUnion(unionsList[unionCounter]);

		//deactivate edge
		int e = unionEdges[unionCounter];
		edgeState[e] = 0;
		activeEdgeCount--;
        currentRelaxedSolution--;

		//Decrease vertices' degrees
		//Vertices connected by edge
		int u = G->edges[e].first;
		int v = G->edges[e].second;
		vertexDegrees[u]--; vertexDegrees[v]--;

		//Update branch vertice count
		if(vertexDegrees[u] == 2) {
			branchVertexCount--;
		}
		if(vertexDegrees[v] == 2) {
			branchVertexCount--;
		}	

		/*
		cout<<"deactivate [";
		for (int i = 0; i < G->mEdges; i++) {
			cout<<i<<": "<<edgeState[i]<<", ";
		}
		cout<<"]"<<endl;	
		*/
	}

	void prohibitEdge(int e) {
		//Vertices connected by edge
		int u = G->edges[e].first;
		int v = G->edges[e].second;
		vertexProhibitedDegrees[u]++;
		vertexProhibitedDegrees[v]++;
		edgeState[e] = -1;

		/*
		cout<<"  prohibit [";
		for (int i = 0; i < G->mEdges; i++) {
			cout<<i<<": "<<edgeState[i]<<", ";
		}
		cout<<"]"<<endl;
		*/
	}

	void undoProhibitEdge(int e) {
		//Vertices connected by edge
		int u = G->edges[e].first;
		int v = G->edges[e].second;
		vertexProhibitedDegrees[u]--;
		vertexProhibitedDegrees[v]--;
		edgeState[e] = 0;

		/*
		cout<<"unprohibit [";
		for (int i = 0; i < G->mEdges; i++) {
			cout<<i<<": "<<edgeState[i]<<", ";
		}
		cout<<"]"<<endl;
		*/
	}

	int getBranchVertexCount() {
		return branchVertexCount;
	}

	int getActiveEdgeCount() {
		return activeEdgeCount;
	}

	int getDisjointSetCount() {
		return vertexConnected->getDisjointSetCount();
	}

    float getRelaxedSolution() {
        return relaxedSolution[currentRelaxedSolution] - constantRelaxedSubtractor;
    }
};