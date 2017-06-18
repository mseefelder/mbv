#include <limits>
#include <iostream>
#include <forward_list>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>
#include <chrono>
#include <random>

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
			cout<<"Already added all edges"<<endl;
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

		cout<<"Sorted!"<<endl;
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

int Backtrack(Graph &G, MBVSolutionUndo &sol, int &minimumSol, auto &start, long long &pruneCount)
{
	int e = 0;
	// Prune by inviability
	if (sol.getBranchVertexCount() >= minimumSol || minimumSol == 0) {
		return minimumSol;
	}
	// A valid solution was found
	else if (sol.getActiveEdgeCount() == G.nVertices-1) {
		if(sol.getBranchVertexCount() < minimumSol) {
			cout<<"Solution "<<sol.getBranchVertexCount()<<"(";
			auto end = chrono::steady_clock::now();
			cout << chrono::duration <double, milli> (end-start).count() << " ms)" << endl;
			minimumSol = sol.getBranchVertexCount();
		}
	}
	// No solution yet
	else{
		for (int i = 0; i<G.mEdges; i++) {
			e = sol.edgeIndex[i];

			if(sol.getEdgeState(e) == 0){
				if(sol.activateEdge(e)) {
                    minimumSol = Backtrack(G, sol, minimumSol, start, pruneCount);
					sol.undoActivateEdge();
				}

				sol.prohibitEdge(e);
				/**/
				
				//sol.sortEdges();
				// cout<<"Before Kruskal"<<endl;
				// cout<<"[";
				// for (int x = 0; x < G.mEdges; x++) {
				// 	cout<<x<<": "<<sol.edgeState[x]<<", ";
				// }
				// cout<<"]"<<endl;

				//Try kruskal ---
				int ke = 0;
				int opcount = 0;
				for (int k = 0; k < G.mEdges; ++k)
				{
					ke = sol.edgeIndex[k];
					if(sol.activateEdge(ke)) 
					{
						opcount++;
					}
					if(sol.getActiveEdgeCount() == G.nVertices - 1 ){
						break;
					}
				}

				// cout<<"[";
				// for (int x = 0; x < G.mEdges; x++) {
				// 	cout<<x<<": "<<sol.edgeState[x]<<", ";
				// }
				// cout<<"]"<<endl;

				if(sol.getBranchVertexCount() < minimumSol && sol.getActiveEdgeCount() == G.nVertices-1) {
					cout<<"Solution on Kruskal "<<sol.getBranchVertexCount()<<"(";
					auto end = chrono::steady_clock::now();
					cout << chrono::duration <double, milli> (end-start).count() << " ms)" << endl;
					minimumSol = sol.getBranchVertexCount();
				}

                //Check if we should branch
                bool doBranch = (sol.getRelaxedSolution() < minimumSol) && (sol.getActiveEdgeCount() == G.nVertices-1);
                if(!doBranch)
                    pruneCount++;

				for (int x = 0; x < opcount; ++x)
				{
					sol.undoActivateEdge();
				}

				// cout<<"[";
				// for (int x = 0; x < G.mEdges; x++) {
				// 	cout<<x<<": "<<sol.edgeState[x]<<", ";
				// }
				// cout<<"]"<<endl;
				//            ---
				/**/

				//Only branch if relaxed bound is < minimumSol
                if(doBranch)
                    minimumSol = Backtrack(G, sol, minimumSol, start, pruneCount);

				sol.undoProhibitEdge(e);
				//sol.sortEdges();

				break;
			}

		}
	}

	return minimumSol;
}

int MBVGrasp (int maxIter, int randomMargin) {
	// Random setup
	std::random_device r;
	std::default_random_engine engine(r());
	std::uniform_int_distribution<int> uDist(0, randomMargin);

	// Create solution
	MBVSolutionUndo gS(G);

	// Make edges connecting vertices with degree 1 obligatory
	// Number of obligatory edges
	int obl = 0;
	for (int i = 0; i < nVertices; ++i)
	{
		if (G.vertexDegrees[i] == 1)
		{
			gS.activateEdge(G.vertices[i].front());
			obl++;
		}
	}

	// Run maxIter iterations
	for (int iter = 0; iter < maxIter; ++iter)
	{
		// Find randomized greedy solution gS
		// Fill vector with edges
		vector<int> edges;
		edges.reserve(G.mEdges);
		for (int e = 0; e < G.mEdges; ++e) {
			edges.push_back(gS.edgeIndex[e]);
		}

		// Create vector for cycle-making edges
		queue<int> cycleEdges;

		// Auxiliar structure for local search (represents the tree)
		int parent[G.nVertices];
		int treeDegrees[G.nVertices];
		int rank[G.nVertices];
		for (int i = 0; i < g.nVertices; ++i)
		{
			parent[i] = -1;
			treeDegrees[i] = gS.vertexDegrees[i];
			rank[i] = 0;
		}

		// Choose semi-randomly (randomMargin) with Kruskal
		int e = 0, u = 0, v = 0;
		for (int i = 0; i < g=G.mEdges; ++i) {
			// get random index (from end of vector to -randomMargin positions)
			e = uDist(engine);
			// if there are less than randomMargin elements, fix index
			if (e >= edges.size()) {
				e = 0;
			}
			// get value of index on edges vector
			e = *(edges.rbegin()+e);
			// try to add edge to solution
			if(!gS.activateEdge(e)) {
				cycleEdges.push(e);
			}
			else {
				// Get vertices connected by edge and their degrees onthe tree so far
				u = G.edges[e].first;
				v = G.edges[e].second;
				treeDegrees[u]++;
				treeDegrees[v]++;
				// Add edge to auxiliar structure for local search (in the form of a parent relation)
				if(treeDegrees[v] > treeDegrees[u]) {
					parent[u] = v;
					rank[u] = rank[v]+1;
				}
				else {
					parent[v] = u;
					rank[v] = rank[u]+1;
				}
			}

			// Remove edge from vector of edges left
			edges.erase(edges.rbegin()+e);

			// If kruskal is over
			if (gS.getActiveEdgeCount == G.nVertices -1) {
				break;
			}
		}

		// Now, we have a solution on gS. 
		// On edges vector we have the edges that weren't used on Kruskal
		// On cycleEdges we have the ones that would form a cycle

		// Current branch vertex count
		int branchVertexCount = gS.getBranchVertexCount(); /// ALWAYS UPDATE THIS WHEN UPDATING TREE

		// Add edges left from edges vector to cycleEdges vector
		while(!edges.empty()) {
		    cycleEdges.push(edges.back());
		    edges.pop_back();
		}

		// Number of iterations that will be executed to look for improvement
		int remaining = cycleEdges.size();
		// Number of edges reinserted on 
		int inserted = 0;

		// Local search starting on gS
		// Iterate over edges to evaluate local search neighbourhood
		e = 0, u = 0, v = 0;
		while(remaining!=0) {
			// Remove last element from vector
			e = cycleEdges.front();
			cycleEdges.pop();

			// Get edge's vertices
			u = G.edges[e].first;
			v = G.edges[e].second;

			bool keepChange = false;
			/// Virtually add e to tree, forming a cycle C
						
			/// Get sequence SEQ of edges in cycle C
			
			/// For each edge e'(i,j) on SEQ:
				/// If removing e' from tree reduces the degree of (i || j) from 3 to 2
					/// Remove e' from tree
					/// Add e' to cycleEdges
					/// (remaining = inserted+1)
					/// (inserted = 0)
					/// Make sure e will be added to tree (keepChange = true)
					/// Break

			/// If no improvement was encountered (keepChange == false):
				/// Undo adding e to tree 
				/// Add e to cycleEdges
				/// (inserted++)

			/// (remaining--)
		}

		/// Save results and reset solution
		/// Undo nVertices - 1 - obl activations
		
	}
}

int main(int argc, char const *argv[])
{
	int nVertices, mEdges, trash;

	cin >> nVertices >> mEdges >> trash;
	cout << nVertices << " vertices and " << mEdges << " edges."<<endl;

	Graph G(nVertices, mEdges);

	int a, b;
	for (int e = 0; e <= mEdges; e++) {
		cin >> a >> b >> trash;
		G.addEdge(a-1, b-1);
	}

	int minimumSol = numeric_limits<int>::max();
	MBVSolutionUndo S(G);
	MBVSolutionUndo initialKruskal(G);

	//make bridges obligatory
	for (int i = 0; i < nVertices; ++i)
	{
		if (G.vertexDegrees[i] == 1)
		{
			S.activateEdge(G.vertices[i].front());
		}
	}

	//Run Kruskal to get a headstart
	int e = 0;
	for (int i = 0; i < mEdges; ++i)
	{
		e = initialKruskal.edgeIndex[i];
		bool activated = initialKruskal.activateEdge(e);

        cout<<e<<": "<<initialKruskal.edgeWeight[e]<<"(";
        cout<<G.vertexDegrees[G.edges[e].first]<<","<<G.vertexDegrees[G.edges[e].second];
        cout<<") "<<activated<<endl;

		if(initialKruskal.getActiveEdgeCount() == nVertices - 1 ){
            cout<<endl;
            cout<<"Relaxed solution "<<initialKruskal.getRelaxedSolution()<<endl;
			break;
		}
	}

	minimumSol = initialKruskal.getBranchVertexCount();
	cout<<"Initial heuristic: "<<minimumSol<<endl;

    for (int i = 0; i < nVertices-1; i++) {
        initialKruskal.undoActivateEdge();
    }
    cout<<"Relaxed solution "<<initialKruskal.getRelaxedSolution()<<endl;

    long long pruneCount = 0;

	auto start = chrono::steady_clock::now();
	minimumSol = Backtrack(G, S, minimumSol, start, pruneCount);
	auto end = chrono::steady_clock::now();
	cout << "Solution has " << minimumSol << " branch vertices! (";
	cout << chrono::duration <double, milli> (end-start).count() << " ms), ";
    cout << pruneCount << " prunes." << endl;
	return 0;
}