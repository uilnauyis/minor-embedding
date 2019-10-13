/*
  Heuristic search Minor embedding algorithm <zyan632@aucklanduni.ac.nz> 
  Basic version with new improvements (random selection)
  usage: embedding.exe guest 
  put host graph (host.in) in the same folder
*/

#include "Graph.h"
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <ctime>
#include <fstream>
#include <bits/stdc++.h>

// random generator function:
vector<int> randomizeVertexOrder(int size)
{
  vector<int> myvector;
  for (int i = 0; i < size; ++i)
  {
    myvector.push_back(i); //set value 0 1 2 3 ... size-1
  }
  random_shuffle(myvector.begin(), myvector.end()); // using built-in random generator:
  return myvector;
}

// sort vertex order by size
class Vertex
{
public:
  int index;
  int size;
};
bool sortByVertexSize(const Vertex &v1, const Vertex &v2)
{
  return v1.size > v2.size;
}

vector<int> sortVertexOrder(Graph &H)
{
  vector<Vertex> toSort;
  vector<int> myvector;
  for (int i = 0; i < H.order(); i++)
  {
    Vertex v;
    v.index = i;
    v.size = H.deg(i);
    toSort.push_back(v);
  }
  sort(toSort.begin(), toSort.end(), sortByVertexSize);
  for (int i = 0; i < H.order(); i++)
    myvector.push_back(toSort.at(i).index);
  return myvector;
}

// initialize the mapping
vector<vector<int>> initializeMapping(int size)
{
  vector<vector<int>> mapping;
  for (int i = 0; i < size; i++)
  {
    vector<int> empty;
    mapping.push_back(empty);
  }
  return mapping;
}

// improving 1, minimize max number of overlap
int getMaxOverlap(vector<vector<int>> qubitsToVerticesMapping)
{
  int max = 0;
  for (int i = 0; i < qubitsToVerticesMapping.size(); i++)
  {
    int thisOverlap = qubitsToVerticesMapping.at(i).size();
    max = thisOverlap > max ? thisOverlap : max;
  }
  return max;
}

// improving 2, minimize total number of physical qubits
int getTotalSize(vector<vector<int>> mapping)
{
  int total = 0;
  for (int i = 0; i < mapping.size(); i++)
  {
    total += mapping.at(i).size();
  }
  return total;
}

// define vertex weight - weight of a vertex in G grows exponentialy with num of vertices of H represented there
int getQubitWeight(int overlapNum, bool isIn)
{
  int DIAMETER = 30;
  if (isIn)
    return pow(DIAMETER, overlapNum - 1);
  return pow(DIAMETER, overlapNum);
}

// check all neighbors empty
bool checkAllNeighbors(Graph &H, vector<vector<int>> mapping, int currentVertex)
{
  for (int i = 0; i < H.deg(currentVertex); i++)
  {
    if (mapping.at(H.adj[currentVertex].at(i)).size() != 0)
      return false;
  }
  return true;
}

// find the weight for shortest path from A to B
void dijkstraComputePaths(int source, Graph &adjacency_list, vector<int> &previous,
                          vector<int> &min_distance,
                          vector<vector<int>> qubitToVertexMapping,
                          vector<vector<int>> vertexToQubitMapping,
                          int currentVertex)
{
  int n = adjacency_list.order();
  min_distance.clear();
  min_distance.resize(n, 1 << 30);
  int sourceQubitOverlappedNum = qubitToVertexMapping.at(source).size();
  min_distance[source] = 0;
  previous.clear();
  previous.resize(n, -1);
  set<pair<int, int>> vertex_queue;
  vertex_queue.insert(make_pair(min_distance[source], source));

  while (!vertex_queue.empty())
  {
    int dist = vertex_queue.begin()->first;
    int u = vertex_queue.begin()->second;
    vertex_queue.erase(vertex_queue.begin());

    // The shortest path from root to current vertex is already find.
    if (dist > min_distance[adjacency_list.adj.size() - 1])
    {
      return;
    }

    // Visit each edge exiting u
    for (int i = 0; i < adjacency_list.deg(u); i++)
    {
      int v = adjacency_list.adj[u].at(i);
      int weight;
      if (v >= qubitToVertexMapping.size()) // v is dummy qubit
      {
        weight = 0;
      }
      else
      {
        int overlappedVertexNum = qubitToVertexMapping.at(v).size();
        bool isIn = false;
        vector<int> currentVertexQubits = vertexToQubitMapping.at(currentVertex);
        if (std::find(currentVertexQubits.begin(),
                      currentVertexQubits.end(), v) == currentVertexQubits.end())
        {
          isIn = true;
        }
        weight = getQubitWeight(overlappedVertexNum + 1, isIn);
      }
      int distance_through_u = dist + weight;
      if (distance_through_u < min_distance[v])
      {
        vertex_queue.erase(make_pair(min_distance[v], v));
        min_distance[v] = distance_through_u;
        previous[v] = u;
        vertex_queue.insert(make_pair(min_distance[v], v));
      }
    }
  }
}

// find the shortest
vector<int> DijkstraGetShortestPathTo(
    int vertex, const std::vector<int> &previous)
{
  vector<int> path;
  for (; vertex != -1; vertex = previous[vertex])
    path.push_back(vertex);
  return path;
}

//create a dummy graph, add a dummy vertex which adjacent to every vertex in the model
void createDummyGraph(Graph &dummyGraph, Graph &G, Graph &H,
                      vector<vector<int>> vertexToQubitsMapping,
                      int currentVertex, int index)
{
  dummyGraph = G;
  dummyGraph.adj.resize(G.order() + 1);
  dummyGraph.adj[G.order()].resize(0);

  for (int k = 0; k < vertexToQubitsMapping.at(H.adj[currentVertex].at(index)).size(); k++)
  {
    dummyGraph.addArc(G.order(), vertexToQubitsMapping.at(H.adj[currentVertex].at(index)).at(k));
    dummyGraph.addArc(vertexToQubitsMapping.at(H.adj[currentVertex].at(index)).at(k), G.order());
  }
}

void expand(vector<vector<int>> &qubitToVertexMapping,
            vector<vector<int>> &verticesToQubitsMapping,
            vector<vector<int>> &verticesToAvailableEdgeMapping,
            set<int> &affectedVertices,
            int currentQubit,
            int currentVertex,
            Graph G)
{
  // Add random qubit to the current vertex model.
  qubitToVertexMapping.at(currentQubit).push_back(currentVertex);
  verticesToQubitsMapping.at(currentVertex).push_back(currentQubit);

  vector<int> randomQubitNeighborQubits = G.adj.at(currentQubit);
  for (int neighborIndex = 0; neighborIndex < randomQubitNeighborQubits.size();
       neighborIndex++)
  {
    int neighborQubit = randomQubitNeighborQubits.at(neighborIndex);
    vector<int> neighborQubitVertices = qubitToVertexMapping.at(neighborQubit);

    // Mark the random qubit as unavailable for all neighbor vertex model
    // as it has now been added to current vertex model.
    if (neighborQubitVertices.size() > 0)
    {
      for (int neighborQubitVertexIndex = 0;
           neighborQubitVertexIndex < neighborQubitVertices.size();
           neighborQubitVertexIndex++)
      {
        int neighborQubitVertex = neighborQubitVertices.at(neighborQubitVertexIndex);
        vector<int> &neighborQubitVertexAvailableEdges =
            verticesToAvailableEdgeMapping.at(neighborQubitVertex);
        neighborQubitVertexAvailableEdges.erase(
            std::remove(neighborQubitVertexAvailableEdges.begin(),
                        neighborQubitVertexAvailableEdges.end(),
                        currentQubit));

        // As this vertix has been affected by the extension on current vertex,
        // we save this neighbor vertex and extend it later
        if (neighborQubitVertex != currentVertex)
        {
          affectedVertices.insert(neighborQubitVertex);
        }
      }
    }
    // If the neighbor qubit is not assigned to any vertex model yet,
    // it is added to the 'verticesToAvailableEdgeMapping' for the current vertex model.
    else
    {
      vector<int> &currentVertexAvailableEdges =
          verticesToAvailableEdgeMapping.at(currentVertex);
      currentVertexAvailableEdges.push_back(neighborQubit);
    }
  }
}

/* 
    Extend a partial embedding to satisfy the ratio.
*/
bool extendChain(Graph &G, Graph &H, vector<vector<int>> &verticesToQubitsMapping,
                 vector<vector<int>> &qubitToVertexMapping,
                 vector<vector<int>> &verticesToAvailableEdgeMapping,
                 int currentVertex, float ratio)
{
  std::set<int> affectedVertices;
  while (true)
  {
    // Calculate unembedded neighbors of the current vertex.
    int unembeddedNeighborVerticesNum = 0;
    vector<int> currentNeighborVertices = H.adj.at(currentVertex);
    for (int neighIndex = 0; neighIndex < currentNeighborVertices.size(); neighIndex++)
    {
      int neigh = currentNeighborVertices.at(neighIndex);
      if (verticesToQubitsMapping.at(neigh).size() > 0)
      {
        continue;
      }
      unembeddedNeighborVerticesNum++;
    }
    vector<int> currentVertexAvailableEdges = verticesToAvailableEdgeMapping.at(currentVertex);
    int availableEdgesNum = currentVertexAvailableEdges.size();

    // If there are sufficient available edges, there is no need to extend current
    // vertex model.
    if (availableEdgesNum >= std::ceil(ratio * unembeddedNeighborVerticesNum))
    {
      break;
    }

    // Get an available that adjacent to current vertex model randomly.
    int randomQubit = currentVertexAvailableEdges.at(rand() % availableEdgesNum);

    expand(qubitToVertexMapping, verticesToQubitsMapping,
           verticesToAvailableEdgeMapping, affectedVertices,
           randomQubit, currentVertex, G);
  }

  std::set<int>::iterator it;
  for (it = affectedVertices.begin(); it != affectedVertices.end(); ++it)
  {
    int neighborQubit = *it;
    vector<int> neighborQubitVertices = qubitToVertexMapping.at(neighborQubit);
    if (neighborQubitVertices.size() > 0)
    {
      for (int neighborQubitVertexIndex = 0;
           neighborQubitVertexIndex < neighborQubitVertices.size();
           neighborQubitVertexIndex++)
      {
        int neighborQubitVertex = neighborQubitVertices.at(neighborQubitVertexIndex);
        extendChain(G, H, verticesToQubitsMapping, qubitToVertexMapping,
                    verticesToAvailableEdgeMapping, neighborQubitVertex, ratio);
      }
    }
  }
}

void extendAll(Graph G, Graph H, float ratio,
               std::set<int> affectedVertices, vector<vector<int>> qubitToVertexMapping,
               vector<vector<int>> verticesToQubitsMapping,
               vector<vector<int>> verticesToAvailableEdgeMapping)
{
  std::set<int>::iterator it;
  for (it = affectedVertices.begin(); it != affectedVertices.end(); ++it)
  {
    extendChain(G, H, verticesToQubitsMapping, qubitToVertexMapping,
                verticesToAvailableEdgeMapping, *it, ratio);
  }
}

//find the union paths from root g* to each vertex model
void updateModels(Graph &G, Graph &H, int root, int currentVertex,
                  vector<vector<int>> &qubitToVertexMapping,
                  vector<vector<int>> &vertexToQubitMapping,
                  vector<vector<int>> &vertexToAvailableEdges,
                  float ratio)
{
  map<int, int> vertex_usage;
  int dist[G.order()];

  vector<int> path;
  for (int i = 0; i < H.deg(currentVertex); i++)
  {
    path.clear();
    if (vertexToQubitMapping.at(H.adj[currentVertex].at(i)).size() != 0)
    {
      Graph dummyGraph(0);
      createDummyGraph(dummyGraph, G, H, vertexToQubitMapping, currentVertex, i);
      vector<int> min_distance;
      vector<int> previous;
      dijkstraComputePaths(root, dummyGraph, previous, min_distance,
                           qubitToVertexMapping, vertexToQubitMapping, currentVertex);
      path = DijkstraGetShortestPathTo(G.order(), previous);
    }
    for (int j = 2; j < path.size(); j++)
    { // don't want include the dummy vertex and the vertex in neighbor's map
      if (vertex_usage.find(path.at(j)) == vertex_usage.end())
      {
        vertex_usage.insert(pair<int, int>(path.at(j), 1)); // if not exist, create a new one
        dist[path.at(j)] = H.adj[currentVertex].at(i);      // destination
      }
      else
      {
        vertex_usage[path.at(j)]++; // update the usage
      }
    }
  }

  std::set<int> affectedVertices;
  //update all vertex models

  for (map<int, int>::iterator ii = vertex_usage.begin(); ii != vertex_usage.end(); ++ii)
  {
    if ((*ii).first == root)
    {
      expand(qubitToVertexMapping, vertexToQubitMapping,
             vertexToAvailableEdges, affectedVertices, (*ii).first,
             currentVertex, G);
      continue;
    }
    if ((*ii).second > 1)
    { // appears in more than one path, assign to new model ****************************
      expand(qubitToVertexMapping, vertexToQubitMapping,
             vertexToAvailableEdges, affectedVertices, (*ii).first,
             currentVertex, G);
    }

    else
    { // appears once, update
      if (find(vertexToQubitMapping.at(dist[(*ii).first]).begin(),
               vertexToQubitMapping.at(dist[(*ii).first]).end(),
               (*ii).first) == vertexToQubitMapping.at(dist[(*ii).first]).end())
      {
        expand(qubitToVertexMapping, vertexToQubitMapping,
               vertexToAvailableEdges, affectedVertices, (*ii).first,
               dist[(*ii).first], G);
      }
      else
      {
        cout << "meow";
      }
    }
  }
  extendAll(G, H, ratio, affectedVertices, qubitToVertexMapping,
            vertexToQubitMapping, vertexToAvailableEdges);
}

void tearChain(Graph &G, Graph &H, vector<vector<int>> &vertexToQubitMapping,
               vector<vector<int>> &qubitToVertexMapping,
               vector<vector<int>> &verticesToAvailableEdgeMapping,
               int currentVertex, float ratio)
{
  vector<int> currentVertexQubits = vertexToQubitMapping.at(currentVertex);
  for (int currentVertexQubitIndex = 0;
       currentVertexQubitIndex < currentVertexQubits.size();
       currentVertexQubitIndex++)
  {
    int currentVertexQubit = currentVertexQubits.at(currentVertexQubitIndex);

    // This qubit is not overlapped.
    if (qubitToVertexMapping.at(currentVertexQubit).size() <= 1)
    {
      vector<int> neighborQubits = G.adj.at(currentVertexQubit);
      for (int neighborQubitIndex = 0;
           neighborQubitIndex < neighborQubits.size();
           neighborQubitIndex++)
      {
        int neighborQubit = neighborQubits.at(neighborQubitIndex);
        vector<int> neighborQubitVertices = qubitToVertexMapping.at(neighborQubit);
        if (qubitToVertexMapping.at(neighborQubit).size() > 0)
        // This neighbor qubit not in current vertex and is mapped already.
        {
          for (int neighborQubitVertexIndex = 0;
               neighborQubitVertexIndex < neighborQubitVertices.size();
               neighborQubitVertexIndex++)
          {
            int neighborQubitVertex = neighborQubitVertices.at(neighborQubitVertexIndex);
            if (neighborQubitVertex == currentVertex)
            {
              continue;
            }
            vector<int> &neighborQubitVertexAvailableEdges =
                verticesToAvailableEdgeMapping.at(neighborQubitVertex);
            neighborQubitVertexAvailableEdges.push_back(currentVertexQubit);
          }
        }
      }
    }
    vector<int> &currentMappedVertices = qubitToVertexMapping.at(currentVertexQubit);
    currentMappedVertices.erase(
        std::remove(currentMappedVertices.begin(),
                    currentMappedVertices.end(), currentVertex),
        currentMappedVertices.end());
  }
  vertexToQubitMapping.at(currentVertex).clear();
  verticesToAvailableEdgeMapping.at(currentVertex).clear();
}

//find the minimal vertex model
void findMinimalVertexModel(Graph &G, Graph &H, vector<vector<int>> &vertexToQubitMapping,
                            vector<vector<int>> &qubitToVertexMapping,
                            vector<vector<int>> &verticesToAvailableEdgeMapping,
                            int currentVertex, float ratio)
{
  // if all neighbors are empty
  if (checkAllNeighbors(H, vertexToQubitMapping, currentVertex))
  {
    vector<int> qubitOrder = randomizeVertexOrder(G.order());
    for (int qubitIndex = 0; qubitIndex < qubitOrder.size(); qubitIndex++)
    {
      int qubit = qubitOrder.at(qubitIndex);
      int vOverlappedNum = qubitToVertexMapping.at(qubit).size();
      if (vOverlappedNum == 0 && G.adj[qubit].size() != 0) // Pick unisolated qubits
      {
        std::set<int> affectedVertices;
        affectedVertices.insert(currentVertex);
        expand(qubitToVertexMapping, vertexToQubitMapping, verticesToAvailableEdgeMapping,
               affectedVertices, qubit, currentVertex, G);

        // After randomly picking one qubit for this mapping,
        // we extend the mapping to maintain the ratio.
        extendAll(G, H, ratio, affectedVertices, qubitToVertexMapping,
                  vertexToQubitMapping, verticesToAvailableEdgeMapping);
        return;
      }
    }
  }

  // tear the previous embedding for current vertex.
  tearChain(G, H, vertexToQubitMapping, qubitToVertexMapping, verticesToAvailableEdgeMapping, currentVertex, ratio);

  // new improvement:  choose Candidate from first level
  vector<int> rootCandidate;
  vector<int> rootCandidateWeight;
  for (int neighborVertexIndex = 0;
       neighborVertexIndex < H.deg(currentVertex);
       neighborVertexIndex++)
  {
    int neighborVertex = H.adj[currentVertex].at(neighborVertexIndex);
    vector<int> neighborVertexQubits = vertexToQubitMapping.at(neighborVertex);
    for (int neighborVertexQubitIndex = 0;
         neighborVertexQubitIndex < neighborVertexQubits.size();
         neighborVertexQubitIndex++)
    {
      int neighborVertexQubit = neighborVertexQubits.at(neighborVertexQubitIndex);
      vector<int> firstLevelQubits = G.adj.at(neighborVertexQubit);
      for (int potentialCandidateIndex = 0;
           potentialCandidateIndex < firstLevelQubits.size();
           potentialCandidateIndex++)
      {
        int potentialCandidate = firstLevelQubits.at(potentialCandidateIndex);
        // ------------ Could be improved. -----------------
        if (find(rootCandidate.begin(),
                 rootCandidate.end(),
                 potentialCandidate) == rootCandidate.end() &&
            qubitToVertexMapping.at(potentialCandidate).size() == 0) // ?? Should I restrict this?
        {                                                            // if not contain
          rootCandidate.push_back(potentialCandidate);
        }
      }
    }
  }

  int cost[rootCandidate.size()][H.deg(currentVertex)];
  int minCost = 1 << 30;
  int root = 0;
  for (int i = 0; i < rootCandidate.size(); i++)
  {
    int sumCost = 0;
    int currentRootCandidate = rootCandidate.at(i);
    for (int j = 0; j < H.deg(currentVertex); j++)
    {
      // This vertex not imbedded yet.
      if (vertexToQubitMapping.at(H.adj[currentVertex].at(j)).size() == 0)
      {
        cost[i][j] = 0;
      }
      // Root candidate qubit is overlapped with this vertex.
      else if (find(vertexToQubitMapping.at(H.adj[currentVertex].at(j)).begin(),
                    vertexToQubitMapping.at(H.adj[currentVertex].at(j)).end(),
                    currentRootCandidate) !=
               vertexToQubitMapping.at(H.adj[currentVertex].at(j)).end())
      {
        cost[i][j] = 0;
      }
      else
      {
        // add a dummy vertex which adjacent to every vertex in the model - compute the distance from a vertex to a vertex-model
        Graph dummyGraph(0);
        createDummyGraph(dummyGraph, G, H, vertexToQubitMapping,
                         currentVertex, j);
        vector<int> min_distance;
        vector<int> previous;
        dijkstraComputePaths(currentRootCandidate, dummyGraph, previous, min_distance,
                             qubitToVertexMapping, vertexToQubitMapping, currentVertex);
        cost[i][j] = min_distance.at(G.order());
      }
      sumCost += cost[i][j];
      if (sumCost > minCost)
        break;
    }
    // Add the cost of root to sum of weights.
    sumCost += getQubitWeight(
        qubitToVertexMapping.at(currentRootCandidate).size() + 1, false);
    rootCandidateWeight.push_back(sumCost);
    // find root g* with minimum cost
    if (sumCost <= minCost)
    {
      minCost = sumCost;
    }
  }

  // random get the root from the first level rootCandidate
  vector<int> candidates;
  for (int i = 0; i < rootCandidate.size(); i++)
  {
    if (minCost == rootCandidateWeight.at(i))
      candidates.push_back(rootCandidate.at(i));
  }

  // random get the root from the candidates
  if (candidates.size() == 0)
  {
    vector<int> qubitOrder = randomizeVertexOrder(G.order());
    for (int i = 0; i < G.order(); i++)
    {
      if (G.deg(qubitOrder.at(i)) == 0)
      { // isolated vertex, ingore
        continue;
      }
      if (qubitToVertexMapping.at(qubitOrder.at(i)).size() <= 0)
      { // if no candidates, random choose one
        root = qubitOrder.at(i);
        break;
      }
    }
  }
  else
  {
    random_shuffle(candidates.begin(), candidates.end()); // using built-in random generator:
    root = candidates.at(0);
  }

  //find the union paths from root g* to each vertex model
  updateModels(G, H, root, currentVertex, qubitToVertexMapping,
               vertexToQubitMapping, verticesToAvailableEdgeMapping, ratio);
}

// main method for finding Minor embedding
bool findMinorEmbedding(Graph &G, Graph &H, float ratio)
{
  srand(unsigned(time(0))); // for randomizing
#if 1
  vector<int> vertexOrder = randomizeVertexOrder(H.order());
#else
  vector<int> vertexOrder = sortVertexOrder(H); // sort order by size
#endif
  vector<vector<int>> verticesToQubitsMapping = initializeMapping(H.order());
  vector<vector<int>> qubitsToVerticesMapping = initializeMapping(G.order());
  vector<vector<int>> verticesToAvailableEdgeMapping =
      initializeMapping(H.order());
  int previousChainLength;
  int previousOverlap;
  int currentChainLength;
  int currentOverlap;
  int stage = 1;
  while (true)
  {
    currentChainLength = stage <= 1 ? INT_MAX : getTotalSize(verticesToQubitsMapping);
    currentOverlap = stage <= 1 ? INT_MAX : getMaxOverlap(qubitsToVerticesMapping);
    if (previousChainLength <= currentChainLength &&
        previousOverlap <= currentOverlap &&
        stage > 2) // No improvement in last update.
    {
      break;
    }
    previousChainLength = currentChainLength;
    previousOverlap = currentOverlap;
    for (int i = 0; i < vertexOrder.size(); i++)
    {
      int currentVertex = vertexOrder.at(i);
      findMinimalVertexModel(G, H, verticesToQubitsMapping, qubitsToVerticesMapping,
                             verticesToAvailableEdgeMapping, currentVertex, ratio);
    }
    stage++;
  }

  //check the usage for each vertex in G, if all <=1, then we find the embedding. otherwise "failure"
  bool embedding = true;
  //print the embedding, if find
  if (embedding)
  {
    int totalQubit = 0;
    int max = 0;
    cout << "[";
    for (int i = 0; i < verticesToQubitsMapping.size(); i++)
    {
      cout << "[";
      totalQubit = totalQubit + verticesToQubitsMapping.at(i).size();
      if (verticesToQubitsMapping.at(i).size() > max)
        max = verticesToQubitsMapping.at(i).size();
      sort(verticesToQubitsMapping.at(i).begin(), verticesToQubitsMapping.at(i).end());
      for (int j = 0; j < verticesToQubitsMapping.at(i).size() - 1; j++)
      {
        cout << verticesToQubitsMapping.at(i).at(j) << ", ";
      }
      cout << verticesToQubitsMapping.at(i).at(verticesToQubitsMapping.at(i).size() - 1);
      if (i == verticesToQubitsMapping.size() - 1)
      {
        cout << "]";
      }
      else
      {
        cout << "], " << endl;
      }
    }
    cout << "]" << endl;
    cout << "Total number of Qubits used is: " << totalQubit << endl;
    cout << "Max chain length is: " << max << endl;
    return true;
  }
  else
  {
    return false;
  }
}

int main(int argc, char **argv)
{
  clock_t start = clock();
  ifstream host("host.in");
  Graph G(0);
  host >> G;
  Graph H(0);
  cin >> H;
  float ratio = std::stof(argv[1]);
  bool foundEmbedding;
  cout << "Guest graph: " << endl;
  cout << H << endl;
  for (int i = 1; i <= 2; i++)
  {
    cout << "try " << i << ":" << endl;
    foundEmbedding = findMinorEmbedding(G, H, ratio);
    if (foundEmbedding)
    {
      break;
    }
  }
  if (!foundEmbedding)
    cout << "No embedding found" << endl;

  clock_t ends = clock();
  //cout << "Running Time : " << (double)(ends - start) / CLOCKS_PER_SEC << endl;
  return 0;
}
