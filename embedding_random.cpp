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
#include <queue>
#include <utility>
#include <functional>
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

vector<set<int>> initializeSetMapping(int size)
{
  vector<set<int>> mapping;
  for (int i = 0; i < size; i++)
  {
    set<int> empty;
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

// improving 2, minimize total number of physical qubits
int getTotalSetSize(vector<set<int>> mapping)
{
  int total = 0;
  for (int i = 0; i < mapping.size(); i++)
  {
    total += mapping.at(i).size();
  }
  return total;
}

// define vertex weight - weight of a vertex in G grows exponentialy with num of vertices of H represented there
int getQubitWeight(int overlapNum)
{
  int DIAMETER = 30;
  return pow(DIAMETER, overlapNum - 1);
}

// check all neighbors empty
bool checkAllNeighbors(Graph &H, vector<set<int>> mapping, int currentVertex)
{
  for (int i = 0; i < H.deg(currentVertex); i++)
  {
    if (mapping.at(H.adj[currentVertex].at(i)).size() != 0)
      return false;
  }
  return true;
}

auto selectrandomFromSet(set<int> set)
{
  if (set.size() == 0)
  {
    return -1; // No neighbor qubit available
  }
  auto it = set.begin();
  auto rounds = rand() % set.size();
  std::advance(it, rounds);
  return *it;
}

// Compute shortest paths from current vertex to all neighbor vertices.
// Return the total cost as result.
// 'previous' is the vector for tracing back the path.
int dijkstraComputePaths(int source, int currentVertex,
                         Graph &G, Graph &H,
                         vector<int> &previous,
                         vector<tuple<int, int, int>> &vertexTuples, // <vertex, qubit, distance>
                         vector<vector<int>> &qubitToVertexMapping,
                         vector<set<int>> &vertexToQubitsMapping)
{
  int n = G.order();
  vector<int> min_distance;
  min_distance.resize(n, 1 << 30);
  min_distance[source] = 0;
  previous.clear();
  previous.resize(n, -1);
  set<pair<int, int>> vertex_queue;
  set<int> visited;

  for (int j = 0; j < H.deg(currentVertex); j++)
  {
    // This vertex not imbedded yet.
    if (vertexToQubitsMapping.at(H.adj[currentVertex].at(j)).size() == 0)
    {
      continue;
    }
    // Root candidate qubit is overlapped with this vertex.
    else if (find(vertexToQubitsMapping.at(H.adj[currentVertex].at(j)).begin(),
                  vertexToQubitsMapping.at(H.adj[currentVertex].at(j)).end(),
                  source) !=
             vertexToQubitsMapping.at(H.adj[currentVertex].at(j)).end())
    {
      continue;
    }
    else
    {
      int thisVertex = H.adj[currentVertex].at(j);
      vertexTuples.push_back(make_tuple(thisVertex, -1, INT_MAX));
    }
  }

  vertex_queue.insert(make_pair(min_distance[source], source));

  while (!vertex_queue.empty())
  {
    int dist = vertex_queue.begin()->first;
    int u = vertex_queue.begin()->second;
    vertex_queue.erase(vertex_queue.begin());
    visited.insert(u);

    // End the function if dist is already greater than or equal to
    // the distances to all neighbor vertices.
    bool canReturn = true;
    for (vector<tuple<int, int, int>>::iterator vertexTupleIt = vertexTuples.begin();
         vertexTupleIt != vertexTuples.end(); ++vertexTupleIt)
    {
      if (get<2>(*vertexTupleIt) > dist)
      {
        canReturn = false;
        break;
      }
    }
    if (canReturn)
    {
      break;
    }

    // Visit each edge exiting u
    for (int i = 0; i < G.deg(u); i++)
    {
      int v = G.adj[u].at(i);

      // If v has been visited, ignore
      if (visited.find(v) != visited.end())
      {
        continue;
      }

      // Check if this qubit is mapped to any neighbor vertex.
      // If it is, we update the shortest distance to this vertex
      vector<int> overlappedVertices = qubitToVertexMapping.at(v);
      for (vector<tuple<int, int, int>>::iterator vertexTuplesIt = vertexTuples.begin();
           vertexTuplesIt != vertexTuples.end(); ++vertexTuplesIt)
      {
        tuple<int, int, int> &vertexTuple = *vertexTuplesIt;
        if (find(overlappedVertices.begin(),
                 overlappedVertices.end(),
                 get<0>(vertexTuple)) != overlappedVertices.end() &&
            get<2>(vertexTuple) > dist)
        {
          get<1>(*vertexTuplesIt) = v;
          get<2>(*vertexTuplesIt) = dist;
        }
      }

      int weight = getQubitWeight(overlappedVertices.size() + 1);
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

  int sumCost = 0;
  for (vector<tuple<int, int, int>>::iterator vertexTuplesIt = vertexTuples.begin();
       vertexTuplesIt != vertexTuples.end(); ++vertexTuplesIt)
  {
    sumCost += get<2>(*vertexTuplesIt);
  }

  // Add the cost of root to sum of weights.
  sumCost += getQubitWeight(
      qubitToVertexMapping.at(source).size() + 1);
  return sumCost;
}

void expand(vector<vector<int>> &qubitToVertexMapping,
            vector<set<int>> &vertexToQubitsMapping,
            vector<set<int>> &vertexToAvailableEdgesMapping,
            set<int> &affectedVertices,
            int currentQubit,
            int currentVertex,
            Graph G)
{
  // Add random qubit to the current vertex model.
  qubitToVertexMapping.at(currentQubit).push_back(currentVertex);
  vertexToQubitsMapping.at(currentVertex).insert(currentQubit);

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
        set<int> &neighborQubitVertexAvailableEdges =
            vertexToAvailableEdgesMapping.at(neighborQubitVertex);
        neighborQubitVertexAvailableEdges.erase(currentQubit);

        // As this vertix has been affected by the extension on current vertex,
        // we save this neighbor vertex and extend it later
        if (neighborQubitVertex != currentVertex)
        {
          affectedVertices.insert(neighborQubitVertex);
        }
      }
    }
    // If the neighbor qubit is not assigned to any vertex model yet,
    // it is added to the 'vertexToAvailableEdgesMapping' for the current vertex model.
    else
    {
      set<int> &currentVertexAvailableEdges =
          vertexToAvailableEdgesMapping.at(currentVertex);
      currentVertexAvailableEdges.insert(neighborQubit);
    }
  }
}

/* 
    Extend a partial embedding to satisfy the ratio.
*/
void extendChain(Graph &G, Graph &H, vector<set<int>> &vertexToQubitsMapping,
                 vector<vector<int>> &qubitToVertexMapping,
                 vector<set<int>> &vertexToAvailableEdgesMapping,
                 int currentVertex, float ratio, std::set<int> &affectedVertices)
{
  while (true)
  {
    // Calculate unembedded neighbors of the current vertex.
    int unembeddedNeighborVerticesNum = 0;
    vector<int> currentNeighborVertices = H.adj.at(currentVertex);
    for (int neighIndex = 0; neighIndex < currentNeighborVertices.size(); neighIndex++)
    {
      int neigh = currentNeighborVertices.at(neighIndex);
      if (vertexToQubitsMapping.at(neigh).size() > 0)
      {
        continue;
      }
      unembeddedNeighborVerticesNum++;
    }
    set<int> currentVertexAvailableEdges = vertexToAvailableEdgesMapping.at(currentVertex);
    int availableEdgesNum = currentVertexAvailableEdges.size();

    // If there are sufficient available edges, there is no need to extend current
    // vertex model.
    int avgChainLength = G.adj.size() / H.adj.size();
    if (availableEdgesNum >= std::ceil(ratio * unembeddedNeighborVerticesNum)
    || avgChainLength < vertexToQubitsMapping[currentVertex].size())
    {
      break;
    }

    // Get an available that adjacent to current vertex model randomly.
    int randomQubit = selectrandomFromSet(currentVertexAvailableEdges);
    if (randomQubit == -1)
    {
      return;
    }
    expand(qubitToVertexMapping, vertexToQubitsMapping,
           vertexToAvailableEdgesMapping, affectedVertices,
           randomQubit, currentVertex, G);
  }
}

void extendAll(Graph G, Graph H, float ratio,
               std::set<int> &affectedVertices, vector<vector<int>> &qubitToVertexMapping,
               vector<set<int>> &vertexToQubitsMapping,
               vector<set<int>> &vertexToAvailableEdgesMapping)
{
  while (affectedVertices.size() > 0)
  {
    int affectedVertex = *affectedVertices.begin();
    extendChain(G, H, vertexToQubitsMapping, qubitToVertexMapping,
                vertexToAvailableEdgesMapping, affectedVertex,
                ratio, affectedVertices);
    affectedVertices.erase(affectedVertex);
  }
}

vector<pair<int, int>> restorePath(int root, vector<int> previous,
                                   vector<tuple<int, int, int>> vertexTuples)
{
  vector<pair<int, int>> path;
  for (vector<tuple<int, int, int>>::iterator vertexTuplesIt = vertexTuples.begin();
       vertexTuplesIt != vertexTuples.end(); ++vertexTuplesIt)
  {
    int qubit = get<1>(*vertexTuplesIt);
    bool ignore = true; // ignore the first qubit as it is in neighbor vertex.
    while (true)
    {
      if (!ignore)
      {
        path.push_back(make_pair(qubit, get<0>(*vertexTuplesIt)));
      }
      ignore = false;
      if (qubit == root)
      {
        break;
      }
      qubit = previous[qubit];
    }
  }
  return path;
}

//find the union paths from root g* to each vertex model
void updateModels(Graph &G, Graph &H, int root, int currentVertex,
                  vector<int> previous,
                  vector<tuple<int, int, int>> vertexTuples,
                  vector<vector<int>> &qubitToVertexMapping,
                  vector<set<int>> &vertexToQubitsMapping,
                  vector<set<int>> &vertexToAvailableEdges,
                  float ratio)
{
  map<int, pair<int, int>> vertex_usage;
  vector<vector<int>> paths;
  vector<pair<int, int>> vertexToDistance;
  // qubit, leaf of path
  vector<pair<int, int>> path = restorePath(root, previous, vertexTuples);

  for (vector<pair<int, int>>::iterator pathIt = path.begin();
       pathIt != path.end(); pathIt++)
  { // don't want include the dummy vertex and the vertex in neighbor's map
    if (vertex_usage.find((*pathIt).first) == vertex_usage.end())
    {
      vertex_usage.insert(
          pair<int, pair<int, int>>((*pathIt).first,
                                    pair<int, int>(1, (*pathIt).second))); // if not exist, create a new one
    }
    else
    {
      vertex_usage[(*pathIt).first].first++; // update the usage
    }
  }

  std::set<int> affectedVertices;
  //update all vertex models

  for (map<int, pair<int, int>>::iterator ii = vertex_usage.begin();
       ii != vertex_usage.end(); ++ii)
  {
    if ((*ii).first == root)
    {
      expand(qubitToVertexMapping, vertexToQubitsMapping,
             vertexToAvailableEdges, affectedVertices, (*ii).first,
             currentVertex, G);
      continue;
    }
    if ((*ii).second.first > 1)
    { // appears in more than one path, assign to new model ****************************
      expand(qubitToVertexMapping, vertexToQubitsMapping,
             vertexToAvailableEdges, affectedVertices, (*ii).first,
             currentVertex, G);
    }

    else
    { // appears once, update

      expand(qubitToVertexMapping, vertexToQubitsMapping,
             vertexToAvailableEdges, affectedVertices, (*ii).first,
             (*ii).second.second, G);
    }
  }
  extendAll(G, H, ratio, affectedVertices, qubitToVertexMapping,
            vertexToQubitsMapping, vertexToAvailableEdges);
}

void tearChain(Graph &G, Graph &H, vector<set<int>> &vertexToQubitsMapping,
               vector<vector<int>> &qubitToVertexMapping,
               vector<set<int>> &vertexToAvailableEdgesMapping,
               int currentVertex, float ratio)
{
  set<int> currentVertexQubits = vertexToQubitsMapping.at(currentVertex);
  for (std::set<int>::iterator currentVertexQubitsIt = currentVertexQubits.begin();
       currentVertexQubitsIt != currentVertexQubits.end(); currentVertexQubitsIt++)
  {
    int currentVertexQubit = *currentVertexQubitsIt;

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
            set<int> &neighborQubitVertexAvailableEdges =
                vertexToAvailableEdgesMapping.at(neighborQubitVertex);
            neighborQubitVertexAvailableEdges.insert(currentVertexQubit);
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
  vertexToQubitsMapping.at(currentVertex).clear();
  vertexToAvailableEdgesMapping.at(currentVertex).clear();
}

//find the minimal vertex model
void findMinimalVertexModel(Graph &G, Graph &H, vector<set<int>> &vertexToQubitsMapping,
                            vector<vector<int>> &qubitToVertexMapping,
                            vector<set<int>> &vertexToAvailableEdgesMapping,
                            int currentVertex, float ratio)
{
  // if all neighbors are empty
  if (checkAllNeighbors(H, vertexToQubitsMapping, currentVertex))
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
        expand(qubitToVertexMapping, vertexToQubitsMapping, vertexToAvailableEdgesMapping,
               affectedVertices, qubit, currentVertex, G);

        // After randomly picking one qubit for this mapping,
        // we extend the mapping to maintain the ratio.
        extendAll(G, H, ratio, affectedVertices, qubitToVertexMapping,
                  vertexToQubitsMapping, vertexToAvailableEdgesMapping);
        return;
      }
    }
  }

  // tear the previous embedding for current vertex.
  tearChain(G, H, vertexToQubitsMapping, qubitToVertexMapping, vertexToAvailableEdgesMapping, currentVertex, ratio);

  // new improvement:  choose Candidate from first level
  vector<int> rootCandidates;

  // vector<cost, root, previous, vertextuple
  vector<tuple<int, int, vector<int>, vector<std::tuple<int, int, int>>>> rootCandidateDetails;
  for (int neighborVertexIndex = 0;
       neighborVertexIndex < H.deg(currentVertex);
       neighborVertexIndex++)
  {
    int neighborVertex = H.adj[currentVertex].at(neighborVertexIndex);
    set<int> neighborVertexQubits = vertexToQubitsMapping.at(neighborVertex);
    for (std::set<int>::iterator neighborVertexQubitsIt = neighborVertexQubits.begin();
         neighborVertexQubitsIt != neighborVertexQubits.end(); neighborVertexQubitsIt++)
    {
      int neighborVertexQubit = *neighborVertexQubitsIt;
      vector<int> firstLevelQubits = G.adj.at(neighborVertexQubit);
      for (int potentialCandidateIndex = 0;
           potentialCandidateIndex < firstLevelQubits.size();
           potentialCandidateIndex++)
      {
        int potentialCandidate = firstLevelQubits.at(potentialCandidateIndex);
        // ------------ Could be improved. -----------------
        if (find(rootCandidates.begin(),
                 rootCandidates.end(),
                 potentialCandidate) == rootCandidates.end() &&
            qubitToVertexMapping.at(potentialCandidate).size() == 0) // ?? Should I restrict this?
        {                                                            // if not contain
          rootCandidates.push_back(potentialCandidate);
        }
      }
    }
  }

  int minCost = 1 << 30;
  int root = 0;
  for (int i = 0; i < rootCandidates.size(); i++)
  {
    int currentRootCandidate = rootCandidates.at(i);
    vector<int> previous;
    vector<tuple<int, int, int>> vertexTuples;
    int sumCost = dijkstraComputePaths(currentRootCandidate, currentVertex, G, H,
                                       previous, vertexTuples,
                                       qubitToVertexMapping,
                                       vertexToQubitsMapping);
    rootCandidateDetails.push_back(make_tuple(
        sumCost, currentRootCandidate, previous, vertexTuples));
    // find root g* with minimum cost
    if (sumCost <= minCost)
    {
      minCost = sumCost;
    }
  }

  // random get the root from the first level rootCandidates
  vector<tuple<int, int, NLIST, vector<tuple<int, int, int>>>> candidates;
  for (int i = 0; i < rootCandidates.size(); i++)
  {
    if (minCost == get<0>(rootCandidateDetails.at(i)))
      candidates.push_back(rootCandidateDetails.at(i));
  }

  vector<int> previous;
  vector<tuple<int, int, int>> vertexTuples;
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
        vector<int> previous;
        vector<tuple<int, int, int>> vertexTuples;
        dijkstraComputePaths(root, currentVertex, G, H,
                             previous, vertexTuples, qubitToVertexMapping,
                             vertexToQubitsMapping);
        break;
      }
    }
  }
  else
  {
    random_shuffle(candidates.begin(), candidates.end()); // using built-in random generator:
    root = get<1>(candidates.at(0));
    previous = get<2>(candidates.at(0));
    vertexTuples = get<3>(candidates.at(0));
  }

  //find the union paths from root g* to each vertex model
  updateModels(G, H, root, currentVertex, previous, vertexTuples, qubitToVertexMapping,
               vertexToQubitsMapping, vertexToAvailableEdgesMapping, ratio);
}

// main method for finding Minor embedding
bool findMinorEmbedding(Graph &G, Graph &H, float ratio)
{
  //srand(unsigned(time(0))); // for randomizing
#if 1
  vector<int> vertexOrder = randomizeVertexOrder(H.order());
#else
  vector<int> vertexOrder = sortVertexOrder(H); // sort order by size
#endif
  vector<set<int>> vertexToQubitsMapping = initializeSetMapping(H.order());
  vector<vector<int>> qubitsToVerticesMapping = initializeMapping(G.order());
  vector<set<int>> vertexToAvailableEdgesMapping = initializeSetMapping(H.order());
  int previousChainLength;
  int previousOverlap;
  int currentChainLength;
  int currentOverlap;
  int stage = 1;
  while (true)
  {
    currentChainLength = stage <= 1 ? INT_MAX : getTotalSetSize(vertexToQubitsMapping);
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
      findMinimalVertexModel(G, H, vertexToQubitsMapping, qubitsToVerticesMapping,
                             vertexToAvailableEdgesMapping, currentVertex, ratio);
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
    for (int i = 0; i < vertexToQubitsMapping.size(); i++)
    {
      cout << "[";
      totalQubit = totalQubit + vertexToQubitsMapping.at(i).size();
      if (vertexToQubitsMapping.at(i).size() > max)
        max = vertexToQubitsMapping.at(i).size();
      for (std::set<int>::iterator it = vertexToQubitsMapping.at(i).begin();
           it != vertexToQubitsMapping.at(i).end();
           ++it)
      {
        cout << *it << ", ";
      }
      if (i == vertexToQubitsMapping.size() - 1)
      {
        cout << "\b\b]";
      }
      else
      {
        cout << "\b\b], ";
      }
    }
    cout << "]";
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
  cout << "Running Time : " << (double)(ends - start) / CLOCKS_PER_SEC << endl;
  return 0;
}
