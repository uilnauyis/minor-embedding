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

enum RootFindingTechnique
{
  Bfs,
  Dijkstra,
  Hybrid,
  Last
};

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

int getOverlapNum(vector<vector<int>> qubitsToVerticesMapping)
{
  int count = 0;
  for (int i = 0; i < qubitsToVerticesMapping.size(); i++)
  {
    int thisOverlap = qubitsToVerticesMapping.at(i).size();
    count = thisOverlap > 1 ? count + 1 : count;
  }
  return count;
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
int getOverlaps(vector<vector<int>> qubitsToVerticesMapping)
{
  int overlapNum = 0;
  for (int i = 0; i < qubitsToVerticesMapping.size(); i++)
  {
    int verticesNum = qubitsToVerticesMapping.at(i).size();
    overlapNum += verticesNum > 0 ? verticesNum - 1 : 0;
  }
  return overlapNum;
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
int getQubitWeight(int mappingNum)
{
  int DIAMETER = 30;
  return pow(DIAMETER, mappingNum - 1);
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

int dijkstraFindRoot(int currentVertex, Graph G, Graph H,
                     vector<vector<int>> &qubitToVertexMapping,
                     vector<set<int>> &vertexToQubitsMapping)
{
  int root = 0;
  int rootMinDistance = INT_MAX;
  int n = G.order();
  vector<int> distance;
  distance.resize(n, 0);

  for (int index = 0; index < G.order(); index++)
  {
    if (G.adj[index].size() == 0)
    {
      distance.at(index) = INT_MAX;
    }
  }

  for (vector<int>::iterator vertexIndex = H.adj.at(currentVertex).begin();
       vertexIndex != H.adj.at(currentVertex).end(); vertexIndex++)
  //for (int vertex = 0; vertex < H.order(); vertex++)
  {
    int vertex = *vertexIndex;
    set<int> qubits = vertexToQubitsMapping.at(vertex);
    set<pair<int, int>> pq; // first: distance, second: qubit

    vector<int> min_distance;
    min_distance.resize(n, -1);

    for (set<int>::iterator qubitsIt = qubits.begin();
         qubitsIt != qubits.end(); qubitsIt++)
    {
      int qubit = *qubitsIt;

      for (vector<int>::iterator neighborIt = G.adj.at(qubit).begin();
           neighborIt != G.adj.at(qubit).end(); neighborIt++)
      {
        int neighborQubit = *neighborIt;
        if (find(qubits.begin(), qubits.end(), neighborQubit) == qubits.end())
        {
          int overlapCount = qubitToVertexMapping.at(neighborQubit).size();
          int weight = getQubitWeight(qubitToVertexMapping.at(neighborQubit).size() + 1);
          min_distance[neighborQubit] = weight;
          pq.insert(pair<int, int>(min_distance[neighborQubit], neighborQubit));
        }
      }
    }

    while (pq.size() > 0)
    {
      int dist = pq.begin()->first;
      int qubit = pq.begin()->second;
      pq.erase(pq.begin());

      for (int i = 0; i < G.deg(qubit); i++)
      {
        int v = G.adj[qubit].at(i);
        int overlapCount = qubitToVertexMapping.at(v).size();
        int weight = getQubitWeight(overlapCount + 1);
        int newDistance = dist + weight;

        if (min_distance[v] != -1 && min_distance[v] <= newDistance)
        {
          continue;
        }

        pq.erase(make_pair(min_distance[v], v));
        min_distance[v] = newDistance;
        pq.insert(pair<int, int>(min_distance[v], v));
      }
    }

    for (int i = 0; i < min_distance.size(); i++)
    {
      distance[i] += min_distance[i] == -1 ? 0 : min_distance[i];
    }
  }

  for (int i = 0; i < distance.size(); i++)
  {
    if (distance[i] < rootMinDistance)
    {
      rootMinDistance = distance[i];
      root = i;
    }
  }
  return root;
}

int bfsFindRoot(int currentVertex, Graph G, Graph H,
                vector<vector<int>> &qubitToVertexMapping,
                vector<set<int>> &vertexToQubitsMapping)
{
  int root = 0;
  int rootMinDistance = INT_MAX;
  int n = G.order();
  vector<int> distance;
  distance.resize(n, 0);

  for (int index = 0; index < G.order(); index++)
  {
    if (G.adj[index].size() == 0)
    {
      distance.at(index) = INT_MAX;
    }
  }

  for (vector<int>::iterator vertexIndex = H.adj.at(currentVertex).begin();
       vertexIndex != H.adj.at(currentVertex).end(); vertexIndex++)
  //for (int vertex = 0; vertex < H.order(); vertex++)
  {
    int vertex = *vertexIndex;
    set<int> qubits = vertexToQubitsMapping.at(vertex);
    vector<queue<pair<int, int>>> queues; // first: distance, second: qubit

    vector<int> min_distance;
    min_distance.resize(n, -1);

    for (set<int>::iterator qubitsIt = qubits.begin();
         qubitsIt != qubits.end(); qubitsIt++)
    {
      int qubit = *qubitsIt;

      for (vector<int>::iterator neighborIt = G.adj.at(qubit).begin();
           neighborIt != G.adj.at(qubit).end(); neighborIt++)
      {
        int neighborQubit = *neighborIt;
        if (find(qubits.begin(), qubits.end(), neighborQubit) == qubits.end())
        {
          int overlapCount = qubitToVertexMapping.at(neighborQubit).size();
          int weight = getQubitWeight(qubitToVertexMapping.at(neighborQubit).size() + 1);
          min_distance[neighborQubit] = weight;
          while (overlapCount + 1 > queues.size())
          {
            queue<pair<int, int>> queue;
            queues.push_back(queue);
          }
          queues.at(overlapCount).push(pair<int, int>(min_distance[neighborQubit], neighborQubit));
        }
      }
    }

    while (true)
    {
      int dist = -1, qubit = -1, poppedQueueIndex = -1;
      for (int i = 0; i < queues.size(); i++)
      {
        queue<pair<int, int>> queue = queues.at(i);
        if (queue.size() > 0)
        {
          pair<int, int> front = queue.front();
          if (dist == -1 || dist > front.first)
          {
            dist = front.first;
            qubit = front.second;
            poppedQueueIndex = i;
          }
        }
      }

      if (dist == -1)
      { // all queues are empty.
        break;
      }
      queue<pair<int, int>> &thisQueue = queues.at(poppedQueueIndex);
      thisQueue.pop();

      //
      //int distance_through_u = dist + weight;
      // Visit each edge exiting u
      for (int i = 0; i < G.deg(qubit); i++)
      {
        int v = G.adj[qubit].at(i);
        int overlapCount = qubitToVertexMapping.at(v).size();
        int weight = getQubitWeight(overlapCount + 1);

        // If v has been visited or has been mapped, ignore
        if (min_distance[v] != -1)
        {
          continue;
        }

        min_distance[v] = dist + weight;
        while (overlapCount + 1 > queues.size())
        {
          queue<pair<int, int>> queue;
          queues.push_back(queue);
        }
        queues.at(overlapCount).push(pair<int, int>(min_distance[v], v));
      }
    }

    for (int i = 0; i < min_distance.size(); i++)
    {
      distance[i] += min_distance[i] == -1 ? 0 : min_distance[i];
    }
  }

  for (int i = 0; i < distance.size(); i++)
  {
    if (distance[i] < rootMinDistance)
    {
      rootMinDistance = distance[i];
      root = i;
    }
  }
  return root;
}

int hybridFindRoot(int currentVertex, Graph G, Graph H,
                   vector<vector<int>> &qubitToVertexMapping,
                   vector<set<int>> &vertexToQubitsMapping)
{
  int root = 0;
  int rootMinDistance = INT_MAX;
  int n = G.order();
  vector<int> distance;
  distance.resize(n, 0);

  for (int index = 0; index < G.order(); index++)
  {
    if (G.adj[index].size() == 0)
    {
      distance.at(index) = INT_MAX;
    }
  }

  for (vector<int>::iterator vertexIndex = H.adj.at(currentVertex).begin();
       vertexIndex != H.adj.at(currentVertex).end(); vertexIndex++)
  {
    int vertex = *vertexIndex;
    set<int> qubits = vertexToQubitsMapping.at(vertex);
    queue<pair<int, int>> q; // first: distance, second: qubit
    set<pair<int, int>> pq;   // first: distance, second: qubit

    vector<int> min_distance;
    min_distance.resize(n, -1);

    for (set<int>::iterator qubitsIt = qubits.begin();
         qubitsIt != qubits.end(); qubitsIt++)
    {
      int qubit = *qubitsIt;

      for (vector<int>::iterator neighborIt = G.adj.at(qubit).begin();
           neighborIt != G.adj.at(qubit).end(); neighborIt++)
      {
        int neighborQubit = *neighborIt;
        if (find(qubits.begin(), qubits.end(), neighborQubit) == qubits.end())
        {
          int overlapCount = qubitToVertexMapping.at(neighborQubit).size();
          int weight = getQubitWeight(qubitToVertexMapping.at(neighborQubit).size() + 1);
          min_distance.at(neighborQubit) = weight;
          if (weight == 1)
            q.push(pair<int, int>(1, neighborQubit));
          else
          {
            pq.insert(pair<int, int>(min_distance[neighborQubit], neighborQubit));
          }
        }
      }
    }

    while (q.size() != 0 || pq.size() != 0)
    {
      int dist, qubit;
      if (q.size() == 0 && pq.size() != 0)
      {
        dist = pq.begin()->first;
        qubit = pq.begin()->second;
        pq.erase(pq.begin());
      }
      else if (pq.size() == 0 && q.size() != 0)
      {
        pair<int, int> front = q.front();
        dist = front.first;
        qubit = front.second;
        q.pop();
      }
      else
      {
        pair<int, int> front = q.front();
        if (front.first < pq.begin()->first)
        {
          dist = front.first;
          qubit = front.second;
          q.pop();
        }
        else
        {
          dist = pq.begin()->first;
          qubit = pq.begin()->second;
          pq.erase(pq.begin());
        }
      }

      for (int i = 0; i < G.deg(qubit); i++)
      {
        int v = G.adj[qubit].at(i);
        int overlapCount = qubitToVertexMapping.at(v).size();
        int weight = getQubitWeight(overlapCount + 1);
        int newDistance = dist + weight;

        if (min_distance[v] != -1)
        {
          continue;
        }

        min_distance[v] = newDistance;
        if (weight == 1)
        {
          q.push(pair<int, int>(min_distance[v], v));
        }
        else
        {
          pq.insert(pair<int, int>(min_distance[v], v));
        }
      }
    }

    for (int i = 0; i < min_distance.size(); i++)
    {
      distance[i] += min_distance[i] == -1 ? 0 : min_distance[i];
    }
  }

  for (int i = 0; i < distance.size(); i++)
  {
    if (distance[i] < rootMinDistance)
    {
      rootMinDistance = distance[i];
      root = i;
    }
  }
  return root;
}

// Compute shortest paths from current vertex to all neighbor vertices using
// multiple FIFO queueS.
// Return the total cost as result.
// 'previous' is the vector for tracing back the path.
int bfsComputePaths(int source, int currentVertex,
                    Graph &G, Graph &H,
                    vector<int> &previous,
                    vector<tuple<int, int, int>> &vertexTuples, // <vertex, qubit, distance>
                    vector<vector<int>> &qubitToVertexMapping,
                    vector<set<int>> &vertexToQubitsMapping)
{
  int n = G.order();
  vector<int> min_distance;
  min_distance.resize(n, -1); // using '-1' to distinguish unvisited qubits
  min_distance[source] = 0;
  previous.clear();
  vertexTuples.clear();
  previous.resize(n, -1);
  vector<pair<int, int>> q; // first: distance, second: qubit
  set<pair<int, int>> pq;   // first: distance, second: qubit

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

  q.push_back(pair<int, int>(0, source));
  while (q.size() != 0 || pq.size() != 0)
  {
    int dist, qubit;
    if (q.size() == 0 && pq.size() != 0)
    {
      dist = pq.begin()->first;
      qubit = pq.begin()->second;
      pq.erase(pq.begin());
    }
    else if (pq.size() == 0 && q.size() != 0)
    {
      dist = q.begin()->first;
      qubit = q.begin()->second;
      q.erase(q.begin());
    }
    else
    {
      if (q.begin()->first < pq.begin()->first)
      {
        dist = q.begin()->first;
        qubit = q.begin()->second;
        q.erase(q.begin());
      }
      else
      {
        dist = pq.begin()->first;
        qubit = pq.begin()->second;
        pq.erase(pq.begin());
      }
    }

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
    for (int i = 0; i < G.deg(qubit); i++)
    {
      int v = G.adj[qubit].at(i);

      // If v has been visited, ignore
      if (min_distance[v] != -1)
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
      min_distance[v] = distance_through_u;
      previous[v] = qubit;
      if (weight == 1)
      {
        q.push_back(pair<int, int>(distance_through_u, v));
      }
      else
      {
        pq.insert(pair<int, int>(distance_through_u, v));
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
            int currentQubit,
            int currentVertex,
            Graph G)
{
  // Add random qubit to the current vertex model.
  qubitToVertexMapping.at(currentQubit).push_back(currentVertex);
  vertexToQubitsMapping.at(currentVertex).insert(currentQubit);
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

  for (map<int, pair<int, int>>::iterator ii = vertex_usage.begin();
       ii != vertex_usage.end(); ++ii)
  {
    if ((*ii).first == root || (*ii).second.first > 1)
    { // appears in more than one path, assign to new model ****************************
      expand(qubitToVertexMapping, vertexToQubitsMapping,
             (*ii).first, currentVertex, G);
    }

    else
    { // appears once, update

      expand(qubitToVertexMapping, vertexToQubitsMapping,
             (*ii).first, (*ii).second.second, G);
    }
  }
}

void tearChain(Graph &G, Graph &H, vector<set<int>> &vertexToQubitsMapping,
               vector<vector<int>> &qubitToVertexMapping,
               int currentVertex, float ratio)
{
  set<int> currentVertexQubits = vertexToQubitsMapping.at(currentVertex);
  for (std::set<int>::iterator currentVertexQubitsIt = currentVertexQubits.begin();
       currentVertexQubitsIt != currentVertexQubits.end(); currentVertexQubitsIt++)
  {
    int currentVertexQubit = *currentVertexQubitsIt;
    vector<int> &currentMappedVertices = qubitToVertexMapping.at(currentVertexQubit);
    currentMappedVertices.erase(
        std::remove(currentMappedVertices.begin(),
                    currentMappedVertices.end(), currentVertex),
        currentMappedVertices.end());
  }
  vertexToQubitsMapping.at(currentVertex).clear();
}

//find the minimal vertex model
void findMinimalVertexModel(Graph &G, Graph &H, vector<set<int>> &vertexToQubitsMapping,
                            vector<vector<int>> &qubitToVertexMapping,
                            int currentVertex, float ratio, int rootFindingTechnique)
{
  // tear the previous embedding for current vertex.
  tearChain(G, H, vertexToQubitsMapping, qubitToVertexMapping, currentVertex, ratio);
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
        expand(qubitToVertexMapping, vertexToQubitsMapping,
               qubit, currentVertex, G);
        return;
      }
    }
  }

  int root;
  if (rootFindingTechnique == Dijkstra)
  {
    root = dijkstraFindRoot(currentVertex, G, H, qubitToVertexMapping, vertexToQubitsMapping);
  }
  else if (rootFindingTechnique == Bfs)
  {
    root = bfsFindRoot(currentVertex, G, H, qubitToVertexMapping, vertexToQubitsMapping);
  }
  else
  {
    root = hybridFindRoot(currentVertex, G, H, qubitToVertexMapping, vertexToQubitsMapping);
  }
  //if (rootD != rootH || rootD != rootB || rootB != rootH)
  //{
  //  throw runtime_error("different root found");
  //}

  vector<int> previous;
  vector<tuple<int, int, int>> vertexTuples;
  bfsComputePaths(root, currentVertex, G, H, previous, vertexTuples,
                  qubitToVertexMapping, vertexToQubitsMapping);

  //find the union paths from root g* to each vertex model
  updateModels(G, H, root, currentVertex, previous, vertexTuples, qubitToVertexMapping,
               vertexToQubitsMapping, ratio);
}

// main method for finding Minor embedding
bool findMinorEmbedding(Graph &G, Graph &H, float ratio, int &qubitNumSum, int &maxChainSum, int rootFindingTechnique)
{
  srand(unsigned(time(0))); // for randomizing
#if 1
  vector<int> vertexOrder = randomizeVertexOrder(H.order());
#else
  vector<int> vertexOrder = sortVertexOrder(H); // sort order by size
#endif
  vector<set<int>> vertexToQubitsMapping = initializeSetMapping(H.order());
  vector<vector<int>> qubitsToVerticesMapping = initializeMapping(G.order());
  vector<set<int>> vertexToAvailableEdgesMapping = initializeSetMapping(H.order());
  vector<set<int>> previousVertexToQubitsMapping;
  int previousChainLength = INT_MAX;
  int previousOverlap = INT_MAX;
  int previousOverlapNum = INT_MAX;
  int currentChainLength = INT_MAX;
  int currentOverlap = INT_MAX;
  int currentOverlapNum = INT_MAX;
  int stage = 1;
  while (true)
  {
    if (previousChainLength <= currentChainLength &&
        (previousOverlap <= currentOverlap) &&
        (previousOverlapNum <= currentOverlapNum) &&
        stage > 2) // No improvement in last update.
    {
      break;
    }

    previousVertexToQubitsMapping = vertexToQubitsMapping;

    //std::cout << "Current chain length: " << currentChainLength
    //          << " current overlap: " << currentOverlap
    //          << " current overlap num: " << currentOverlapNum << endl;

    previousChainLength = currentChainLength;
    previousOverlap = currentOverlap;
    previousOverlapNum = currentOverlapNum;
    for (int i = 0; i < vertexOrder.size(); i++)
    {
      int currentVertex = vertexOrder.at(i);
      findMinimalVertexModel(G, H, vertexToQubitsMapping, qubitsToVerticesMapping,
                             currentVertex, ratio, rootFindingTechnique);
    }
    currentChainLength = stage <= 1 ? INT_MAX : getTotalSetSize(vertexToQubitsMapping);
    currentOverlap = stage <= 1 ? INT_MAX : getMaxOverlap(qubitsToVerticesMapping);
    currentOverlapNum = getOverlapNum(qubitsToVerticesMapping);
    stage++;
  }

  //check the usage for each vertex in G, if all <=1, then we find the embedding. otherwise "failure"
  bool embedding = previousOverlap <= 1;
  //print the embedding, if find
  if (embedding)
  {
    int totalQubit = 0;
    int max = 0;
    ofstream writeFile;
    writeFile.open("./alists/output.alist", std::ios_base::app);
    //std::cout << "[";
    writeFile << "[";
    for (int i = 0; i < previousVertexToQubitsMapping.size(); i++)
    {
      //std::cout << "[";
      writeFile << "[";

      totalQubit = totalQubit + previousVertexToQubitsMapping.at(i).size();
      if (previousVertexToQubitsMapping.at(i).size() > max)
        max = previousVertexToQubitsMapping.at(i).size();
      for (std::set<int>::iterator it = previousVertexToQubitsMapping.at(i).begin();
           it != previousVertexToQubitsMapping.at(i).end();
           ++it)
      {
        if (it == previousVertexToQubitsMapping.at(i).begin())
        {
          //std::cout << *it;
          writeFile << *it;
        }
        else
        {
          //std::cout << ", " << *it;
          writeFile << ", " << *it;
        }
      }
      if (i == previousVertexToQubitsMapping.size() - 1)
      {
        //std::cout << "]";
        writeFile << "]";
      }
      else
      {
        //std::cout << "], ";
        writeFile << "], ";
      }
    }
    //std::cout << "]" << endl;
    writeFile << "]" << endl;
    //std::cout << "Total number of Qubits used is: " << totalQubit << endl;
    //std::cout << "Max chain length is: " << max << endl;
    writeFile.close();

    qubitNumSum += totalQubit;
    maxChainSum += max;

    return true;
  }
  else
  {
    return false;
  }
}

int main(int argc, char **argv)
{
  clock_t successRunTime = 0;

  ifstream host(argv[2]);
  Graph G(0);
  host >> G;
  Graph H(0);
  std::cin >> H;
  float ratio = std::stof(argv[1]);
  bool foundEmbedding;
  std::cout << "Guest graph: " << endl;
  std::cout << H << endl;
  for (int rootFindingTechnique = Bfs; rootFindingTechnique < Last; rootFindingTechnique++)
  {
    int qubitNumSum = 0;
    int maxChainSum = 0;
    int successCount = 0;
    int successRunTime = 0;
    int ATTEMPSNUM = 100;
    for (int i = 1; i <= ATTEMPSNUM; i++)
    {
      clock_t start = clock();
      std::cout << "try " << i << ":" << endl;
      foundEmbedding = findMinorEmbedding(G, H, ratio, qubitNumSum, maxChainSum, rootFindingTechnique);
      clock_t end = clock();
      if (foundEmbedding)
      {
        successRunTime += end - start;
        successCount += 1;
      }
    }

    std::stringstream sstm;
    sstm << "./meow/" << rootFindingTechnique;

    ofstream writeFile(sstm.str(), std::ios_base::app);
    writeFile << argv[3] << " "
              << (double)successRunTime / CLOCKS_PER_SEC / successCount << ", "
              << (double)qubitNumSum / successCount << ", "
              << (double)maxChainSum / successCount << ", "
              << (double)successCount / ATTEMPSNUM << "\n";
    writeFile.close();
    cout << argv[3];
  }
}
