#!/usr/bin/python
# check minor embedding properties  (python 2.7; networkx 1.9; 11-July-2015)
# usage: check_embed.py guest.alist host.alist [embed.txt]

import sys, itertools
import networkx as nx

def read_graph(infile=sys.stdin):  # CS220 adjacency list format
    n=int(infile.readline().strip())
    G=nx.empty_graph(n,create_using=nx.Graph()) 
    for u in range(n):
        neighbors=infile.readline().split()
        for v in neighbors: G.add_edge(u,int(v))
    return G

def print_graph(G): # CS220 adjacency list format
    n=G.order()
    print n
    for u in range(n):
        for v in G[u]: print v,
        print

assert 3<=len(sys.argv)<=4
print sys.argv

gfile=open(sys.argv[1],'r')
Guest=read_graph(gfile)
hfile=open(sys.argv[2],'r')
Host=read_graph(hfile)

if len(sys.argv)==4: # assume embedding is at end of guest file, otherwise own file
    gfile.close()
    gfile=open(sys.argv[3],'r')
embedding=eval(gfile.readline())

gfile.close(); hfile.close()

nG=Guest.order()
assert nG == len(embedding)

#print_graph(Guest) #; print_graph(Host)
print 'embedding=', embedding 
k=max([len(embedding[i]) for i in range(nG)])
print 'max map size=', k
a=sum([len(embedding[i]) for i in range(nG)])/(nG+0.0)
print 'avg map size=', a

# check embedding disjoint
flat=[]
for i in range(nG):
  flat.extend(embedding[i])
for x in flat:
  if flat.count(x) >= 2:
    print x
assert len([x for x in flat if flat.count(x) >= 2])==0
print 'embedding is disjoint'

# check embedding connected pieces
for i in range(nG):
   SG=Host.subgraph(embedding[i])
   assert nx.is_connected(SG)
print 'embedding is connected pieces'

# check path or set of pieces
#isPath = False
#for i in range(nG):
#   isPath = False
#   l=len(embedding[i]) 
#   if l<4: 
#      isPath=True; continue
#   SG=Host.subgraph(embedding[i])
#   for p in itertools.permutations(embedding[i]):
#       curPath=True
#       for j in range(l-1):
#          if SG.has_edge(p[j],p[j+1])==False: 
#             curPath=False; break
#       if curPath:
#          isPath=True; break
#   if isPath==False: 
#      print 'non-path', embedding[i], SG.edges()
#      break
#           
#if isPath: print 'each set is a path'
#else:      print 'at least one set not a path'

# check guest adjacencies in host
isConnected=False
for (u,v) in Guest.edges():
   isConnected=False
   for (x,y) in itertools.product(embedding[u],embedding[v]):
       if Host.has_edge(x,y): 
          isConnected=True; break
   if isConnected==False: 
      print 'non-connected', u,v
      break

if isConnected: print 'is minor embedding'
else:           print 'is not minor embedding'
