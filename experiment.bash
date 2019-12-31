#!/bin/bash
mkdir -p ./meow
/usr/bin/g++ -O3 -g /home/siyuan/Repo/minor-embedding/src/embedding.cpp /home/siyuan/Repo/minor-embedding/src/Graph.h /home/siyuan/Repo/minor-embedding/src/Graph.cpp -o /home/siyuan/Repo/minor-embedding/meow/embedding
for ((i=45; i <=50; i++))
do
  ./meow/embedding 0.5 ./alists/dw2x.alist ${i} < ./alists/clique${i}.alist
done

