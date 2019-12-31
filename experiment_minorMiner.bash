#!/bin/bash
mkdir -p ./meow-api
for ((i=30; i <=50; i++))
do
  python ./sapi2_embed.py ./alists/clique${i}.alist ./alists/dw2x.alist
done

