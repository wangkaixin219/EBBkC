#!/usr/bin/zsh

for k in {3..8}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
    ./BBkC e /home/kaixin/Dataset/large/citeseer.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/citeseer.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/citeseer.index "$k" 0 1 >> large.txt
	./BBkC v /home/kaixin/Dataset/large/citeseer.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/citeseer.index "$k" 0 3 >> large.txt
    echo
done


for k in {3..8}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
    ./BBkC e /home/kaixin/Dataset/large/dblp.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/dblp.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/dblp.index "$k" 0 1 >> large.txt
	./BBkC v /home/kaixin/Dataset/large/dblp.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/dblp.index "$k" 0 3 >> large.txt
    echo
done


for k in {3..8}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
    ./BBkC e /home/kaixin/Dataset/large/digg.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/digg.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/digg.index "$k" 0 1 >> large.txt
	./BBkC v /home/kaixin/Dataset/large/digg.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/digg.index "$k" 0 3 >> large.txt
    echo
done


for k in {3..8}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
    ./BBkC e /home/kaixin/Dataset/large/sk.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/sk.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/sk.index "$k" 0 1 >> large.txt
	./BBkC v /home/kaixin/Dataset/large/sk.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/sk.index "$k" 0 3 >> large.txt
    echo
done


for k in {3..8}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
    ./BBkC e /home/kaixin/Dataset/large/skitter.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/skitter.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/skitter.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/skitter.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/skitter.index "$k" 0 3 >> large.txt
	echo
done


for k in {3..8}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
    ./BBkC e /home/kaixin/Dataset/large/stanford.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/stanford.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/stanford.index "$k" 0 1 >> large.txt
	./BBkC v /home/kaixin/Dataset/large/stanford.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/stanford.index "$k" 0 3 >> large.txt
    echo
done


for k in {3..8}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
    ./BBkC e /home/kaixin/Dataset/large/diel.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/diel.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/diel.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/diel.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/diel.index "$k" 0 3 >> large.txt
    echo
done


for k in {3..8}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
    ./BBkC e /home/kaixin/Dataset/large/orkut.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/orkut.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/orkut.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/orkut.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/orkut.index "$k" 0 3 >> large.txt
    echo
done







for k in {82..87}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
	./BBkC e /home/kaixin/Dataset/large/citeseer.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/citeseer.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/citeseer.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/citeseer.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/citeseer.index "$k" 0 3 >> large.txt
    echo
done


for k in {109..114}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
	./BBkC e /home/kaixin/Dataset/large/dblp.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/dblp.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/dblp.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/dblp.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/dblp.index "$k" 0 3 >> large.txt
    echo
done


for k in {45..50}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
	./BBkC e /home/kaixin/Dataset/large/digg.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/digg.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/digg.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/digg.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/digg.index "$k" 0 3 >> large.txt
    echo
done


for k in {77..82}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
	./BBkC e /home/kaixin/Dataset/large/sk.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/sk.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/sk.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/sk.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/sk.index "$k" 0 3 >> large.txt
    echo
done


for k in {62..67}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
	./BBkC e /home/kaixin/Dataset/large/skitter.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/skitter.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/skitter.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/skitter.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/skitter.index "$k" 0 3 >> large.txt
    echo
done


for k in {56..61}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
	./BBkC e /home/kaixin/Dataset/large/stanford.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/stanford.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/stanford.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/stanford.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/stanford.index "$k" 0 3 >> large.txt
    echo
done


for k in {40..45}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
    ./BBkC e /home/kaixin/Dataset/large/diel.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/diel.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/diel.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/diel.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/diel.index "$k" 0 3 >> large.txt
    echo
done


for k in {42..47}
do
    echo "k = $k"
    echo "k = $k" >> large.txt
    ./BBkC e /home/kaixin/Dataset/large/orkut.index "$k" 2 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/orkut.index "$k" 0 0 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/orkut.index "$k" 0 1 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/orkut.index "$k" 0 2 >> large.txt
    ./BBkC v /home/kaixin/Dataset/large/orkut.index "$k" 0 3 >> large.txt
    echo
done
