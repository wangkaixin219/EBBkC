#!/usr/bin/zsh

#./BBkC e /home/kaixin/Dataset/small/baidu.edges >> small.txt
#./BBkC e /home/kaixin/Dataset/small/facebook-wosn.edges >> small.txt
#./BBkC e /home/kaixin/Dataset/small/nasasrb.edges >> small.txt
#./BBkC e /home/kaixin/Dataset/small/sc-shipsec5.edges >> small.txt
#./BBkC e /home/kaixin/Dataset/small/soc-pokec.edges >> small.txt
#./BBkC e /home/kaixin/Dataset/small/wiki.edges >> small.txt
#
#./BBkC e /home/kaixin/Dataset/large/citeseer.edges >> small.txt
#./BBkC e /home/kaixin/Dataset/large/dblp.edges >> small.txt
#./BBkC e /home/kaixin/Dataset/large/digg.edges >> small.txt
#./BBkC e /home/kaixin/Dataset/large/sk.edges >> small.txt
#./BBkC e /home/kaixin/Dataset/large/skitter.edges >> small.txt
#./BBkC e /home/kaixin/Dataset/large/stanford.edges >> small.txt



for k in {3..31}
do
    echo "k = $k"
    echo "k = $k" >> small.txt
	for l in {0..5}
	do
		echo "l = $l"
		./BBkC e /home/kaixin/Dataset/small/baidu.index "$k" "$l" 2 >> small.txt
	done
	./BBkC v /home/kaixin/Dataset/small/baidu.index "$k" 0 0 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/baidu.index "$k" 0 1 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/baidu.index "$k" 0 2 >> small.txt
	./BBkC v /home/kaixin/Dataset/small/baidu.index "$k" 0 3 >> small.txt
    echo
done


for k in {3..30}
do
    echo "k = $k"
    echo "k = $k" >> small.txt
	for l in {0..5}
    do
		echo "l = $l"
        ./BBkC e /home/kaixin/Dataset/small/facebook.index "$k" "$l" 2 >> small.txt
    done
    ./BBkC v /home/kaixin/Dataset/small/facebook.index "$k" 0 0 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/facebook.index "$k" 0 1 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/facebook.index "$k" 0 2 >> small.txt
	./BBkC v /home/kaixin/Dataset/small/facebook.index "$k" 0 3 >> small.txt
	echo
done


for k in {3..24}
do
    echo "k = $k"
    echo "k = $k" >> small.txt
	for l in {0..5}
    do
		echo "l = $l"
        ./BBkC e /home/kaixin/Dataset/small/nasasrb.index "$k" "$l" 2 >> small.txt
    done
    ./BBkC v /home/kaixin/Dataset/small/nasasrb.index "$k" 0 0 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/nasasrb.index "$k" 0 1 >> small.txt
	./BBkC v /home/kaixin/Dataset/small/nasasrb.index "$k" 0 2 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/nasasrb.index "$k" 0 3 >> small.txt 
   	echo
done


for k in {3..24}
do
    echo "k = $k"
    echo "k = $k" >> small.txt
	for l in {0..5}
    do
		echo "l = $l"
        ./BBkC e /home/kaixin/Dataset/small/shipsec.index "$k" "$l" 2 >> small.txt
    done
    ./BBkC v /home/kaixin/Dataset/small/shipsec.index "$k" 0 0 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/shipsec.index "$k" 0 1 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/shipsec.index "$k" 0 2 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/shipsec.index "$k" 0 3 >> small.txt
	echo
done


for k in {3..25}
do
    echo "k = $k"
    echo "k = $k" >> small.txt
	for l in {0..5}
    do
		echo "l = $l"
        ./BBkC e /home/kaixin/Dataset/small/wiki.index "$k" "$l" 2 >> small.txt
    done
    ./BBkC v /home/kaixin/Dataset/small/wiki.index "$k" 0 0 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/wiki.index "$k" 0 1 >> small.txt
	./BBkC v /home/kaixin/Dataset/small/wiki.index "$k" 0 2 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/wiki.index "$k" 0 3 >> small.txt
    echo
done


for k in {3..29}
do
    echo "k = $k"
    echo "k = $k" >> small.txt
	for l in {0..5}
    do
		echo "l = $l"
        ./BBkC e /home/kaixin/Dataset/small/pokec.index "$k" "$l" 2 >> small.txt
    done
    ./BBkC v /home/kaixin/Dataset/small/pokec.index "$k" 0 0 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/pokec.index "$k" 0 1 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/pokec.index "$k" 0 2 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/pokec.index "$k" 0 3 >> small.txt
	echo
done


for k in {3..33}
do
    echo "k = $k"
    echo "k = $k" >> small.txt
	for l in {0..5}
    do
		echo "l = $l"
        ./BBkC e /home/kaixin/Dataset/small/cn.index "$k" "$l" 2 >> small.txt
    done
    ./BBkC v /home/kaixin/Dataset/small/cn.index "$k" 0 0 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/cn.index "$k" 0 1 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/cn.index "$k" 0 2 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/cn.index "$k" 0 3 >> small.txt
    echo
done


for k in {3..17}
do
    echo "k = $k"
    echo "k = $k" >> small.txt
	for l in {0..5}
    do
		echo "l = $l"
        ./BBkC e /home/kaixin/Dataset/small/youtube.index "$k" "$l" 2 >> small.txt
    done
    ./BBkC v /home/kaixin/Dataset/small/youtube.index "$k" 0 0 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/youtube.index "$k" 0 1 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/youtube.index "$k" 0 2 >> small.txt
    ./BBkC v /home/kaixin/Dataset/small/youtube.index "$k" 0 3 >> small.txt
    echo
done
