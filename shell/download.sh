#!/bin/zsh

DATASET="/home/kaixin/DESCol/dataset/"
ZIPFILE="$DATASET/zipfile"
LARGE="$DATASET/large"
SMALL="$DATASET/small"
SYN="$DATASET/syn"

cd $DATASET

if [ ! -d "$ZIPFILE" ]; then
	mkdir zipfile
fi

if [ ! -d "$LARGE" ]; then
	mkdir large
fi

if [ ! -d "$SMALL" ]; then 
	mkdir small
fi

if [ ! -d "$SYN" ]; then
	mkdir syn
fi

cd $ZIPFILE
# nasasrb
wget https://nrvis.com/download/data/sc/sc-nasasrb.zip
unzip -n sc-nasasrb.zip
# fbwosn
wget https://nrvis.com/download/data/socfb/socfb-wosn-friends.zip
unzip -n socfb-wosn-friends.zip
# wiki-trust
wget https://nrvis.com/download/data/ia/ia-wiki-trust-dir.zip
unzip -n ia-wiki-trust-dir.zip
# shipsec5
wget https://nrvis.com/download/data/sc/sc-shipsec5.zip
unzip -n sc-shipsec5.zip
# pokec
wget https://nrvis.com/download/data/soc/soc-pokec.zip
unzip -n soc-pokec.zip
# baidu-baike
wget https://nrvis.com/download/data/web/web-baidu-baike.zip
unzip -n web-baidu-baike.zip

# web-sk
wget https://nrvis.com/download/data/web/web-sk-2005.zip
unzip -n web-sk-2005.zip
# citeseer
wget https://nrvis.com/download/data/ca/ca-citeseer.zip
unzip -n ca-citeseer.zip
# stanford
wget https://nrvis.com/download/data/web/web-Stanford.zip
unzip -n web-Stanford.zip
# dblp
wget https://nrvis.com/download/data/misc/com-dblp.zip
unzip -n com-dblp.zip
# digg
wget https://nrvis.com/download/data/soc/soc-digg.zip
unzip -n soc-digg.zip
# skitter
wget https://nrvis.com/download/data/tech/tech-as-skitter.zip
unzip -n tech-as-skitter.zip

rm readme.*

