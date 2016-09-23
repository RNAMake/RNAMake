to install RNAfold and RNAcofold from Vienna:

./configure
make
cd Progs/
gcc -DHAVE_CONFIG_H -I. -I.. -I./../H    -g -O2 -MT RNAcofold.o -MD -MP -MF .deps/RNAcofold.Tpo -c -o RNAcofold.o RNAcofold.c
gcc  -g -O2   -o RNAcofold RNAcofold.o ../lib/libRNA.a -lm 

