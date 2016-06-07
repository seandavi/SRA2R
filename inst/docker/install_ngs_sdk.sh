#/bin/bash
cd /tmp
mkdir ncbi
cd ncbi
git clone https://github.com/ncbi/ngs.git
git clone https://github.com/ncbi/ncbi-vdb.git
cd ngs
./configure
make -C ngs-sdk
make -C ngs-sdk install
make -C ngs-java
make -C ngs-java install
cd ngs-bam
./configure
cd ..
make -C ngs-bam
make -C ngs-bam install

# NCBI-VDB
cd ncbi-vdb/
./configure
make
make install
