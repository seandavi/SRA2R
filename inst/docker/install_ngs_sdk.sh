#/bin/bash
cd /tmp
mkdir ncbi
cd ncbi
git clone https://github.com/ncbi/ngs.git
git clone https://github.com/ncbi/ncbi-vdb.git
cd ngs
echo "********************************************"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "*          Making NGS-SDK                  *"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "********************************************"

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

echo "********************************************"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "*          Making NCBI-VCB                 *"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "*                                          *"
echo "********************************************"


# NCBI-VDB
cd ../ncbi-vdb/
./configure
make
make install
