#ifndef READ_H
#define READ_H
#include <Rcpp.h>
#include <ncbi-vdb/NGS.hpp>
#include <ngs-bam/ngs-bam.hpp>
#include <ngs/ErrorMsg.hpp>
#include <ngs/ReadCollection.hpp>
#include <ngs/ReadIterator.hpp>
#include <ngs/Read.hpp>

#include <math.h>
#include <iostream>

long getFastqCount(Rcpp::String acc, bool forward_to_r);
#endif
