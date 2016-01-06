// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// readCount
long readCount(Rcpp::String acc);
RcppExport SEXP SRA2R_readCount(SEXP accSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::String >::type acc(accSEXP);
    __result = Rcpp::wrap(readCount(acc));
    return __result;
END_RCPP
}
// reads
List reads(Rcpp::String acc, int n, SEXP lkup);
RcppExport SEXP SRA2R_reads(SEXP accSEXP, SEXP nSEXP, SEXP lkupSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::String >::type acc(accSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lkup(lkupSEXP);
    __result = Rcpp::wrap(reads(acc, n, lkup));
    return __result;
END_RCPP
}
// getPileUp
DataFrame getPileUp(Rcpp::String acc, Rcpp::String refname, int start, int stop, int MinPileUpDepth);
RcppExport SEXP SRA2R_getPileUp(SEXP accSEXP, SEXP refnameSEXP, SEXP startSEXP, SEXP stopSEXP, SEXP MinPileUpDepthSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::String >::type acc(accSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type refname(refnameSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type stop(stopSEXP);
    Rcpp::traits::input_parameter< int >::type MinPileUpDepth(MinPileUpDepthSEXP);
    __result = Rcpp::wrap(getPileUp(acc, refname, start, stop, MinPileUpDepth));
    return __result;
END_RCPP
}
// getFastqCount
long getFastqCount(Rcpp::String acc);
RcppExport SEXP SRA2R_getFastqCount(SEXP accSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::String >::type acc(accSEXP);
    __result = Rcpp::wrap(getFastqCount(acc));
    return __result;
END_RCPP
}
// getFastqReads
Rcpp::List getFastqReads(Rcpp::String acc);
RcppExport SEXP SRA2R_getFastqReads(SEXP accSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::String >::type acc(accSEXP);
    __result = Rcpp::wrap(getFastqReads(acc));
    return __result;
END_RCPP
}
// getFastqReadsWithQuality
Rcpp::List getFastqReadsWithQuality(Rcpp::String acc);
RcppExport SEXP SRA2R_getFastqReadsWithQuality(SEXP accSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::String >::type acc(accSEXP);
    __result = Rcpp::wrap(getFastqReadsWithQuality(acc));
    return __result;
END_RCPP
}
// getSRAReadsWithRegion
Rcpp::List getSRAReadsWithRegion(Rcpp::String acc, Rcpp::String refname, long start, long stop);
RcppExport SEXP SRA2R_getSRAReadsWithRegion(SEXP accSEXP, SEXP refnameSEXP, SEXP startSEXP, SEXP stopSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::String >::type acc(accSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type refname(refnameSEXP);
    Rcpp::traits::input_parameter< long >::type start(startSEXP);
    Rcpp::traits::input_parameter< long >::type stop(stopSEXP);
    __result = Rcpp::wrap(getSRAReadsWithRegion(acc, refname, start, stop));
    return __result;
END_RCPP
}
// getBamReadsWithRegion
Rcpp::List getBamReadsWithRegion(Rcpp::String acc, Rcpp::String refname, long start, long stop);
RcppExport SEXP SRA2R_getBamReadsWithRegion(SEXP accSEXP, SEXP refnameSEXP, SEXP startSEXP, SEXP stopSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::String >::type acc(accSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type refname(refnameSEXP);
    Rcpp::traits::input_parameter< long >::type start(startSEXP);
    Rcpp::traits::input_parameter< long >::type stop(stopSEXP);
    __result = Rcpp::wrap(getBamReadsWithRegion(acc, refname, start, stop));
    return __result;
END_RCPP
}
// getReference
DataFrame getReference(Rcpp::String acc);
RcppExport SEXP SRA2R_getReference(SEXP accSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::String >::type acc(accSEXP);
    __result = Rcpp::wrap(getReference(acc));
    return __result;
END_RCPP
}
