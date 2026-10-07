#ifndef EXTREMARANK_COUNT_READER_H
#define EXTREMARANK_COUNT_READER_H
#include <Rcpp.h>
#include <cmath>
#include <cstdint>
#include <limits>

namespace ercounts {
constexpr double max_exact = 9007199254740991.0;
struct Entry { uint32_t gene, unit, lo, hi; };
// A borrowed view: no conversion of an entire count matrix or stored-value vector.
class Reader {
    SEXP data;
    const int *p = nullptr, *rows = nullptr, *cols = nullptr;
    int kind = 0, outer = 0;
    R_xlen_t position = 0;
public:
    int nr, nc;
    R_xlen_t length;
    explicit Reader(SEXP x) {
        SEXP dim;
        if (Rf_isMatrix(x) && (TYPEOF(x)==REALSXP || TYPEOF(x)==INTSXP)) {
            dim=Rf_getAttrib(x,R_DimSymbol); data=x;
        } else if (Rf_isS4(x)) {
            Rcpp::S4 z(x); dim=z.slot("Dim"); data=z.slot("x");
            if (Rf_inherits(x,"dgCMatrix")) kind=1;
            else if (Rf_inherits(x,"dgRMatrix")) kind=2;
            else if (Rf_inherits(x,"dgTMatrix")) kind=3;
            else if (!Rf_inherits(x,"dgeMatrix")) Rcpp::stop("unsupported native count representation");
            if(kind==1 || kind==2) {
                SEXP pp=z.slot("p"), ii=z.slot(kind==1?"i":"j");
                if(TYPEOF(pp)!=INTSXP || TYPEOF(ii)!=INTSXP) Rcpp::stop("invalid sparse indices");
                p=INTEGER(pp); rows=INTEGER(ii);
                if(XLENGTH(ii)!=XLENGTH(data)) Rcpp::stop("invalid sparse index length");
            } else if(kind==3) {
                SEXP ii=z.slot("i"), jj=z.slot("j");
                if(TYPEOF(ii)!=INTSXP || TYPEOF(jj)!=INTSXP || XLENGTH(ii)!=XLENGTH(data) || XLENGTH(jj)!=XLENGTH(data)) Rcpp::stop("invalid triplet indices");
                rows=INTEGER(ii); cols=INTEGER(jj);
            }
        } else Rcpp::stop("unsupported native count representation");
        if(TYPEOF(dim)!=INTSXP || XLENGTH(dim)!=2) Rcpp::stop("invalid count dimensions");
        nr=INTEGER(dim)[0]; nc=INTEGER(dim)[1]; length=XLENGTH(data);
        if(nr<1 || nc<1 || (TYPEOF(data)!=REALSXP && TYPEOF(data)!=INTSXP)) Rcpp::stop("invalid numeric count matrix");
        if(kind==0 && length!=static_cast<R_xlen_t>(nr)*nc) Rcpp::stop("invalid dense count length");
        if(kind==1 || kind==2) {
            Rcpp::S4 z(x); SEXP pp=z.slot("p"); int n=kind==1?nc:nr;
            if(XLENGTH(pp)!=static_cast<R_xlen_t>(n)+1 || p[0]!=0 || p[n]!=length) Rcpp::stop("invalid compressed count pointers");
            for(int j=0;j<n;++j) if(p[j]<0 || p[j]>p[j+1]) Rcpp::stop("invalid compressed count pointers");
        }
    }
    bool next(int &gene,int &cell,double &value) {
        if(position==length) return false;
        if(kind==0) { gene=position%nr; cell=position/nr; }
        else if(kind==3) { gene=rows[position]; cell=cols[position]; }
        else {
            int n=kind==1?nc:nr;
            while(outer+1<n && position>=p[outer+1]) ++outer;
            if(kind==1) { gene=rows[position]; cell=outer; }
            else { gene=outer; cell=rows[position]; }
        }
        if(gene<0 || gene>=nr || cell<0 || cell>=nc) Rcpp::stop("count index out of bounds");
        value=TYPEOF(data)==REALSXP?REAL(data)[position]:INTEGER(data)[position]; ++position;
        if(!std::isfinite(value) || value<0 || value>max_exact || std::floor(value)!=value)
            Rcpp::stop("raw counts must be finite, nonnegative integers no larger than 2^53 - 1");
        return true;
    }
};
inline void check_units(const Reader &r,const Rcpp::IntegerVector &unit,int n) {
    if(n<1 || unit.size()!=r.nc) Rcpp::stop("invalid count grouping");
    for(int u:unit) if(u==NA_INTEGER || u<0 || u>n) Rcpp::stop("count group outside declared units");
}
}
#endif
