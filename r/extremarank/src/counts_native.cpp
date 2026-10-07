#include "count_reader.h"
using namespace Rcpp;
// [[Rcpp::export]]
bool cpp_validate_counts(SEXP x) {
    ercounts::Reader r(x); int g,c; double v;
    R_xlen_t n=0;
    while(r.next(g,c,v)) if((++n & 1048575)==0) checkUserInterrupt();
    return true;
}
// [[Rcpp::export]]
List cpp_pseudobulk_native(SEXP x,IntegerVector unit,int n_units) {
    ercounts::Reader r(x); ercounts::check_units(r,unit,n_units);
    NumericMatrix result(r.nr,n_units); int g,c; double v; R_xlen_t n=0;
    while(r.next(g,c,v)) {
        if((++n & 1048575)==0) checkUserInterrupt();
        int u=unit[c]; if(!u || !v) continue;
        double &a=result[static_cast<R_xlen_t>(r.nr)*(u-1)+g];
        if(a>ercounts::max_exact-v) stop("pseudobulk sum exceeds exact integer range 2^53 - 1");
        a+=v;
    }
    return List::create(_["counts"]=result,_["backend"]="native",_["gpu_buffer_bytes"]=0.0,
        _["gpu_command_seconds"]=0.0,_["gpu_commands"]=0,_["stored_values_scanned"]=static_cast<double>(n));
}

#include <array>
#include <unordered_map>
struct ErGroupHash {
    size_t operator()(const std::array<int,3>&k) const {
        size_t h=0; for(int v:k) h^=std::hash<int>{}(v)+0x9e3779b9+(h<<6)+(h>>2); return h;
    }
};
// [[Rcpp::export]]
List cpp_group_codes(List codes) {
    if(codes.size()!=3) stop("grouping needs three categorical code vectors");
    IntegerVector a=codes[0],b=codes[1],c=codes[2];
    if(a.size()!=b.size() || a.size()!=c.size()) stop("grouping code lengths differ");
    IntegerVector unit(a.size()); std::vector<int> first;
    std::unordered_map<std::array<int,3>,int,ErGroupHash> map;
    for(R_xlen_t i=0;i<a.size();++i) {
        if((i & 1048575)==0) checkUserInterrupt();
        if(a[i]<1 || b[i]<1 || c[i]<1) stop("grouping codes must be positive");
        std::array<int,3> key={a[i],b[i],c[i]}; auto found=map.find(key);
        if(found==map.end()) { int id=first.size()+1; map.emplace(key,id); first.push_back(i+1); unit[i]=id; }
        else unit[i]=found->second;
    }
    return List::create(_["unit"]=unit,_["first"]=wrap(first));
}
