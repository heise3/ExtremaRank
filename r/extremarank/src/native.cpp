#include <Rcpp.h>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <functional>
#include <numeric>
#include <map>
#include <tuple>
#include <set>
#include <string>
#include <vector>

using namespace Rcpp;
using Z = boost::multiprecision::cpp_int;
using Rat = boost::multiprecision::cpp_rational;
using Quad = boost::multiprecision::cpp_bin_float_quad;
using VI = std::vector<int>;
int sgn(const Z& x) { return x > 0 ? 1 : x < 0 ? -1 : 0; }
int sgn(const Rat& x) { return x > 0 ? 1 : x < 0 ? -1 : 0; }
Rat ratio(const Z& n, const Z& d) { return Rat(n, d); }

// Decode IEEE-754 bits, retaining subnormals and all cancellation exactly.
struct Feature {
  std::vector<Z> x; VI sorted, squareSorted;
  Z den = 1, s[2], q[2];
  Feature(NumericMatrix a, int j, const VI& labels) {
    int n=a.ncol(), emin=0;
    std::vector<uint64_t> mant(n); VI exponent(n), negative(n);
    for(int i=0;i<n;++i) {
      double value=a(j,i); if(!std::isfinite(value)) stop("all prepared inputs must be finite");
      uint64_t bits; std::memcpy(&bits,&value,sizeof(bits));
      int e=int((bits>>52)&2047); negative[i]=int(bits>>63);
      mant[i]=bits&((uint64_t(1)<<52)-1);
      if(e) { mant[i] |= uint64_t(1)<<52; exponent[i]=e-1023-52; }
      else exponent[i]=-1074;
      if(mant[i]) emin=std::min(emin,exponent[i]);
    }
    den <<= -emin; x.resize(n);
    for(int i=0;i<n;++i) {
      if(mant[i]) { x[i]=mant[i]; x[i]<<=exponent[i]-emin; if(negative[i]) x[i]=-x[i]; }
      s[labels[i]]+=x[i]; q[labels[i]]+=x[i]*x[i];
    }
    sorted.resize(n); std::iota(sorted.begin(),sorted.end(),0);
    std::sort(sorted.begin(),sorted.end(),[&](int a,int b){return x[a]!=x[b]?x[a]<x[b]:a<b;});
    squareSorted=sorted; std::sort(squareSorted.begin(),squareSorted.end(),[&](int a,int b){Z aa=x[a]*x[a],bb=x[b]*x[b];return aa!=bb?aa<bb:a<b;});
  }
};
struct Score { int sign=0; Z num=0, den=1; Rat effect=0; };
Score paired(const Z& s,const Z& q,int n,const Z& scale) {
  Score a; a.sign=sgn(s); a.num=(n-1)*s*s; a.den=n*q-s*s; a.effect=ratio(s,Z(n)*scale); return a;
}
Score statistic(const Rat& d,const Rat& v) {
  Score a; a.sign=sgn(d); a.effect=d;
  if(numerator(v)==0) { a.num=a.sign?1:0; a.den=0; return a; }
  Rat r=d*d; r/=v; a.num=numerator(r); a.den=denominator(r); return a;
}
int compare(const Score& a,const Score& b,bool absolute=false) {
  int sa=absolute?std::abs(a.sign):a.sign, sb=absolute?std::abs(b.sign):b.sign;
  if(sa!=sb) return sa>sb?1:-1; if(!sa) return 0;
  int c;
  if(a.den==0 || b.den==0) c=(a.den==0)-(b.den==0);
  else { Z lhs=a.num*b.den, rhs=b.num*a.den; c=lhs>rhs?1:lhs<rhs?-1:0; }
  return sa*c;
}
double display(const Rat& x) {
  Quad v(numerator(x).str()); v/=Quad(denominator(x).str()); return v.convert_to<double>();
}
double displayT(const Score& a) {
  if(!a.sign) return 0; if(a.den==0) return a.sign*R_PosInf;
  Quad v(a.num.str()); v/=Quad(a.den.str()); v=sqrt(v); return a.sign*v.convert_to<double>();
}
struct Bound { Score lo,hi; };
struct Ext { Bound bound; VI lo,hi; int candidates=0; };
Ext pairedBound(const Feature& f,const VI& inc,const VI& opt,int m,bool witness=false) {
  int u=opt.size(), n=inc.size()+m; if(m<0 || m>u || n<2) stop("invalid retained cardinality");
  std::vector<bool> take(f.x.size(),false);for(int i:opt)take[i]=true; VI order;for(int i:f.sorted)if(take[i])order.push_back(i);
  std::vector<Z> ps(u+1),pq(u+1); Z sf=0,qf=0;
  for(int i:inc) { sf+=f.x[i]; qf+=f.x[i]*f.x[i]; }
  for(int i=0;i<u;++i) { ps[i+1]=ps[i]+f.x[order[i]]; pq[i+1]=pq[i]+f.x[order[i]]*f.x[order[i]]; }
  Ext out; bool first=true; int lokind=0,loarg=0,hikind=0,hiarg=0;
  auto inspect=[&](const Z& s,const Z& q,int kind,int arg) {
    Score a=paired(sf+s,qf+q,n,f.den);
    if(first || compare(a,out.bound.lo)<0) {out.bound.lo=a;lokind=kind;loarg=arg;}
    if(first || compare(a,out.bound.hi)>0) {out.bound.hi=a;hikind=kind;hiarg=arg;}
    first=false; ++out.candidates;
  };
  if(m==0 || m==u) inspect(ps[m],pq[m],0,0);
  else {
    for(int start=0;start<=u-m;++start) inspect(Z(ps[start+m]-ps[start]),Z(pq[start+m]-pq[start]),0,start);
    for(int j=1;j<m;++j) inspect(Z(ps[j]+ps[u]-ps[u-(m-j)]),Z(pq[j]+pq[u]-pq[u-(m-j)]),1,j);
  }
  if(witness) {
    auto unpack=[&](int kind,int arg) {VI ids=inc;
      if(!kind) ids.insert(ids.end(),order.begin()+arg,order.begin()+arg+m);
      else {ids.insert(ids.end(),order.begin(),order.begin()+arg);ids.insert(ids.end(),order.end()-(m-arg),order.end());}
      std::sort(ids.begin(),ids.end()); return ids;
    }; out.lo=unpack(lokind,loarg);out.hi=unpack(hikind,hiarg);
  }
  return out;
}
struct Box { Rat slo,shi,vlo,vhi; };
Box groupBox(const Feature& f,const VI& inc,const VI& opt,int m) {
  int n=inc.size()+m; if(n<2 || m<0 || m>int(opt.size())) stop("invalid Welch node");
  Z sf=0,qf=0; std::vector<Z> xs,squares,allowed;
  std::vector<bool> fixed(f.x.size(),false),optional(f.x.size(),false);
  for(int i:inc) {sf+=f.x[i];qf+=f.x[i]*f.x[i];fixed[i]=true;}
  for(int i:opt)optional[i]=true;
  for(int i:f.sorted){if(optional[i])xs.push_back(f.x[i]);if(optional[i]||fixed[i])allowed.push_back(f.x[i]);}
  for(int i:f.squareSorted)if(optional[i])squares.push_back(f.x[i]*f.x[i]);
  Z slo=sf,shi=sf,qlo=qf,qhi=qf;
  for(int i=0;i<m;++i) {slo+=xs[i];shi+=xs[xs.size()-1-i];qlo+=squares[i];qhi+=squares[squares.size()-1-i];}
  Z s=0,q=0;for(int i=0;i<n;++i){s+=allowed[i];q+=allowed[i]*allowed[i];}
  Z floor=Z(n)*q-s*s;
  for(int start=1;start<=int(allowed.size())-n;++start) {
    Z a=allowed[start-1],b=allowed[start+n-1];s=s-a+b;q=q-a*a+b*b;Z v=Z(n)*q-s*s;if(v<floor)floor=v;
  }
  // Pairwise decomposition: fixed-fixed + fixed-optional + optional-optional.
  // Each separate minimum is a lower bound even when not jointly attainable.
  int fixedCount=inc.size();Z constrained=Z(fixedCount)*qf-sf*sf;
  if(fixedCount>0){std::vector<Z> cross;for(const Z& x:xs)cross.push_back(Z(fixedCount)*x*x-2*sf*x+qf);
  std::sort(cross.begin(),cross.end());for(int i=0;i<m;++i)constrained+=cross[i];
  if(m>1){Z os=0,oq=0;for(int i=0;i<m;++i){os+=xs[i];oq+=xs[i]*xs[i];}Z best=Z(m)*oq-os*os;
    for(int i=1;i<=int(xs.size())-m;++i){Z a=xs[i-1],b=xs[i+m-1];os=os-a+b;oq=oq-a*a+b*b;Z v=Z(m)*oq-os*os;if(v<best)best=v;}constrained+=best;}
  floor=std::max(floor,constrained);}
  Z a=slo*slo,b=shi*shi;Z minSquare=slo<=0 && shi>=0?Z(0):std::min(a,b);
  Z lo=Z(n)*qlo-std::max(a,b);lo=std::max(Z(0),std::max(floor,lo));
  Z hi=Z(n)*qhi-minSquare;hi=std::max(Z(0),hi);Z vd=Z(n)*n*(n-1);
  return {ratio(slo,Z(n)),ratio(shi,Z(n)),ratio(lo,vd),ratio(hi,vd)};
}
Z chooseZ(int n,int r) {if(r<0 || r>n)return 0;r=std::min(r,n-r);Z x=1;for(int i=1;i<=r;++i){x*=n-r+i;x/=i;}return x;}
// Streaming combinations: false from visitor stops before materializing a universe.
bool combinations(const VI& pool,int r,const std::function<bool(const VI&)>& visit) {
  if(r<0 || r>int(pool.size())) return true;
  VI positions(r);std::iota(positions.begin(),positions.end(),0);
  if(!r) return visit(VI());
  while(true) { VI value;for(int i:positions)value.push_back(pool[i]);if(!visit(value))return false;
    int j=r-1;while(j>=0 && positions[j]==int(pool.size())-r+j)--j;if(j<0)break;
    ++positions[j];for(int i=j+1;i<r;++i)positions[i]=positions[i-1]+1;
  } return true;
}
struct Query {int feature,property;};
class Problem {
public:
  int n,g,k,ng;bool isPaired;std::string direction;VI labels,counts,pools[2],baseOrder,baseSign;std::vector<bool> baseMember;
  std::vector<std::string> ids;std::vector<Feature> data;std::vector<Score> baseline;
  Problem(NumericMatrix x,CharacterVector names,IntegerVector groups,bool pairedDesign,std::string dir,int kk):
    n(x.ncol()),g(x.nrow()),k(kk),ng(pairedDesign?1:2),isPaired(pairedDesign),direction(dir) {
    if(g<1 || n<2 || names.size()!=g || groups.size()!=n || k<1 || k>g) stop("invalid native matrix/IDs/top-K");
    if(dir!="up" && dir!="down" && dir!="absolute")stop("invalid ranking direction");
    counts.assign(ng,0);std::set<std::string> seen;
    for(int j=0;j<g;++j){if(CharacterVector::is_na(names[j]))stop("missing feature ID");std::string id=as<std::string>(names[j]);if(id.empty() || !seen.insert(id).second)stop("duplicate/empty feature IDs");ids.push_back(id);}
    for(int i=0;i<n;++i){int label=groups[i];if(label<0 || label>=ng)stop("invalid group label");labels.push_back(label);++counts[label];pools[label].push_back(i);}
    for(int count:counts)if(count<2)stop("each native group needs at least two units");
    data.reserve(g);for(int j=0;j<g;++j){if(j%1024==0)checkUserInterrupt();data.emplace_back(x,j,labels);}
    baseline=evaluate(VI());baseOrder=rank(baseline,g);baseSign.resize(g);baseMember.assign(g,false);
    for(int j=0;j<g;++j)baseSign[j]=baseline[j].sign;for(int i=0;i<k;++i)baseMember[baseOrder[i]]=true;
  }
  std::vector<Score> evaluate(const VI& removed) const {
    VI dcounts(ng,0);std::set<int> seen;for(int i:removed){if(i<0 || i>=n || !seen.insert(i).second)stop("invalid deletion positions");++dcounts[labels[i]];}
    for(int h=0;h<ng;++h)if(counts[h]-dcounts[h]<2)stop("deletion leaves fewer than two observations");
    std::vector<Score> out;out.reserve(g);
    for(int j=0;j<g;++j){if(j%1024==0)checkUserInterrupt();const Feature& f=data[j];Z s[2]={f.s[0],f.s[1]},q[2]={f.q[0],f.q[1]};
      for(int i:removed){s[labels[i]]-=f.x[i];q[labels[i]]-=f.x[i]*f.x[i];}
      int nt=counts[0]-dcounts[0];if(isPaired){out.push_back(paired(s[0],q[0],nt,f.den));continue;}
      int nr=counts[1]-dcounts[1];Z d=s[0]*nr-s[1]*nt;Score a;a.sign=sgn(d);a.num=d*d*(nt-1)*(nr-1);
      Z vt=nt*q[0]-s[0]*s[0],vr=nr*q[1]-s[1]*s[1];a.den=vt*nr*nr*(nr-1)+vr*nt*nt*(nt-1);a.effect=ratio(d,Z(nt)*nr*f.den);out.push_back(a);
    }return out;
  }
  bool better(int a,const Score& sa,int b,const Score& sb)const {
    int c=compare(sa,sb,direction=="absolute");if(direction=="down")c=-c;return c>0 || (!c && ids[a]<ids[b]);
  }
  VI rank(const std::vector<Score>& scores,int limit)const {
    VI order(g);std::iota(order.begin(),order.end(),0);auto cmp=[&](int a,int b){return better(a,scores[a],b,scores[b]);};
    if(limit==g)std::sort(order.begin(),order.end(),cmp);else{std::partial_sort(order.begin(),order.begin()+limit,order.end(),cmp);order.resize(limit);}return order;
  }
  std::vector<VI> allocations(int r)const {
    if(isPaired)return {{n-r}};std::vector<VI> out;
    for(int a=std::max(0,r-counts[1]+2);a<=std::min(r,counts[0]-2);++a)out.push_back({counts[0]-a,counts[1]-(r-a)});return out;
  }
  Z feasible(int r)const {Z count=0;for(const VI& kept:allocations(r)){Z a=1;for(int h=0;h<ng;++h)a*=chooseZ(counts[h],kept[h]);count+=a;}return count;}
  bool deletions(int r,const std::function<bool(const VI&)>& visit)const {
    if(isPaired)return combinations(pools[0],r,visit);
    for(const VI& kept:allocations(r)){
      bool result=combinations(pools[0],counts[0]-kept[0],[&](const VI& dt){return combinations(pools[1],counts[1]-kept[1],[&](const VI& dr){VI d=dt;d.insert(d.end(),dr.begin(),dr.end());std::sort(d.begin(),d.end());return visit(d);});});
      if(!result)return false;
    }return true;
  }
  std::vector<Bound> bounds(const VI& inc,const VI& opt,const VI& need)const {
    std::vector<Bound> out;out.reserve(g);VI fi[2],fo[2];for(int i:inc)fi[labels[i]].push_back(i);for(int i:opt)fo[labels[i]].push_back(i);
    for(int j=0;j<g;++j){if(j%1024==0)checkUserInterrupt();const Feature& f=data[j];
      if(isPaired){out.push_back(pairedBound(f,inc,opt,need[0]).bound);continue;}
      Box a=groupBox(f,fi[0],fo[0],need[0]),b=groupBox(f,fi[1],fo[1],need[1]);Rat lo=a.slo-b.shi,hi=a.shi-b.slo,vl=a.vlo+b.vlo,vh=a.vhi+b.vhi;
      try {out.push_back({statistic(lo,lo<0?vl:vh),statistic(hi,hi>0?vl:vh)});} catch(const std::exception& e){stop(std::string("Welch score bound: ")+e.what());}
    }return out;
  }
  Bound target(Bound a)const {
    if(direction=="down")std::swap(a.lo,a.hi);
    else if(direction=="absolute"){
      int c=compare(a.lo,a.hi,true);Score low=c<=0?a.lo:a.hi,high=c<=0?a.hi:a.lo;
      if(a.lo.sign<=0 && a.hi.sign>=0)low=Score();a={low,high};
    }return a;
  }
  bool proves(Query q,const std::vector<Bound>& bound)const {
    int j=q.feature;if(q.property==1){int sign=baseSign[j];return sign>0?bound[j].lo.sign>0:sign<0?bound[j].hi.sign<0:bound[j].lo.sign==0 && bound[j].hi.sign==0;}
    Bound a=target(bound[j]);int possible=0,definite=0;
    for(int l=0;l<g;++l)if(l!=j){Bound b=target(bound[l]);possible+=better(l,b.hi,j,a.lo);definite+=better(l,b.lo,j,a.hi);}
    return baseMember[j]?possible<k:definite>=k;
  }
};
struct Prover {
  const Problem& p;const std::vector<Bound>& b;bool indexed;std::vector<Score> lo,hi;VI lows,highs;
  Prover(const Problem& pp,const std::vector<Bound>& bb,const VI& active,const std::vector<Query>& qs):p(pp),b(bb),indexed(false){
    int members=0;for(int q:active)members+=qs[q].property==0;if(members<32)return;indexed=true;
    for(int j=0;j<p.g;++j){Bound z=p.target(b[j]);lo.push_back(z.lo);hi.push_back(z.hi);}
    lows=p.rank(lo,p.g);highs=p.rank(hi,p.g);
  }
  int count(const VI& order,const std::vector<Score>& values,int j,const Score& score)const {
    int left=0,right=order.size();while(left<right){int mid=(left+right)/2;if(p.better(order[mid],values[order[mid]],j,score))left=mid+1;else right=mid;}return left;
  }
  bool operator()(Query q)const {
    if(!indexed || q.property==1)return p.proves(q,b);int j=q.feature;
    int possible=count(highs,hi,j,lo[j])-int(p.better(j,hi[j],j,lo[j]));int definite=count(lows,lo,j,hi[j]);return p.baseMember[j]?possible<p.k:definite>=p.k;
  }
};
struct Witness { VI deleted,order;int sign=0;bool found=false;};
struct State {int lower=1,upper=0;Witness witness;};
struct Node {VI inc,opt,need,queries;};
SEXP nullable(int value){return value?wrap(value):R_NilValue;}
VI onebased(VI values){for(int& i:values)++i;return values;}
CharacterVector namedIds(const VI& positions,const std::vector<std::string>& ids){CharacterVector out(positions.size());for(size_t i=0;i<positions.size();++i)out[i]=ids[positions[i]];return out;}
List packWitness(const Witness& w,const Problem& p,CharacterVector unitIds,bool sign=true){
  CharacterVector units(w.deleted.size());for(size_t i=0;i<w.deleted.size();++i)units[i]=unitIds[w.deleted[i]];
  List out=List::create(_["deleted_indices"]=onebased(w.deleted),_["deleted_units"]=units,_["topk"]=namedIds(w.order,p.ids));if(sign)out["effect_sign"]=w.sign;return out;
}
List packProperty(const State& a,int budget,const Problem& p,CharacterVector unitIds){
  std::string status=a.upper?"REFUTED":a.lower>budget?"CERTIFIED":"UNRESOLVED";
  return List::create(_["status"]=status,_["minimum_change_lower_bound"]=a.lower,_["minimum_change_upper_bound"]=nullable(a.upper),_["exact_minimum_change"]=nullable(a.upper && a.lower==a.upper?a.upper:0),_["certified_through"]=std::min(budget,a.lower-1),_["witness"]=a.witness.found?SEXP(packWitness(a.witness,p,unitIds)):R_NilValue);
}

// [[Rcpp::export]]
List cpp_audit(NumericMatrix x,CharacterVector ids,CharacterVector unit_ids,IntegerVector labels,bool paired_design,int k,int budget,std::string direction,IntegerVector selected,int max_nodes,int max_scenarios,int witness_trials=8,bool use_bounds=true,bool influence_order=false){
  Problem p(x,ids,labels,paired_design,direction,k);if(unit_ids.size()!=p.n || budget<0 || budget>p.n-2*p.ng || max_nodes<1 || max_scenarios<1 || witness_trials<0)stop("invalid budget, unit IDs or limits");
  std::vector<Query> queries;VI chosen;std::set<int> selectedSet;
  VI requested;for(int value:selected)requested.push_back(value);if(requested.empty())for(int i=0;i<k;++i)requested.push_back(p.baseOrder[i]+1);
  for(int value:requested){int j=value-1;if(j<0 || j>=p.g || !selectedSet.insert(j).second)stop("invalid queried features");chosen.push_back(j);queries.push_back({j,0});queries.push_back({j,1});}if(chosen.empty())stop("no queried features");
  std::vector<State> state(queries.size());std::set<VI> evaluated;Witness first;int nodes=0,checked=0,pruned=0;bool stopped=false;
  auto inspect=[&](const VI& raw,int r){VI d=raw;std::sort(d.begin(),d.end());if(!evaluated.insert(d).second)return;++checked;auto scores=p.evaluate(d);VI order=p.rank(scores,k);std::vector<bool> members(p.g,false);for(int j:order)members[j]=true;
    if(members!=p.baseMember && !first.found)first={d,order,0,true};
    for(size_t q=0;q<queries.size();++q){int j=queries[q].feature;bool changed=queries[q].property==0?members[j]!=p.baseMember[j]:scores[j].sign!=p.baseSign[j];
      if(changed && (!state[q].upper || r<state[q].upper)){state[q].upper=r;state[q].witness={d,order,scores[j].sign,true};}}
  };
  VI all(p.n);std::iota(all.begin(),all.end(),0);
  if(influence_order){std::vector<double> weight(p.n,0);for(int j:chosen){Z largest=0;for(const Z& x:p.data[j].x){Z a=abs(x);if(a>largest)largest=a;}if(largest!=0)for(int i=0;i<p.n;++i)weight[i]+=display(ratio(abs(p.data[j].x[i]),largest));}
    std::stable_sort(all.begin(),all.end(),[&](int a,int b){return weight[a]>weight[b];});for(int h=0;h<p.ng;++h){p.pools[h].clear();for(int i:all)if(p.labels[i]==h)p.pools[h].push_back(i);}}
  for(int r=1;r<=budget;++r){
    VI active;for(size_t q=0;q<queries.size();++q)if(!state[q].upper)active.push_back(q);if(active.empty())break;
    int trials=0;p.deletions(r,[&](const VI& d){if(trials>=witness_trials || checked>=max_scenarios)return false;inspect(d,r);++trials;return true;});
    active.erase(std::remove_if(active.begin(),active.end(),[&](int q){return state[q].upper!=0;}),active.end());
    std::set<int> pending;auto allocations=p.allocations(r);
    if(p.feasible(r)<=2000){
      std::set<int> needs;if(!use_bounds)needs.insert(active.begin(),active.end());
      if(use_bounds)for(const VI& keep:allocations){
        if(nodes>=max_nodes){stopped=true;needs.insert(active.begin(),active.end());break;}
        ++nodes;auto b=p.bounds(VI(),all,keep);Prover proof(p,b,active,queries);
        for(int q:active)if(proof(queries[q]))++pruned;else needs.insert(q);
      }
      if(stopped){for(int q:needs)if(!state[q].upper)pending.insert(q);}
      else if(!needs.empty())p.deletions(r,[&](const VI& d){
        bool unresolved=false;for(int q:needs)unresolved|=!state[q].upper;if(!unresolved)return false;if(evaluated.count(d))return true;
        if(checked>=max_scenarios || nodes>=max_nodes){stopped=true;for(int q:needs)if(!state[q].upper)pending.insert(q);return false;}
        ++nodes;inspect(d,r);return true;
      });
    }else{
      std::vector<Node> stack;for(const VI& keep:allocations)stack.push_back({VI(),all,keep,active});
      while(!stack.empty()){
        if(nodes>=max_nodes || checked>=max_scenarios){for(const Node& node:stack)for(int q:node.queries)if(!state[q].upper)pending.insert(q);stopped=true;break;}
        Node node=std::move(stack.back());stack.pop_back();node.queries.erase(std::remove_if(node.queries.begin(),node.queries.end(),[&](int q){return state[q].upper!=0;}),node.queries.end());if(node.queries.empty())continue;
        bool feasible=true;
        for(int h=0;h<p.ng;++h){VI pool;for(int i:node.opt)if(p.labels[i]==h)pool.push_back(i);int need=node.need[h];if(need<0 || need>int(pool.size())){feasible=false;break;}
          if(need==0 || need==int(pool.size())){if(need)node.inc.insert(node.inc.end(),pool.begin(),pool.end());node.opt.erase(std::remove_if(node.opt.begin(),node.opt.end(),[&](int i){return p.labels[i]==h;}),node.opt.end());node.need[h]=0;}}
        if(!feasible)continue;std::sort(node.inc.begin(),node.inc.end());++nodes;
        if(node.opt.empty()){VI d;for(int i:all)if(!std::binary_search(node.inc.begin(),node.inc.end(),i))d.push_back(i);std::sort(d.begin(),d.end());inspect(d,r);continue;}
        if(use_bounds){auto b=p.bounds(node.inc,node.opt,node.need);Prover proof(p,b,node.queries,queries);VI remain;for(int q:node.queries)if(proof(queries[q]))++pruned;else remain.push_back(q);node.queries=remain;}
        if(node.queries.empty())continue;int pivot=node.opt[0],group=p.labels[pivot];node.opt.erase(node.opt.begin());stack.push_back(node);
        if(node.need[group]>0){node.inc.push_back(pivot);--node.need[group];stack.push_back(std::move(node));}
      }
    }
    for(int q:active)if(!state[q].upper && !pending.count(q))state[q].lower=r+1;if(stopped)break;
  }
  VI ranks(p.g);for(int i=0;i<p.g;++i)ranks[p.baseOrder[i]]=i+1;List items(chosen.size());bool allOriginal=true,allCertified=true;int lower=budget+1;
  for(int i=0;i<k;++i)if(!selectedSet.count(p.baseOrder[i]))allOriginal=false;
  for(size_t q=0;q<chosen.size();++q){int j=chosen[q];items[q]=List::create(_["feature_id"]=p.ids[j],_["baseline_rank"]=ranks[j],_["baseline_topk"]=bool(p.baseMember[j]),_["baseline_effect_sign"]=p.baseSign[j],_["baseline_t"]=displayT(p.baseline[j]),_["baseline_effect"]=display(p.baseline[j].effect),_["membership"]=packProperty(state[2*q],budget,p,unit_ids),_["direction"]=packProperty(state[2*q+1],budget,p,unit_ids));
    if(p.baseMember[j]){lower=std::min(lower,state[2*q].lower);allCertified&=!state[2*q].upper && state[2*q].lower>budget;}}
  if(!allOriginal)lower=1;int upper=first.found?first.deleted.size():0;Z total=0;for(int r=1;r<=budget;++r)total+=p.feasible(r);
  return List::create(_["status"]=first.found?"REFUTED":allOriginal && allCertified?"CERTIFIED":"UNRESOLVED",_["baseline_topk"]=namedIds(VI(p.baseOrder.begin(),p.baseOrder.begin()+k),p.ids),_["features"]=items,_["budget"]=budget,_["top_k"]=k,_["design"]=paired_design?"paired":"welch",_["direction"]=direction,_["unit_ids"]=unit_ids,_["first_topk_witness"]=first.found?SEXP(packWitness(first,p,unit_ids,false)):R_NilValue,_["topk_minimum_change"]=List::create(_["lower_bound"]=lower,_["upper_bound"]=nullable(upper),_["exact"]=nullable(lower==upper?upper:0)),_["nodes"]=nodes,_["scenarios_checked"]=checked,_["query_bounds_pruned"]=pruned,_["feasible_deletion_sets"]=total.str(),_["limits"]=List::create(_["max_nodes"]=max_nodes,_["max_scenarios"]=max_scenarios),_["arithmetic"]="native C++ arbitrary-precision dyadic integers and rational bounds; exact finite binary64 ranking decisions",_["index_base"]=1);
}

// [[Rcpp::export]]
List cpp_extrema(NumericVector values,NumericVector included,int choose){
  int f=included.size(),u=values.size();NumericMatrix x(1,f+u);for(int i=0;i<f;++i)x(0,i)=included[i];for(int i=0;i<u;++i)x(0,f+i)=values[i];
  VI labels(f+u,0),inc(f),opt(u);std::iota(inc.begin(),inc.end(),0);std::iota(opt.begin(),opt.end(),f);Feature feature(x,0,labels);Ext e=pairedBound(feature,inc,opt,choose,true);
  auto pack=[&](const Score& s,const VI& kept){Z sum=0,squares=0;VI optional;for(int i:kept){sum+=feature.x[i];squares+=feature.x[i]*feature.x[i];if(i>=f)optional.push_back(i-f+1);}
    return List::create(_["t"]=displayT(s),_["effect"]=display(s.effect),_["exact_sum"]=ratio(sum,feature.den).str(),_["exact_sum_squares"]=ratio(squares,feature.den*feature.den).str(),_["squared_t_numerator"]=s.num.str(),_["squared_t_denominator"]=s.den.str(),_["optional_indices"]=optional,_["included_count"]=f,_["retained_count"]=f+choose);};
  return List::create(_["minimum"]=pack(e.bound.lo,e.lo),_["maximum"]=pack(e.bound.hi,e.hi),_["candidate_count"]=e.candidates,_["total_subsets"]=chooseZ(u,choose).str(),_["method"]="conditional signed two-family exact extrema");
}

// [[Rcpp::export]]
List cpp_evaluate(NumericMatrix x,CharacterVector ids,IntegerVector labels,bool paired_design,std::string direction,IntegerVector deleted){
  Problem p(x,ids,labels,paired_design,direction,1);VI d;for(int i:deleted)d.push_back(i-1);auto scores=p.evaluate(d);VI order=p.rank(scores,p.g);List records(p.g);
  for(int j=0;j<p.g;++j){Z num=scores[j].num,den=scores[j].den;if(den!=0){Rat canonical=ratio(num,den);num=numerator(canonical);den=denominator(canonical);}
    records[j]=List::create(_["sign"]=scores[j].sign,_["numerator"]=num.str(),_["denominator"]=den.str(),_["effect"]=scores[j].effect.str(),_["t"]=displayT(scores[j]));}
  return List::create(_["order"]=onebased(order),_["scores"]=records);
}

// Internal validation surface: exact signed scalar enclosures at conditional
// search nodes. Public users work through the validated R audit interface.
// [[Rcpp::export]]
List cpp_bounds(NumericMatrix x,CharacterVector ids,IntegerVector labels,bool paired_design,
                IntegerVector included,IntegerVector optional,IntegerVector needed){
  Problem p(x,ids,labels,paired_design,"up",1);VI inc,opt,need;std::set<int> seen;
  for(int i:included){if(i<1 || i>p.n || !seen.insert(i).second)stop("invalid included units");inc.push_back(i-1);}
  for(int i:optional){if(i<1 || i>p.n || !seen.insert(i).second)stop("invalid optional units");opt.push_back(i-1);}
  for(int i:needed)need.push_back(i);if(int(need.size())!=p.ng)stop("needed counts disagree with design");
  std::vector<Bound> b;try{b=p.bounds(inc,opt,need);}catch(const std::exception& e){stop(std::string("conditional bounds: ")+e.what());}List out(p.g);
  auto pack=[](const Score& a){Z n=a.num,d=a.den;if(d!=0){Rat r=ratio(n,d);n=numerator(r);d=denominator(r);}
    return List::create(_["sign"]=a.sign,_["numerator"]=n.str(),_["denominator"]=d.str());};
  for(int j=0;j<p.g;++j)out[j]=List::create(_["lower"]=pack(b[j].lo),_["upper"]=pack(b[j].hi));return out;
}

// [[Rcpp::export]]
List cpp_diagnostics(NumericMatrix x,CharacterVector ids,CharacterVector unit_ids,IntegerVector labels,bool paired_design,int k,std::string direction,int max_diagnostics=1000000){
  Problem p(x,ids,labels,paired_design,direction,k);if(unit_ids.size()!=p.n || max_diagnostics<0)stop("invalid diagnostic limit or unit IDs");
  VI ranks(p.g),best(p.g),worst(p.g),selected(p.g,0),changedSign(p.g,0),skipped,unrun;List scenarios;int nScenarios=0,nChanged=0;
  for(int i=0;i<p.g;++i)ranks[p.baseOrder[i]]=i+1;best=worst=ranks;
  for(int i=0;i<p.n;++i){if(p.counts[p.labels[i]]<=2){skipped.push_back(i);continue;}if(nScenarios>=max_diagnostics){unrun.push_back(i);continue;}
    auto scores=p.evaluate({i});VI order=p.rank(scores,p.g);std::vector<bool> after(p.g,false);for(int r=0;r<p.g;++r){int j=order[r];best[j]=std::min(best[j],r+1);worst[j]=std::max(worst[j],r+1);if(r<k)after[j]=true;changedSign[j]+=scores[j].sign!=p.baseSign[j];}
    VI exits,entries;for(int r=0;r<k;++r){int j=p.baseOrder[r];if(!after[j])exits.push_back(j);j=order[r];if(!p.baseMember[j])entries.push_back(j);}for(int j=0;j<p.g;++j)selected[j]+=after[j];
    bool changed=after!=p.baseMember;++nScenarios;nChanged+=changed;double overlap=k-exits.size();
    scenarios.push_back(List::create(_["deleted_unit"]=as<std::string>(unit_ids[i]),_["deleted_index"]=i+1,_["topk_changed"]=changed,_["jaccard"]=overlap/(2*k-overlap),_["exited_features"]=namedIds(exits,p.ids),_["entered_features"]=namedIds(entries,p.ids)));
  }
  NumericVector t(p.g),fraction(p.g);LogicalVector original(p.g);for(int j=0;j<p.g;++j){t[j]=displayT(p.baseline[j]);fraction[j]=nScenarios?double(selected[j])/nScenarios:NA_REAL;original[j]=p.baseMember[j];}
  DataFrame features=DataFrame::create(_["feature_id"]=ids,_["baseline_rank"]=ranks,_["baseline_t"]=t,_["best_observed_rank"]=best,_["worst_observed_rank"]=worst,_["loo_selected_count"]=selected,_["loo_scenarios"]=IntegerVector(p.g,nScenarios),_["loo_selection_fraction"]=fraction,_["loo_sign_changes"]=changedSign,_["baseline_topk"]=original);
  return List::create(_["features"]=features,_["scenarios"]=scenarios,_["n_scenarios"]=nScenarios,_["n_changed"]=nChanged,_["skipped_units"]=namedIds(skipped,as<std::vector<std::string>>(unit_ids)),_["not_run_units"]=namedIds(unrun,as<std::vector<std::string>>(unit_ids)),_["enumeration_capped"]=!unrun.empty(),_["feasible_scenarios"]=p.n-int(skipped.size()),_["scope"]="baseline plus evaluated feasible one-unit deletions; observed fractions are not probabilities or multi-deletion bounds");
}

// [[Rcpp::export]]
List cpp_plan(IntegerVector labels,bool paired_design,int budget){
  int ng=paired_design?1:2;VI counts(ng,0);for(int i:labels){if(i<0||i>=ng)stop("invalid labels");++counts[i];}
  if(budget<0 || budget>labels.size()-2*ng)stop("invalid budget");
  CharacterVector by(budget);Z total=0;
  for(int r=1;r<=budget;++r){Z v=0;if(paired_design)v=chooseZ(labels.size(),r);else for(int a=std::max(0,r-counts[1]+2);a<=std::min(r,counts[0]-2);++a)v+=chooseZ(counts[0],a)*chooseZ(counts[1],r-a);by[r-1]=v.str();total+=v;}
  return List::create(_["feasible_deletions"]=total.str(),_["by_cardinality"]=by,_["units"]=labels.size());
}

// [[Rcpp::export]]
List cpp_block_plan(IntegerVector a,IntegerVector b,int budget,bool paired_design){
 if(a.size()!=b.size()||budget<0||budget>a.size())stop("invalid block plan");
 int na=0,nb=0;for(int i=0;i<a.size();++i){if(a[i]<0||b[i]<0)stop("negative block counts");na+=a[i];nb+=b[i];}
 using Key=std::tuple<int,int,int>;std::map<Key,Z> dp;dp[Key(0,0,0)]=1;
 for(int i=0;i<a.size();++i){checkUserInterrupt();auto next=dp;for(const auto& z:dp){int r=std::get<0>(z.first),x=std::get<1>(z.first)+a[i],y=std::get<2>(z.first)+b[i];if(r<budget&&na-x>=2&&(paired_design||nb-y>=2))next[Key(r+1,x,y)]+=z.second;}dp=std::move(next);if(dp.size()>1000000)stop("Grouped plan exceeds one million states; reduce the block budget");}
 std::vector<Z> totals(budget+1);for(const auto& z:dp)totals[std::get<0>(z.first)]+=z.second;CharacterVector by(budget);Z total=0;for(int r=1;r<=budget;++r){by[r-1]=totals[r].str();total+=totals[r];}
 return List::create(_["feasible_deletions"]=total.str(),_["by_cardinality"]=by,_["units"]=a.size());
}
// [[Rcpp::export]]
std::string cpp_remaining(std::string planned,int completed){Z z(planned);z-=completed;if(z<0)stop("completed count exceeds plan");return z.str();}
