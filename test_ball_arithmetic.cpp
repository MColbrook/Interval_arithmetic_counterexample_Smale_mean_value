// Exact-rational regression tests for the floating-center disk kernel.
// These tests supplement (not replace) its roundoff-error proof.
#define main unused_verifier_main
#include "verify.cpp"
#undef main
#include <random>
using Q=mpq_class;
struct GQ {Q re,im;};
Q exact(LD x){
 if(x==0)return Q(0);
 require(std::isfinite(x),"Nonfinite exact-conversion input");
 int e;LD f=std::frexp(x,&e);LD n=std::ldexp(f,LDBL_MANT_DIG);
 std::ostringstream out;out<<std::fixed<<std::setprecision(0)<<n;
 Q v(mpz_class(out.str()));
 if(e>=LDBL_MANT_DIG)v*=Q(mpz_class(1)<<(e-LDBL_MANT_DIG));
 else v/=Q(mpz_class(1)<<(LDBL_MANT_DIG-e));
 return v;
}
GQ add(GQ a,GQ b){return {a.re+b.re,a.im+b.im};}
GQ mul(GQ a,GQ b){return {a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re};}
bool contains(const DBall&d,GQ a){Q dx=a.re-exact(d.c.re.lo),dy=a.im-exact(d.c.im.lo),r=exact(d.radius);return dx*dx+dy*dy<=r*r;}
GQ point(const DBall&d,Q re,Q im){return {exact(d.c.re.lo)+exact(d.radius)*re,exact(d.c.im.lo)+exact(d.radius)*im};}
int main(){try{
 rigorous::init();std::mt19937_64 rng(0x534d414c45444953ULL);
 std::vector<std::pair<Q,Q>> dirs={{0,0},{1,0},{-1,0},{0,1},{0,-1},{Q(3,5),Q(4,5)},{Q(-3,5),Q(4,5)},{Q(3,5),Q(-4,5)},{Q(-3,5),Q(-4,5)}};
 unsigned long long checks=0;
 auto random=[&](int e){return std::ldexp((LD(rng())/std::ldexp(LD(1),63))-LD(1),e);};
 for(unsigned j=0;j<360;j++){
  int e=(int)(rng()%2001)-1000,ee=(int)(rng()%2001)-1000;
  if(j==0){e=-16430;ee=-16430;}if(j==1){e=8000;ee=8000;}if(j==2){e=16000;ee=-15990;}
  LD ar=random(e),ai=random(e),br=random(ee),bi=random(ee);
  C ac{I(ar),I(ai)},bc{I(br),I(bi)};
  LD ma=up(std::abs(ar)+std::abs(ai)),mb=up(std::abs(br)+std::abs(bi));
  LD ra=j%3?up(std::ldexp(ma,-48)):0,rb=j%4?up(std::ldexp(mb,-47)):0;
  DBall a{ac,ra,ma},b{bc,rb,mb};
  DBall sum=a+b,product=a*b;
  for(auto u:dirs)for(auto v:dirs){
   GQ x=point(a,u.first,u.second),y=point(b,v.first,v.second);
   require(contains(sum,add(x,y)),"Addition disk misses exact rational image");
   require(contains(product,mul(x,y)),"Product disk misses exact rational image");checks+=2;
  }
  I t=j%2?I("5/6"):I(-LD(0.25),LD(1.75));DBall scalar=a*t;
  for(auto u:dirs)for(LD te:{t.lo,t.hi}){
   GQ x=point(a,u.first,u.second);Q tt=exact(te);
   require(contains(scalar,{x.re*tt,x.im*tt}),"Scalar disk misses exact rational image");checks++;
  }
 }
 std::cout<<"PASS exact_rational_disk_image_checks="<<checks<<" significand_bits="<<LDBL_MANT_DIG<<" seed=0x534d414c45444953\n";
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}return 0;}
