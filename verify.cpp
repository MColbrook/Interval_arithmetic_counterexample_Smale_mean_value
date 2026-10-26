// Exhaustive critical-point disk verification with outward interval radii.
// The correction uses floating centers with explicitly proved roundoff bounds.
// Compile: g++ -O2 -std=c++17 -fno-fast-math -ffp-contract=off -frounding-math
//          verify_48000.cpp -o verify_48000 -l:libmpfr.so.6 -lgmpxx -lgmp
#include "rigorous_interval.hpp"
#include <chrono>
#include <map>
#include <sstream>
using namespace rigorous;
#ifndef N_PARAM
#define N_PARAM 3000
#endif
#ifndef BLOCKS_PARAM
#define BLOCKS_PARAM 4
#endif
#ifndef TAIL_BLOCKS_PARAM
#define TAIL_BLOCKS_PARAM 2
#endif
#ifndef EXTRA_BLOCKS_PARAM
#define EXTRA_BLOCKS_PARAM 0
#endif
constexpr unsigned extrablocks=EXTRA_BLOCKS_PARAM,extracount=16*extrablocks;
std::array<rigorous::C,extracount> extracoeff;
constexpr unsigned m=16,N=N_PARAM,L=N-2,blocks=BLOCKS_PARAM,K=m*blocks;
constexpr unsigned tailblocks=TAIL_BLOCKS_PARAM,tailcount=m*tailblocks;
std::array<C,tailcount> tailcoeff;
const C ii(I(0),I(1));
const I polynomialScale("6/5"),slopeScale("5/6"),twiceSlopeScale("5/3");
I A,CC,S,rm,eta,abeta,alam,alphaabs,Baabs,babs,haLo;
C a,r,b,beta,Wa,Ba,d0,d1,alpha,hc,hb0,hb1,haFinite,mu;
std::array<C,4> freec;
struct FiniteBall{C center;LD radius;};
std::array<std::vector<FiniteBall>,m-1> finiteCoefficients;
std::array<C,2> xi;
std::array<C,K> dc;
std::vector<C> correction;
std::array<I,15> first,acL,rooterr;
constexpr unsigned TO=8;
I BT;
std::vector<I> BTcoef;
I amplitude,lambda;
C tailPhase;
C normalizationBase,baseValue,baseDerivative;
const I rin("49/50"),rout("1999/2000"),anchorErr("1/100000000000000000000000000000000000000000000000000");
const I normErr("1/1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
I target("10000000001/10000000000");
I annulusTarget("5000001/5000000");
C Wfun(const C&z){C U=C(1)+ii*((pow(z,m)+ii*A)/S);return beta*z*pow(U,-I(1)/I(m));}
std::vector<C> mulpoly(const std::vector<C>&x,const std::vector<C>&y){std::vector<C> v(x.size()+y.size()-1);for(unsigned j=0;j<x.size();j++)for(unsigned k=0;k<y.size();k++)v[j+k]+=x[j]*y[k];return v;}
LD chooseLD(unsigned n,unsigned k){LD s=1;for(unsigned j=1;j<=k;j++)s=s*(n+1-j)/j;return s;}
I chooseI(unsigned n,unsigned k){I s(1);for(unsigned j=1;j<=k;j++)s=s*I(n+1-j)/I(j);return s;}
struct Jets{C f,fp,fpp;};
Jets high_tail(const C&z);
Jets finite_P(const C&z);
C highpow(const C&z,unsigned n);
Jets delta_factored(const C&z);
Jets free_delta(const C&z){
 C t=z-mu, u=z-a;
 Jets f{freec[3],C(0),C(0)};
 for(int j=2;j>=0;j--)f={f.f*t+freec[j],f.fp*t+f.f,I(2)*f.fp+f.fpp*t};
 return {u*u*f.f,I(2)*u*f.f+u*u*f.fp,I(2)*f.f+I(4)*u*f.fp+u*u*f.fpp};
}
// Evaluate the correction with its first two derivatives by centered Taylor
// expansion; all centre jets and the thirteenth-derivative bound are intervals.
using Jet=std::array<C,TO>;
Jet jetmul(const Jet&x,const Jet&y){
 Jet out{};for(unsigned j=0;j<TO;j++)for(unsigned k=0;k<=j;k++)out[j]+=x[k]*y[j-k];return out;
}
Jet jetlinear(Jet v,const C&constant,const I&slope=I(1)){
 for(int j=(int)TO-1;j>=0;j--)v[j]=v[j]*constant+(j?v[j-1]*slope:C(0));return v;
}
// Taylor coefficients are enclosed in the original compressed polynomial
// basis. The remainder is a Cauchy bound for the entire correction polynomial.
Jets delta_legacy(const C&z){
 C z0=z.midpoint(),h=z-z0;I rho=abs(h),R("1/100"),az=abs(z0);
 if(rho.hi>=R.lo)throw Inconclusive();
 Jet qj{};qj[0]=(pow(z0,m)+ii*A)/S;
 for(unsigned j=1;j<TO&&j<=m;j++)qj[j]=chooseI(m,j)*pow(z0,m-j)/S;
 Jet uj{};
 for(int ell=(int)blocks-1;ell>=0;ell--){
  Jet w{};
  for(int k=(int)m-1;k>=0;k--){w=jetlinear(w,z0/polynomialScale,I(5)/I(6));w[0]+=dc[ell*m+k];}
  uj=jetmul(uj,qj);for(unsigned j=0;j<TO;j++)uj[j]+=w[j];
 }
 Jet hj{};hj[0]=C(1);for(C c:{a,a,xi[0],xi[0],xi[1],xi[1]})hj=jetlinear(hj,z0-c);
 Jet kj{};kj[0]=C(1);for(C c:{xi[0],xi[0],xi[1],xi[1]})kj=jetlinear(kj,z0-c);
 Jet bj{};bj[0]=hb0+hb1*(z0-a);bj[1]=hb1;
 Jet herm=jetmul(bj,kj);herm[0]+=C(1);
 Jet js=jetmul(hj,uj);for(unsigned j=0;j<TO;j++)js[j]+=hc*herm[j];
 Jet fj{};
 for(int j=3;j>=0;j--){fj=jetlinear(fj,z0-mu);fj[0]+=freec[j];}
 fj=jetlinear(jetlinear(fj,z0-a),z0-a);
 for(unsigned j=0;j<TO;j++)js[j]+=fj[j];
 Jets center=delta_factored(z0);js[0]=center.f;js[1]=center.fp;js[2]=center.fpp/I(2);
 // For |w-z0|<=R, use the mean-value bound on q(w)-q(z0).
 I Z=az+R,Q=abs(qj[0])+I(m)*pow(Z,m-1)*R/S,U(0);
 for(int ell=(int)blocks-1;ell>=0;ell--){
  I w(0);for(int k=(int)m-1;k>=0;k--)w=w*Z/polynomialScale+abs(dc[ell*m+k]);U=U*Q+w;
 }
 I H(1),Kbound(1);for(C c:{a,a,xi[0],xi[0],xi[1],xi[1]})H*=abs(z0-c)+R;
 for(C c:{xi[0],xi[0],xi[1],xi[1]})Kbound*=abs(z0-c)+R;
 I M=H*U+abs(hc)*(I(1)+(abs(hb0)+abs(hb1)*(abs(z0-a)+R))*Kbound);
 I fbound(0);for(int j=3;j>=0;j--)fbound=fbound*(abs(z0-mu)+R)+abs(freec[j]);
 M+=sqr(abs(z0-a)+R)*fbound;
 I x=rho/R,om=I(1)-x;
 C f=js[TO-1],fp=js[TO-1]*I(TO-1),fpp=js[TO-1]*I((TO-1)*(TO-2));
 for(int j=(int)TO-2;j>=0;j--)f=f*h+js[j];
 for(int j=(int)TO-2;j>=1;j--)fp=fp*h+js[j]*I(j);
 for(int j=(int)TO-2;j>=2;j--)fpp=fpp*h+js[j]*I(j*(j-1));
 return {inflate(f,(M*pow(x,TO)/om).hi),
         inflate(fp,(I(TO)*M/R*pow(x,TO-1)/sqr(om)).hi),
         inflate(fpp,(I(TO*(TO-1))*M/sqr(R)*pow(x,TO-2)/pow(om,3u)).hi)};
}
// Factored correction jets for the two very small interior disks.
// Complex disks avoid exponential rectangular wrapping along the q-Horner chain.
// c is a point rectangle (a dyadic centre); radius is an outward upper bound.
// mag bounds |c|, and is exact-interval computed for the fixed multipliers.
struct DBall { C c; LD radius,mag; };
LD l1center(const C&c){return (abs(c.re)+abs(c.im)).hi;}
DBall frombox(const C&c){C mid=c.midpoint();return {mid,box_radius(c),abs(mid).hi};}
DBall ball_result(const C&box,const I&error){
 C mid=box.midpoint();C dev=box-mid;I enclosureRadius=abs(dev.re)+abs(dev.im);return {mid,(error+enclosureRadius).hi,l1center(mid)};
}
// Certified floating-center disk arithmetic. For IEEE round-to-nearest,
// |fl(t)-t| <= u*|t| + eta/2, with u=2^-p and eta=denorm_min.
// Complex center addition has L1 rounding error <= 4u*(|x|+|y|)+8eta.
// Four multiplications and two additions give a product-center L1 error
// <= 8u*|x|*|y|+32eta. These conservative constants hold for p>=64.
// The disk radii and all magnitude bounds below remain outward intervals.
// Overflow or non-finite results are rejected by interval constructors.
const I centerUnit(std::ldexp(LD(1),-LDBL_MANT_DIG));
// A minimum-NORMAL floor dominates denorm_min and avoids denormal assists.
const I centerTiny(std::numeric_limits<LD>::min());
DBall raw_ball(LD re,LD im,const I&rad){
 C center{I(re),I(im)};
 return {center,rad.hi,up(std::abs(re)+std::abs(im))};
}
DBall operator+(const DBall&x,const DBall&y){
 if(x.mag==0&&x.radius==0)return y;
 if(y.mag==0&&y.radius==0)return x;
 LD re=x.c.re.lo+y.c.re.lo,im=x.c.im.lo+y.c.im.lo;
 I rounding=I(4)*centerUnit*(I(x.mag)+I(y.mag))+I(8)*centerTiny;
 return raw_ball(re,im,I(x.radius)+I(y.radius)+rounding);
}
DBall operator*(const DBall&x,const DBall&y){
 if((x.mag==0&&x.radius==0)||(y.mag==0&&y.radius==0))return {C(0),0,0};
 LD ar=x.c.re.lo,ai=x.c.im.lo,br=y.c.re.lo,bi=y.c.im.lo;
 LD ac=ar*br,bd=ai*bi,ad=ar*bi,bc=ai*br;
 LD re=ac-bd,im=ad+bc;
 I rad=I(x.mag)*I(y.radius)+I(y.mag)*I(x.radius)+I(x.radius)*I(y.radius);
 I rounding=I(8)*centerUnit*I(x.mag)*I(y.mag)+I(32)*centerTiny;
 return raw_ball(re,im,rad+rounding);
}
DBall operator*(const DBall&x,const I&t){
 if(x.mag==0&&x.radius==0)return x;
 if(t.lo==0&&t.hi==0)return {C(0),0,0};
 LD c=t.mid();I delta=t.lo==t.hi?I(0):abs(t-I(c));
 LD re=x.c.re.lo*c,im=x.c.im.lo*c;
 I rad=x.radius==0?I(0):abs(I(c))*I(x.radius);
 if(delta.hi>0)rad+=(I(x.mag)+I(x.radius))*delta;
 I rounding=I(4)*centerUnit*I(x.mag)*abs(I(c))+I(8)*centerTiny;
 return raw_ball(re,im,rad+rounding);
}
C enclosure(const DBall&b){return inflate(b.c,b.radius);}
struct BJ {DBall f,fp,fpp;};
std::vector<DBall> dcballs,extraballs;
std::array<bool,extrablocks> extraactive;
std::array<DBall,4> freeballs;
Jets delta_factored(const C&z){
 const DBall zero{C(0),0,0},one{C(1),0,1};
 C q=(pow(z,m)+ii*A)/S,qp=pow(z,m-1)*I(m)/S,qpp=pow(z,m-2)*I(m*(m-1))/S;
 DBall qb=frombox(q),qpb=frombox(qp),qppb=frombox(qpp),zb=frombox(z/polynomialScale);
 BJ u{zero,zero,zero};
 for(int ell=(int)blocks-1;ell>=0;ell--){
  BJ w{zero,zero,zero};
  for(int k=(int)m-1;k>=0;k--){
   w={w.f*zb+dcballs[ell*m+k],w.fp*zb+w.f*slopeScale,w.fpp*zb+w.fp*twiceSlopeScale};
  }
  u={u.f*qb+w.f,u.fp*qb+u.f*qpb+w.fp,u.fpp*qb+(u.fp*qpb)*I(2)+u.f*qppb+w.fpp};
 }
 BJ h{one,zero,zero};
 for(C center:{a,a,xi[0],xi[0],xi[1],xi[1]}){
  DBall t=frombox(z-center);h={h.f*t,h.fp*t+h.f,h.fpp*t+h.fp*I(2)};
 }
 BJ kk{one,zero,zero};
 for(C center:{xi[0],xi[0],xi[1],xi[1]}){
  DBall t=frombox(z-center);kk={kk.f*t,kk.fp*t+kk.f,kk.fpp*t+kk.fp*I(2)};
 }
 DBall bb=frombox(hb0+hb1*(z-a)),b1=frombox(hb1),hcoef=frombox(hc);
 BJ freepart{zero,zero,zero};DBall t=frombox(z-mu),za=frombox(z-a);
 for(int j=3;j>=0;j--)freepart={freepart.f*t+freeballs[j],freepart.fp*t+freepart.f,freepart.fpp*t+freepart.fp*I(2)};
 for(unsigned j=0;j<2;j++)freepart={freepart.f*za,freepart.fp*za+freepart.f,freepart.fpp*za+freepart.fp*I(2)};
 BJ extra{zero,zero,zero};
 for(int ell=(int)extrablocks-1;ell>=0;ell--){
  BJ w{zero,zero,zero};
  if(extraactive[ell])for(int k=(int)m-1;k>=0;k--)w={w.f*zb+extraballs[ell*m+k],w.fp*zb+w.f*slopeScale,w.fpp*zb+w.fp*twiceSlopeScale};
  extra={extra.f*qb+w.f,extra.fp*qb+extra.f*qpb+w.fp,extra.fpp*qb+(extra.fp*qpb)*I(2)+extra.f*qppb+w.fpp};
 }
 for(unsigned j=0;j<2;j++)extra={extra.f*za,extra.fp*za+extra.f,extra.fpp*za+extra.fp*I(2)};
 freepart={freepart.f+extra.f,freepart.fp+extra.fp,freepart.fpp+extra.fpp};
 return {enclosure(h.f*u.f+hcoef*(one+bb*kk.f)+freepart.f),
         enclosure(h.fp*u.f+h.f*u.fp+hcoef*(b1*kk.f+bb*kk.fp)+freepart.fp),
         enclosure(h.fpp*u.f+(h.fp*u.fp)*I(2)+h.f*u.fpp+hcoef*((b1*kk.fp)*I(2)+bb*kk.fpp)+freepart.fpp)};
}

Jets delta(const C&z){ return delta_factored(z); }

struct Model{C z,q,qp,W,Wp,G,Gp,Gpp,Fp,Fpp,D,Dp,Dpp,E;I s;};
Model model_z(const C&z,bool factored=false,const C*knownq=nullptr,const I*knowns=nullptr,bool errorOnly=false){
 Model v;v.z=z;v.q=(pow(z,m)+ii*A)/S;v.qp=pow(z,m-1)*I(m)/S;v.s=abs(v.q);
 if(knownq){v.q=*knownq;v.s=*knowns;v.qp=I(m)*(v.q-ii*A/S)/z;}
 C U=C(1)+ii*v.q,V=U*CC;v.W=beta*z*pow(U,-I(1)/I(m));v.Wp=(beta/CC)*pow(U,-I(1)-I(1)/I(m));
 C den=C(1)-conj(r)*v.W;
 v.G=(v.W-r)/den/Ba;
 v.Gp=eta*v.Wp/pow(den,2u)/Ba;
 v.Gpp=v.Gp*(-ii*A*I(m+1)*v.qp*S/I(m)/V+conj(r)*v.Wp*I(2)/den);
 Jets dj=errorOnly?Jets{C(0),C(0),C(0)}:(factored?delta_factored(z):delta(z));v.D=dj.f;v.Dp=dj.fp;v.Dpp=dj.fpp;
 v.Fp=v.G+(z-a)*v.Gp+v.Dp;v.Fpp=v.Gp*I(2)+(z-a)*v.Gpp+v.Dpp;
 v.E=v.Fp/v.qp;return v;
}
Model model_polar(unsigned sheet,const I&t,const I&s){
 I theta=I(2)*pi()*t;
 C q=C(cos(theta),sin(theta))*s;
 C z=exp(ii*(I(2)*pi()*I(sheet)/I(m)+theta/I(m)))*pow(S*s,I(1)/I(m))*pow(C(1)-ii*A/(q*S),I(1)/I(m));
 Model v=model_z(z,false,&q,&s);return v;
}
// Lower bound for min_{0<=t<=1}|1+t*v|. Positivity also excludes
// the principal binomial branch cut (-infinity,-1] for every v in the box.
I segment_denominator(const C&v){
 C c=v.midpoint();I rho(box_radius(v)),d;
 if(c.re.lo>=0)d=I(1);
 else{
  I nv=norm(c),neg=-c.re;
  if(nv.hi<neg.lo)d=abs(C(1)+c);
  else if(nv.lo>neg.hi)d=abs(c.im)/sqrt(nv);
  else d=min(abs(C(1)+c),abs(c.im)/sqrt(nv));
 }
 I lower=d-rho;if(lower.lo<=0)throw Inconclusive();return lower;
}
struct Errors{I eh,ehp,eg,egp,eE;};
Errors errors(const Model&v,const I*theta=nullptr){
 I u=v.s,B,d,za,qpa;
 C vv=ii*v.q;
 if(theta){
  I sinth=sin(*theta),oneMinus=I(1)-sinth;oneMinus=I(std::max(LD(0),oneMinus.lo),oneMinus.hi);
  B=sqrt(sqr(I(1)-u)+I(2)*u*oneMinus);
  I tt=alphaabs*u,phase=*theta+arg(ii*alpha);I omc=I(1)-cos(phase);omc=I(std::max(LD(0),omc.lo),omc.hi);
  d=sqrt(sqr(I(1)-tt)+I(2)*tt*omc);
  za=pow(sqr(S*u-A)+I(2)*A*S*u*oneMinus,I(1)/I(2*m));
  qpa=I(m)*pow(za,m-1)/S;
 }else{B=abs(C(1)+vv);d=abs(C(1)-alpha*vv);za=abs(v.z);qpa=abs(v.qp);}
 if(B.lo<=0||d.lo<=0||qpa.lo<=0)throw Inconclusive();
 I t=alphaabs*u,tL=pow(t,L),uL=pow(u,L),uLp=uL*u;
 C J=-r*CC*(C(1)+vv)+pow(conj(r),15u)*(C(1)+ii*A)*(v.q*S-ii*A);
 C Jp=-r*CC*ii*v.qp+pow(conj(r),15u)*(C(1)+ii*A)*v.qp*S;
 I d0a=abs(d0),eh=abs(J)*tL/(d0a*d);
 I ehp=abs(Jp)*tL/(d0a*d)+abs(J)*alphaabs*qpa/d0a*(I(L)*pow(t,L-1)/d+tL/sqr(d));
 bool inside=false;I geom(1),s0(1),sigma(1);
 if(inside){geom=I(1)/(I(1)-u);s0=min(geom,I(2)/B);}
 else sigma=segment_denominator(vv);
 I zk=za,bk=babs,rk(1);
 for(unsigned k=1;k<m;k++){
  I gamma=I(m-k)/I(m);
  I tail,tailp;
  if(inside){
   I fac=min(geom,(I(1)+u*(I(1)+gamma)*s0/I(L+2))/B);
   I facp=min(geom,(I(1)+u*gamma*s0/I(L+1))/B);
   tail=uLp*(first[k-1]*fac+acL[k-1])/d;
   tailp=I(L+1)*uL*(first[k-1]*facp+acL[k-1])/d+alphaabs*tail/d;
  }else{
   // Exact Stieltjes remainder, not a convergent-series assertion outside |v|<1.
   I tail0=first[k-1]*uLp/sigma;
   I tail1=first[k-1]*(I(L+1)*uL/sigma
       +uLp*(I(L+1)-gamma)/I(L+2)/sqr(sigma));
   tail=(tail0+acL[k-1]*uLp)/d;
   tailp=(tail1+I(L+1)*acL[k-1]*uL)/d+alphaabs*tail/d;
  }
  I full=pow(B,gamma)/d,fullp=gamma*pow(B,gamma-I(1))/d+alphaabs*pow(B,gamma)/sqr(d);
  I ek=bk*tail+rooterr[k-1]*full,ekp=bk*tailp+rooterr[k-1]*fullp;
  I weight=eta*CC*rk/d0a;
  eh+=weight*zk*ek;
  ehp+=weight*(I(k)*zk/za*ek+zk*qpa*ekp);
  zk*=za;bk*=babs;rk*=rm;
 }
 I eg=eh/haLo+abs(v.G)*anchorErr/haLo;
 I egp=ehp/haLo+abs(v.Gp)*anchorErr/haLo;
 I eE=(eg+abs(v.z-a)*egp)/qpa;
 return {eh,ehp,eg,egp,eE};
}
I quotient2(const Model&v,const I&theta){
 I s=v.s;C vv=ii*v.q,den=C(1)-conj(r)*v.W;
 I om=I(1)-sin(theta);om=I(std::max(LD(0),om.lo),om.hi);
 I U2=sqr(I(1)-s)+I(2)*s*om;
 I w2=norm(v.W),psum(1),wk(1);for(unsigned k=1;k<m;k++){wk*=w2;psum+=wk;}
 I G2=(I(1)+eta*(sqr(S*s)-CC)/(sqr(CC)*U2*norm(den)*psum))/norm(Ba);
 C K=(v.D-v.q*v.E/I(N))/(v.z-a);
 return G2+I(2)*(conj(v.G)*K).re+norm(K);
}
// Centered first-order mean-value enclosure for the real-valued |T|^2.
// Full-box partial derivatives bound every segment in the polar rectangle.
I centered_quotient2(unsigned sheet,const I&t,const I&s,const Model&v){
 I tm(t.mid()),sm(s.mid()),theta0=I(2)*pi()*tm;
 Model v0=model_polar(sheet,tm,sm);
 I w0=quotient2(v0,theta0);
 C T=v.G+(v.D-v.q*v.E/I(N))/(v.z-a);
 C Ez=v.Fpp/v.qp-v.E*I(m-1)/v.z;
 C Tz=((I(1)-I(1)/I(N))*v.Fp-v.q*Ez/I(N)-T)/(v.z-a);
 C common=conj(T)*Tz*v.q/v.qp;
 I dtheta=I(2)*(ii*common).re;
 I ds=I(2)*common.re/s;
 return w0+dtheta*(I(2)*pi()*(t-tm))+ds*(s-sm);
}

// Exact finite h_L and its first two derivatives, enclosed at a point/box.
// The recurrence coefficients are intervals containing the exact rationals.
Jets finite_h_direct(const C&z){
 C q=(pow(z,m)+ii*A)/S,qp=I(m)*pow(z,m-1)/S,qpp=I(m*(m-1))*pow(z,m-2)/S;
 C v=ii*q,vp=ii*qp,vpp=ii*qpp,w=alpha*v,wp=alpha*vp,wpp=alpha*vpp;
 C den=C(1)-w,t=pow(w,L),tm1=pow(w,L-1),tm2=pow(w,L-2);
 C u=(C(1)-t)/den;
 C uw=-I(L)*tm1/den+(C(1)-t)/pow(den,2u);
 C uww=-I(L*(L-1))*tm2/den-I(2*L)*tm1/pow(den,2u)+I(2)*(C(1)-t)/pow(den,3u);
 C up=uw*wp,upp=uww*wp*wp+uw*wpp;
 C weight=pow(conj(r),m-1)*(C(1)+ii*A);
 C J=-r*CC*(C(1)+v)+weight*pow(z,m);
 C Jp=-r*CC*vp+weight*I(m)*pow(z,m-1);
 C Jpp=-r*CC*vpp+weight*I(m*(m-1))*pow(z,m-2);
 Jets out{J*u/d0,(Jp*u+J*up)/d0,(Jpp*u+I(2)*Jp*up+J*upp)/d0};
 C zk(1),rk(1),bk(1),factor=eta*CC/d0;
 for(unsigned k=1;k<m;k++){
  C vm=v.midpoint();I va=abs(vm);LD vr=box_radius(v);
  auto l1=[](const C&c)->LD{return (abs(c.re)+abs(c.im)).hi;};
  auto step=[&](const FiniteBall&x,const FiniteBall&add,unsigned scale)->FiniteBall{
   C cy=x.center*vm+I(scale)*add.center,mid=cy.midpoint();
   I rounding=abs(cy.re-mid.re)+abs(cy.im-mid.im);
   LD rad=(va*I(x.radius)+I(vr)*(I(l1(x.center))+I(x.radius))+I(scale)*I(add.radius)+rounding).hi;
   return {mid,rad};
  };
  FiniteBall bp=finiteCoefficients[k-1][L],bv{C(0),0},bvv{C(0),0};
  for(int j=(int)L-1;j>=0;j--){bvv=step(bvv,bv,2);bv=step(bv,bp,1);bp=step(bp,finiteCoefficients[k-1][j],1);}
  C p=inflate(bp.center,bp.radius),pv=inflate(bv.center,bv.radius),pvv=inflate(bvv.center,bvv.radius);
  zk*=z;bk*=b;C fk=factor*rk*bk*zk;
  out.f+=fk*p;
  out.fp+=fk*(I(k)*p/z+pv*vp);
  out.fpp+=fk*(I(k*(k-1))*p/(z*z)+I(2*k)*pv*vp/z+pvv*vp*vp+pv*vpp);
  rk*=conj(r);
 }
 return out;
}

// Exact finite beta-moment remainder, bounded by a finite geometric identity.
// This is not an infinite expansion outside its convergence disk.
constexpr unsigned MOMENT_TERMS=32;
Jets binomial_tail_moments(const C&v,unsigned k,const C&vLm1){
 I gamma=I(m-k)/I(m),bb=first[k-1],sig0=segment_denominator(v);
 LD rad=std::min((I(1)/I(L+1)).lo,(sig0/I(4)).lo);
 require(rad>0,"Positive moment Cauchy radius required");I R(rad);
 C wide=inflate(v,R.hi);I U=abs(wide),W=abs(C(1)+wide),sigma=segment_denominator(wide);
 std::array<I,MOMENT_TERMS+1> moment;moment[0]=I(1);
 I term=I(1)/sigma,best=term;unsigned count=0;
 for(unsigned j=1;j<=MOMENT_TERMS;j++){
  moment[j]=moment[j-1]*(gamma+I(j))/I(L+1+j);
  term=term*(U/W)*(gamma+I(j))/I(L+1+j);
  if(term.hi<best.hi){best=term;count=j;}
 }
 C inv=C(1)/(C(1)+v),x=v*inv,Q(0),Qp(0),Qpp(0);
 for(int j=int(count)-1;j>=0;j--){Qpp=Qpp*x+I(2)*Qp;Qp=Qp*x+Q;Q=Q*x+C(moment[unsigned(j)]);}
 C inv2=inv*inv,inv3=inv2*inv,inv4=inv3*inv,inv5=inv4*inv;
 C S0=inv*Q,S1=-inv2*Q+inv3*Qp,S2=I(2)*inv3*Q-I(4)*inv4*Qp+inv5*Qpp;
 I signed_b=(L%2==0)?bb:-bb;
 C r0=signed_b*vLm1*v*v*S0;
 C r1=signed_b*vLm1*(I(L+1)*v*S0+v*v*S1);
 C r2=signed_b*vLm1*(I(L)*I(L+1)*S0+I(2)*I(L+1)*v*S1+v*v*S2);
 I E=bb*pow(U,L+1)*best;
 return {inflate(r0,E.hi),inflate(r1,(E/R).hi),inflate(r2,(I(2)*E/sqr(R)).hi)};
}
Jets finite_h_moments(const C&z){
 std::array<C,m+1> zp;zp[0]=C(1);for(unsigned j=1;j<=m;j++)zp[j]=zp[j-1]*z;
 C q=(zp[m]+ii*A)/S,v=ii*q;
 C vp=ii*I(m)*zp[m-1]/S,vpp=ii*I(m*(m-1))*zp[m-2]/S;
 C av=alpha*v,den=C(1)-av,avLm2=highpow(av,L-2);
 C g0=(C(1)-avLm2*av*av)/den;
 C g1=(alpha*g0-I(L)*alpha*avLm2*av)/den;
 C g2=(I(2)*alpha*g1-I(L)*I(L-1)*alpha*alpha*avLm2)/den;
 C jw=pow(conj(r),m-1)*(C(1)+ii*A);
 C J=-r*CC*(C(1)+v)+jw*zp[m];
 C Jp=-r*CC*vp+jw*I(m)*zp[m-1];
 C Jpp=-r*CC*vpp+jw*I(m*(m-1))*zp[m-2];
 Jets out{J*g0/d0,(Jp*g0+J*g1*vp)/d0,(Jpp*g0+I(2)*Jp*g1*vp+J*(g2*vp*vp+g1*vpp))/d0};
 C vLm1=highpow(v,L-1),vL=vLm1*v,vLp1=vL*v;
 C root=pow(C(1)+v,-I(1)/I(m)),rootk(1),bpower(1),rpower(1);
 for(unsigned k=1;k<m;k++){
  rootk*=root;bpower*=b;if(k>1)rpower*=conj(r);
  I gamma=I(m-k)/I(m);C full=(C(1)+v)*rootk;
  Jets rem=binomial_tail_moments(v,k,vLm1);
  auto last=finiteCoefficients[k-1][L];C end=alpha*inflate(last.center,last.radius);
  C c0=(full-rem.f-end*vLp1)/den;
  C c1=(gamma*rootk-rem.fp-end*I(L+1)*vL+alpha*c0)/den;
  C c2=(gamma*(gamma-I(1))*rootk/(C(1)+v)-rem.fpp-end*I(L)*I(L+1)*vLm1+I(2)*alpha*c1)/den;
  C fac=eta*CC/d0*rpower*bpower;
  out.f+=fac*zp[k]*c0;
  out.fp+=fac*(I(k)*zp[k-1]*c0+zp[k]*c1*vp);
  C second=I(2*k)*zp[k-1]*c1*vp+zp[k]*(c2*vp*vp+c1*vpp);
  if(k>1)second+=I(k*(k-1))*zp[k-2]*c0;
  out.fpp+=fac*second;
 }
 return out;
}
Jets finite_h(const C&z){
 try{
  Jets ans=finite_h_moments(z);
  if(box_radius(ans.f)<LD(1e-12L)&&box_radius(ans.fp)<LD(1e-8L)&&box_radius(ans.fpp)<LD(1e-3L))return ans;
 }catch(const Inconclusive&){}
 return finite_h_direct(z);
}

void init_parameters(const std::string&path){
 init();std::ifstream file(path);require(bool(file),"Cannot open parameters");std::vector<std::string> lines;std::string str;while(std::getline(file,str))if(str.size())lines.push_back(str);require(lines.size()==26+2*K+2*tailcount+2*extracount,"Wrong parameter count");
 normalizationBase=C(I(lines[24+2*K+2*tailcount+2*extracount]),I(lines[25+2*K+2*tailcount+2*extracount]));
 A=I(lines[0]);r=C(I(lines[1]),I(lines[2]));a=C(I(lines[3]),I(lines[4]));b=C(I(lines[5]),I(lines[6]));
 xi[0]=C(I(lines[7]),I(lines[8]));xi[1]=C(I(lines[9]),I(lines[10]));for(unsigned j=0;j<K;j++)dc[j]=C(I(lines[11+2*j]),I(lines[12+2*j]));
 hc=C(I(lines[11+2*K]),I(lines[12+2*K]));
 for(unsigned j=0;j<4;j++)freec[j]=C(I(lines[14+2*K+2*j]),I(lines[15+2*K+2*j]));mu=(xi[0]+xi[1])/I(2);
 for(unsigned j=0;j<tailcount;j++)tailcoeff[j]=C(I(lines[22+2*K+2*j]),I(lines[23+2*K+2*j]));
 for(unsigned j=0;j<extracount;j++)extracoeff[j]=C(I(lines[22+2*K+2*tailcount+2*j]),I(lines[23+2*K+2*tailcount+2*j]));
 tailPhase=C(I(lines[22+2*K+2*tailcount+2*extracount]),I(lines[23+2*K+2*tailcount+2*extracount]));
 {mpq_class pr(lines[22+2*K+2*tailcount+2*extracount]),pi(lines[23+2*K+2*tailcount+2*extracount]);
  require(pr*pr+pi*pi==1,"Tail phase is not exactly unimodular");}
 extraballs.clear();for(auto coeff:extracoeff)extraballs.push_back(frombox(coeff));
 for(unsigned ell=0;ell<extrablocks;ell++){extraactive[ell]=false;for(unsigned k=0;k<m;k++){C c=extracoeff[ell*m+k];if(c.re.lo!=0||c.re.hi!=0||c.im.lo!=0||c.im.hi!=0)extraactive[ell]=true;}}
 dcballs.clear();for(auto coeff:dc)dcballs.push_back(frombox(coeff));for(unsigned j=0;j<4;j++)freeballs[j]=frombox(freec[j]);
 amplitude=I(lines[13+2*K]);require(amplitude.lo>0,"Amplitude must be positive");lambda=I(N)*amplitude;
 CC=I(1)+sqr(A);S=A+I(1)/A;require((S*rin-A).lo>0,"Contours must enclose all m branch roots");rm=abs(r);eta=I(1)-norm(r);require(eta.lo>0,"r not in unit disk");
 I eps("2/1000000000000000000000000000000000000000000000000000000000000");
 // Exact GMP rational root isolation, not floating residual arithmetic.
 using Q=mpq_class;struct CQ{Q re,im;};auto qm=[](CQ x,CQ y)->CQ{return {x.re*y.re-x.im*y.im,x.re*y.im+x.im*y.re};};
 Q aq(lines[0]),ccq=1+aq*aq,br(lines[5]),bi(lines[6]),eq("2/1000000000000000000000000000000000000000000000000000000000000");
 CQ bq{br,bi},bp{1,0};for(unsigned j=0;j<16;j++)bp=qm(bp,bq);
 Q res=abs(bp.re-Q(1)/ccq)+abs(bp.im-aq/ccq),bb=abs(br)+abs(bi),high=0;
 auto qp=[](Q x,unsigned n){Q p(1);for(unsigned j=0;j<n;j++)p*=x;return p;};
 for(unsigned k=2;k<=16;k++){mpz_class bin;mpz_bin_uiui(bin.get_mpz_t(),16,k);high+=Q(bin)*qp(bb,16-k)*qp(eq,k);}
 Q low=16*qp(br,15)*eq;require(res+high<low,"Exact rational beta-root test failed");
 require(br>eq&&bi>eq&&10*(bi+eq)<br-eq,"Cannot identify principal beta branch");
 beta=inflate(b,eps.hi);abeta=abs(beta);babs=abs(b);Wa=Wfun(a);Ba=(Wa-r)/(C(1)-conj(r)*Wa);Baabs=abs(Ba);haLo=Baabs-anchorErr;require(haLo.lo>0,"Invalid anchor normalization");
 d0=C(CC)+pow(conj(r),m)*C(-sqr(A),A);d1=ii*A-pow(conj(r),m)*(C(1)+ii*A);alpha=-d1*CC/(ii*A*d0);alphaabs=abs(alpha);
 require((alphaabs*rout).hi<1,"Rationalized pole may enter outer region");
 for(unsigned k=1;k<m;k++){
  I gamma=I(m-k)/I(m),bn(1);C cm(1);LD rad=0;
  finiteCoefficients[k-1].clear();finiteCoefficients[k-1].push_back({C(1),0});
  for(unsigned n=1;n<=L;n++){
   bn=bn*(gamma-I(n-1))/I(n);
   C cy=alpha*cm+C(bn),newmid=cy.midpoint();
   rad=(alphaabs*I(rad)+I(box_radius(cy))).hi;cm=newmid;finiteCoefficients[k-1].push_back({cm,rad});
  }
  C cball=inflate(cm,rad);acL[k-1]=abs(alpha*cball);first[k-1]=abs(bn*(gamma-I(L))/I(L+1));
  rooterr[k-1]=I(k)*eps*pow(abeta+eps,k-1);
 }
 C ka=pow(a-xi[0],2u)*pow(a-xi[1],2u);
 require(!ka.zero(),"Hermite normalization may vanish");
 hb0=-C(1)/ka;hb1=(I(2)/(a-xi[0])+I(2)/(a-xi[1]))/ka;
 require(K+5<m*N,"Correction degree must be below the leading power");
 require(extracount+1<m*N,"Extra correction degree must be below the leading power");
 haFinite=finite_h(a).f;require(abs(haFinite).lo>0,"Finite anchor may vanish");
 Model am=model_z(a);Errors ae=errors(am);require(ae.eh.hi<anchorErr.lo,"Anchor approximation error failed");require(abs(am.q).hi<I("1/4").lo,"Anchor q not small");
 Jets atail=high_tail(a);I actualNorm=abs(atail.fp),constErr=abs(atail.f);
 require(actualNorm.hi<normErr.lo&&constErr.hi<normErr.lo,"Anchor power error failed");
 Jets baseJets=finite_P(normalizationBase);baseValue=baseJets.f;baseDerivative=baseJets.fp;
 require(abs(baseDerivative).lo>0,"Normalization derivative may vanish");
 std::cout<<"NORMALIZATION base="<<normalizationBase<<" value="<<baseValue<<" derivative="<<baseDerivative<<std::endl;
 std::cout<<"ARITHMETIC MPFR="<<mpfr_get_version()<<" long_double_mantissa="<<LDBL_MANT_DIG<<" transcendental_bits=192\n";
 std::cout<<"INITIALIZATION beta_root=EXACT_PASS |alpha|="<<alphaabs<<" |Ba|="<<Baabs<<" |qa|="<<abs(am.q)<<" anchor_error="<<ae.eh<<" Taylor_order="<<TO<<" amplitude="<<amplitude<<" tail_phase="<<tailPhase<<" correction_enclosure=FLOATING_CENTER_DISKS roundoff_bound=EXPLICIT"<<"\n";
}
struct Box{unsigned j;LD tl,th,sl,sh;unsigned depth=0;};
struct Stats{unsigned long long tested=0,accepted=0,excluded=0,splits=0;unsigned maxdepth=0;LD minmargin=inf;};
std::vector<Box> initial(unsigned sheet,const I&s,unsigned nt){std::vector<Box> out;for(unsigned k=0;k<nt;k++)out.push_back({sheet,(LD)k/nt,(LD)(k+1)/nt,s.lo,s.hi,0});return out;}
void split(std::vector<Box>&stack,const Box&b,bool radial,Stats&st){
 if(b.depth>=55){
  std::cerr<<"FAILED_BOX sheet="<<b.j<<" t=["<<std::setprecision(22)<<b.tl<<","<<b.th<<"] s=["<<b.sl<<","<<b.sh<<"]\n";
  throw std::runtime_error("Subdivision limit reached");
 }Box x=b,y=b;x.depth++;y.depth++;if(radial){LD mid=b.sl+(b.sh-b.sl)/2;require(mid>b.sl&&mid<b.sh,"Radial precision exhausted");x.sh=mid;y.sl=mid;}else{LD mid=b.tl+(b.th-b.tl)/2;require(mid>b.tl&&mid<b.th,"Angular precision exhausted");x.th=mid;y.tl=mid;}stack.push_back(y);stack.push_back(x);st.splits++;
}
void statslog(const char*phase,unsigned sheet,const Stats&st){std::cout<<phase<<" sheet="<<sheet<<" tested="<<st.tested<<" accepted="<<st.accepted<<" excluded="<<st.excluded<<" splits="<<st.splits<<" maxdepth="<<st.maxdepth<<" min_margin="<<std::setprecision(16)<<st.minmargin<<std::endl;}
I boundary(unsigned sheet,bool outer){
 I radius=outer?rout:rin;std::vector<Box> todo=initial(sheet,radius,64);Stats st;I angleSum(0);
 while(!todo.empty()){
  Box b=todo.back();todo.pop_back();st.tested++;st.maxdepth=std::max(st.maxdepth,b.depth);bool ok=false;
  try{
   I t(b.tl,b.th),theta=I(2)*pi()*t;Model v=model_polar(sheet,t,radius);Errors e=errors(v,&theta);
   LD margin;
   if(outer){margin=(lambda*pow(radius,N-1)-abs(v.E)-e.eE).lo;ok=margin>0;}
   else{
    margin=(abs(v.Fp)-e.eE*abs(v.qp)-lambda*pow(radius,N-1)*abs(v.qp)).lo;
    if(margin>0&&!v.Fp.zero()){
     Model vl=model_polar(sheet,I(b.tl),radius),vr=model_polar(sheet,I(b.th),radius);C turn=vr.Fp*conj(vl.Fp);
     if(turn.re.lo>0){I da=arg(turn);angleSum+=da;ok=true;}
    }
   }
   if(ok){st.minmargin=std::min(st.minmargin,margin);st.accepted++;}
  }catch(const Inconclusive&){}
  if(!ok)split(todo,b,false,st);
  if(st.tested%20000==0)statslog(outer?"OUTER_PROGRESS":"INNER_PROGRESS",sheet,st);
 }
 statslog(outer?"OUTER":"INNER",sheet,st);if(!outer)std::cout<<"ANGLE sheet="<<sheet<<" interval="<<angleSum<<std::endl;return angleSum;
}
struct RootDisk{C center;I radius;};
std::vector<RootDisk> rootDisks;
void read_disks(const std::string&path){
 std::ifstream f(path);require(bool(f),"Cannot read root disks");std::string x,y,rr;
 while(f>>x){require(bool(f>>y>>rr),"Incomplete root disk");rootDisks.push_back({C(I(x),I(y)),I(rr)});}
 require(!rootDisks.empty()&&rootDisks.size()<m*N,"Invalid number of root disks");
}
void roots(){
 for(unsigned j=0;j<rootDisks.size();j++){
  C center=rootDisks[j].center;I rho=rootDisks[j].radius;
  require(rho.lo>0,"Root radius must be positive");
  C disk=inflate(center,rho.hi);Model vc=model_z(center,true),vd=model_z(disk,j<2);Errors e=errors(vd);
  I derivativeError=e.eg+abs(vd.z-a)*e.egp+lambda*pow(vd.s,N-1)*abs(vd.qp);
  I lhs=abs(vc.Fp)+rho*abs(vd.Fpp-vc.Fpp)+derivativeError,rhs=abs(vc.Fpp)*rho;
  if(!(lhs.hi<rhs.lo)) std::cout<<"ROOT_DIAGNOSTIC j="<<j<<" center_residual="<<abs(vc.Fp)<<" variation="<<abs(vd.Fpp-vc.Fpp)<<" derivative_error="<<derivativeError<<" lhs="<<lhs<<" rhs="<<rhs<<std::endl;
  require(lhs.hi<rhs.lo,"Interior disk Rouche failed");
  I mag=abs(vd.G+vd.D/(vd.z-a));
  I valueError=e.eg+(lambda/I(N)*pow(vd.s,N)+normErr)/abs(vd.z-a);
  I normalized=(mag-valueError)/(I(1)+normErr),upper=(mag+valueError)/(I(1)-normErr);
  require(normalized.lo>target.hi,"Interior quotient bound failed");
  require(vd.s.hi<rin.lo,"Interior disk not within inner region");
  std::cout<<"ROOT_DISK "<<j<<" radius="<<rho<<" Rouche_lhs="<<lhs<<" Rouche_rhs="<<rhs<<" quotient_lower="<<normalized.lo<<" quotient_upper="<<upper.hi<<"\n";
 }
 for(unsigned j=0;j<rootDisks.size();j++)for(unsigned k=0;k<j;k++)
  require(abs(rootDisks[j].center-rootDisks[k].center).lo>(rootDisks[j].radius+rootDisks[k].radius).hi,"Interior disks overlap");
}
void annulus(unsigned sheet){
 Stats st;std::vector<Box>todo=initial(sheet,I(rin.lo,rout.hi),64);
 while(!todo.empty()){
  Box b=todo.back();todo.pop_back();st.tested++;st.maxdepth=std::max(st.maxdepth,b.depth);bool ok=false,excluded=false;
  try{
   I t(b.tl,b.th),theta=I(2)*pi()*t;
   for(unsigned rep=0;rep<10;rep++){
    I s(b.sl,b.sh);Model v=model_polar(sheet,t,s);Errors e=errors(v,&theta);I eabs=abs(v.E);
    LD lower=std::max(LD(0),(eabs-e.eE).lo),upper=(eabs+e.eE).hi;
    if(upper<0)throw std::runtime_error("Negative magnitude bound");
    I ns;if(lower>0)ns=pow(I(lower,upper)/lambda,I(1)/I(N-1));else ns=I(0,pow(I(upper)/lambda,I(1)/I(N-1)).hi);
    LD sl=std::max(b.sl,ns.lo),sh=std::min(b.sh,ns.hi);
    if(sl>sh){excluded=true;break;}
    bool little=(sh-sl)>.92L*(b.sh-b.sl);b.sl=sl;b.sh=sh;
    if(little||rep==9){
     if(abs(v.z-a).lo<=I("1/100").hi)throw Inconclusive();
     I q2=quotient2(v,theta);
     I centered=centered_quotient2(sheet,t,I(b.sl,b.sh),v);
     q2=I(std::max(q2.lo,centered.lo),std::min(q2.hi,centered.hi));
     if(q2.lo<=0)throw Inconclusive();
     I tr=e.eg+v.s/(I(N)*abs(v.qp))*(e.eg/abs(v.z-a)+e.egp);
     I bound=(sqrt(q2)-tr-normErr/abs(v.z-a))/(I(1)+normErr);
     if(bound.lo>annulusTarget.hi){ok=true;st.minmargin=std::min(st.minmargin,bound.lo-I(1).hi);}break;
    }
   }
  }catch(const Inconclusive&){}
  if(excluded)st.excluded++;else if(ok)st.accepted++;else split(todo,b,(b.sh-b.sl)>3*(b.th-b.tl),st);
  if(st.tested%20000==0)statslog("ANNULUS_PROGRESS",sheet,st);
 }
 statslog("ANNULUS",sheet,st);
}

C highpow(const C&z,unsigned n){
 C center=z.midpoint(),rot(1);unsigned k=0;
 if(std::abs(center.re.mid())>=std::abs(center.im.mid())){if(center.re.mid()<0){rot=C(-1);k=2;}}
 else if(center.im.mid()>0){rot=-ii;k=3;}else{rot=ii;k=1;}
 C w=z*rot;require(w.re.lo>0,"High-power argument enclosure too wide");
 I theta=I(n)*arg(w),rad=pow(abs(z),I(n));C out=C(cos(theta),sin(theta))*rad;
 unsigned turn=(4-(k*n)%4)%4;for(unsigned j=0;j<turn;j++)out*=ii;
 return out;
}
// H=tailPhase*q^(N-J) [amplitude*q^J + sum t_ell,k*(z/(6/5))^k*q^ell].
Jets high_tail(const C&z){
 static_assert(tailblocks+2<N,"Invalid top-band width");
 C q=(pow(z,m)+ii*A)/S,qp=I(m)*pow(z,m-1)/S,qpp=I(m*(m-1))*pow(z,m-2)/S;
 Jets u{C(amplitude),C(0),C(0)};
 for(int ell=(int)tailblocks-1;ell>=0;ell--){
  Jets w{C(0),C(0),C(0)};
  for(int k=(int)m-1;k>=0;k--)
   w={w.f*z/polynomialScale+tailcoeff[ell*m+k],w.fp*z/polynomialScale+w.f/polynomialScale,w.fpp*z/polynomialScale+I(2)*w.fp/polynomialScale};
  u={u.f*q+w.f,u.fp*q+u.f*qp+w.fp,u.fpp*q+I(2)*u.fp*qp+u.f*qpp+w.fpp};
 }
 constexpr unsigned M=N-tailblocks;C t=tailPhase*highpow(q,M-2);
 return {t*q*q*u.f,t*(I(M)*q*qp*u.f+q*q*u.fp),
   t*(I(M)*(I(M-1)*qp*qp+q*qpp)*u.f+I(2*M)*q*qp*u.fp+q*q*u.fpp)};
}

I value_quotient2(const Model&v){
 I s=v.s,theta=arg(v.q);C den=C(1)-conj(r)*v.W;
 I om=I(1)-sin(theta);om=I(std::max(LD(0),om.lo),om.hi);
 I U2=sqr(I(1)-s)+I(2)*s*om;
 I w2=norm(v.W),psum(1),wk(1);for(unsigned k=1;k<m;k++){wk*=w2;psum+=wk;}
 I G2=(I(1)+eta*(sqr(S*s)-CC)/(sqr(CC)*U2*norm(den)*psum))/norm(Ba);
 C k=(v.D+high_tail(v.z).f)/(v.z-a);
 return G2+I(2)*(conj(v.G)*k).re+norm(k);
}
C full_first(const Model&v){return v.Fp+high_tail(v.z).fp;}
C full_second(const Model&v){return v.Fpp+high_tail(v.z).fpp;}

Jets finite_P(const C&z){
 Jets h=finite_h(z),d=delta_factored(z),tail=high_tail(z);
 C q=(pow(z,m)+ii*A)/S,qp=I(m)*pow(z,m-1)/S,qpp=I(m*(m-1))*pow(z,m-2)/S;
 return {(z-a)*h.f/haFinite+d.f+tail.f,
         (h.f+(z-a)*h.fp)/haFinite+d.fp+tail.fp,
         (I(2)*h.fp+(z-a)*h.fpp)/haFinite+d.fpp+tail.fpp};
}
// A finite-centre correction, used only when approximation tails would force
// a disk too large for the simpler analytic-model Rouche test.
void finite_disk(unsigned j,I&rho,I&lhs,I&rhs,I&lower,I&upper){
 C supplied=rootDisks[j].center,center=supplied;
 Jets exact=finite_P(center);
 for(unsigned it=0;it<5;it++){
  if(abs(exact.fpp).lo<=0)std::cout<<"FINITE_DIAGNOSTIC id="<<j<<" center="<<center<<" fp="<<exact.fp<<" fpp="<<exact.fpp<<" ha="<<haFinite<<std::endl;
  require(abs(exact.fpp).lo>0,"Finite Newton derivative inconclusive");
  C step=(exact.fp/exact.fpp).midpoint();
  if(abs(step).hi<LD(1e-16L))break;
  center=(center-step).midpoint();exact=finite_P(center);
 }
 Model vc=model_z(center,true);C modelSlope=full_second(vc);
 LD scale=(abs(exact.fp)/abs(exact.fpp)).hi;int exponent;
 std::cout<<"FINITE_START id="<<j<<" fp="<<abs(exact.fp)<<" fpp="<<abs(exact.fpp)<<" scale="<<I(scale)<<std::endl;
 // Start at the first dyadic radius above the residual scale.
 // Starting four times larger can skip the valid Rouche-radius window.
 std::frexp(std::max(scale,std::ldexp(LD(1),-60)),&exponent);rho=I(std::ldexp(LD(1),exponent));
 bool pass=false;
 for(unsigned attempt=0;attempt<11;attempt++){
  require((abs(center-supplied)+rho).hi<rootDisks[j].radius.lo,"Finite centre leaves its disjoint containing disk");
  Model vd=model_z(inflate(center,rho.hi));Errors ed=errors(vd);
  // Cauchy's second-derivative estimate for e=P'-g on a larger disk.
  I secondError(0);bool bounded=false;
  for(int er=-12;er>=-19;er--){
   I R(std::ldexp(LD(1),er));if(R.lo<=rho.hi*2)continue;
   try{
    Model wide=model_z(inflate(center,R.hi),true,nullptr,nullptr,true);Errors ew=errors(wide);
    I E=ew.eg+abs(wide.z-a)*ew.egp;
    I bound=I(2)*E/sqr(R-rho);
    if(!bounded||bound.hi<secondError.hi){secondError=bound;bounded=true;}
    if((rho*secondError).hi<(abs(exact.fpp)/I(1024)).lo)break;
   }catch(const Inconclusive&){}
  }
  if(!bounded)std::cout<<"CAUCHY_DIAGNOSTIC id="<<j<<" radius="<<rho<<" |q|="<<vc.s<<std::endl;
  require(bounded,"No analytic comparison disk for Cauchy estimate");
  I variation=abs(full_second(vd)-modelSlope)+rho*secondError;
  lhs=abs(exact.fp)+rho*variation;rhs=abs(exact.fpp)*rho;
  if(lhs.hi>=rhs.lo){rho*=I(2);continue;}
  // Finite-centre quotient at an independently chosen, exactly rational base.
  I separation=abs(vd.z-normalizationBase),norm=abs(baseDerivative);
  require(separation.lo>0,"Critical disk contains normalization base");
  I value0=abs((exact.f-baseValue)/((center-normalizationBase)*baseDerivative));
  C Q=((vd.z-a)*vd.G+vd.D+high_tail(vd.z).f-baseValue)/(vd.z-normalizationBase);
  I E1=ed.eg+abs(vd.z-a)*ed.egp,E0=abs(vd.z-a)*ed.eg;
  I valueVariation=rho*(abs(full_first(vd)-Q)+E1+E0/separation)/(separation*norm);
  lower=value0-valueVariation;
  upper=value0+valueVariation;
  require(lower.lo>target.hi,"Finite-centre critical quotient failed");
  pass=true;break;
 }
 require(pass,"Finite-centre Rouche failed");
 std::cout<<"FINITE_CENTER id="<<j<<" center="<<center<<" shift="<<abs(center-supplied)<<" radius="<<rho<<std::endl;
}
void allroots(int selected){
 require(rootDisks.size()==m*N-1,"Root list does not match the degree");
 // The supplied disks are disjoint; accepted disks may be strictly smaller.
 std::vector<unsigned> perm(rootDisks.size());for(unsigned j=0;j<perm.size();j++)perm[j]=j;
 std::sort(perm.begin(),perm.end(),[](unsigned j,unsigned k){return rootDisks[j].center.re.lo<rootDisks[k].center.re.lo;});
 LD maxrho=0;for(auto d:rootDisks){require(d.radius.lo>0,"Nonpositive root disk");maxrho=std::max(maxrho,d.radius.hi);}
 unsigned long long pairs=0;
 for(unsigned p=0;p<perm.size();p++)for(unsigned u=p+1;u<perm.size();u++){
  auto d=rootDisks[perm[p]],e=rootDisks[perm[u]];
  if(e.center.re.lo>up(d.center.re.hi+up(LD(2)*maxrho)))break;
  pairs++;
  require(abs(d.center-e.center).lo>(d.radius+e.radius).hi,"Supplied root disks overlap");
 }
 std::cout<<"DISJOINT count="<<rootDisks.size()<<" candidate_pairs="<<pairs<<" maximum_radius="<<I(maxrho)<<std::endl;
 unsigned done=0;LD minimum=inf,maximumUpper=0;
 for(unsigned j=0;j<rootDisks.size();j++){
  if(selected>=0&&j%16!=(unsigned)selected)continue;
  C center=rootDisks[j].center;I maxradius=rootDisks[j].radius;
  Model vc=model_z(center,true);C g0=full_first(vc),g1=full_second(vc);Errors ec=errors(vc);
  I err0=ec.eg+abs(vc.z-a)*ec.egp;
  require(abs(g1).lo>0,"Cannot bound centre derivative away from zero");
  LD scale=((abs(g0)+err0)/abs(g1)).hi;
  int expo;std::frexp(std::max(scale,std::ldexp(LD(1),-60)),&expo);
  I rho(std::ldexp(LD(1),expo+2));bool pass=false;I lastlhs,lastrhs,lower,upper;
  for(unsigned attempt=0;attempt<8;attempt++){
   if(rho.hi>maxradius.lo || rho.hi>LD(1e-6L))break;
   Model vd=model_z(inflate(center,rho.hi));Errors e=errors(vd);
   C gd=full_first(vd),gpd=full_second(vd);
   I error=e.eg+abs(vd.z-a)*e.egp;
   I lhs=abs(g0)+rho*abs(gpd-g1)+error,rhs=abs(g1)*rho;lastlhs=lhs;lastrhs=rhs;
   if(lhs.hi>=rhs.lo){rho*=I(2);continue;}
   I separation=abs(vd.z-normalizationBase),norm=abs(baseDerivative);
   require(separation.lo>0,"Critical disk contains normalization base");
   C Tc=(vc.z-a)*vc.G+vc.D+high_tail(vc.z).f;
   I value0=abs((Tc-baseValue)/((vc.z-normalizationBase)*baseDerivative));
   C Q=((vd.z-a)*vd.G+vd.D+high_tail(vd.z).f-baseValue)/(vd.z-normalizationBase);
   I variation=rho*abs((gd-Q)/((vd.z-normalizationBase)*baseDerivative));
   I valueError=abs(vd.z-a)*e.eg/(separation*norm);
   lower=value0-variation-valueError;
   upper=value0+variation+valueError;
   if(lower.lo<=target.hi){
    std::cout<<"MODEL_VALUE_INCONCLUSIVE id="<<j<<" center="<<center<<" rho="<<rho<<" value0="<<value0<<" variation="<<variation<<" valueError="<<valueError<<" lower="<<lower<<std::endl;
    break; // Try a finite-centre enclosure instead of rejecting the candidate.
   }
   pass=true;break;
  }
  if(!pass){finite_disk(j,rho,lastlhs,lastrhs,lower,upper);pass=true;}
  if(!pass){std::cout<<"ROOT_ROUCHE_FAIL id="<<j<<" center="<<center<<" rho="<<rho<<" maxradius="<<maxradius<<" g0="<<abs(g0)<<" error0="<<err0<<" g1="<<abs(g1)<<" lhs="<<lastlhs<<" rhs="<<lastrhs<<std::endl;throw std::runtime_error("Critical root isolation failed");}
  std::cout<<"CRITICAL_DISK "<<j<<" radius="<<rho<<" lhs="<<lastlhs<<" rhs="<<lastrhs<<" quotient_lower="<<lower.lo<<" quotient_upper="<<upper.hi<<std::endl;
  if(lower.lo<minimum){minimum=lower.lo;maximumUpper=upper.hi;}
  done++;
  if(done%1000==0)std::cout<<"ROOT_PROGRESS done="<<done<<" lower="<<I(minimum)<<std::endl;
 }
 std::cout<<"ALL_ROOTS checked="<<done<<" total="<<rootDisks.size()<<" minimum_lower="<<I(minimum)<<" associated_upper="<<I(maximumUpper)<<std::endl;
}

void probe(int selected){
 for(unsigned j=0;j<rootDisks.size();j++){
  if(selected>=0&&j%16!=(unsigned)selected)continue;
  C z=rootDisks[j].center,q=(pow(z,m)+ii*A)/S;
  if(abs(q).hi<.9988L)continue;
  for(unsigned it=0;it<8;it++){
   Jets ex=finite_P(z);require(abs(ex.fpp).lo>0,"Probe slope may vanish");
   C step=(ex.fp/ex.fpp).midpoint();z=(z-step).midpoint();
   if(abs(step).hi<1e-17L)break;
  }
  std::cout<<"PROBE "<<j<<" "<<std::setprecision(23)<<z.re.mid()<<" "<<z.im.mid()<<"\n";
 }
}

int main(int argc,char**argv){try{
 std::cout.setf(std::ios::unitbuf);auto start=std::chrono::steady_clock::now();
 std::string mode=argc>1?argv[1]:"all";int selected=argc>2?std::stoi(argv[2]):-1;
 init_parameters(argc>3?argv[3]:"parameters.txt");
 if(mode=="init")return 0;
 require(mode=="allroots"||mode=="disk"||mode=="probe"||mode=="compare","Only allroots, disk, probe, compare, and init are supported");
 read_disks(argc>4?argv[4]:"root_disks.txt");
 if(mode=="disk"){
  require(selected>=0&&(unsigned)selected<rootDisks.size(),"Invalid disk index");
  I rho,lhs,rhs,lower,upper;finite_disk(selected,rho,lhs,rhs,lower,upper);
  std::cout<<"SINGLE_DISK id="<<selected<<" radius="<<rho<<" lhs="<<lhs<<" rhs="<<rhs<<" lower="<<lower<<" upper="<<upper<<std::endl;
 }
 if(mode=="allroots")allroots(selected);
 if(mode=="probe")probe(selected);
 if(mode=="compare"){
  require(selected>=0&&(unsigned)selected<rootDisks.size(),"Bad comparison index");
  C z=rootDisks[selected].center;Jets fast=finite_h_moments(z),direct=finite_h_direct(z);
  for(auto pair:std::array<std::pair<C,C>,3>{{{fast.f,direct.f},{fast.fp,direct.fp},{fast.fpp,direct.fpp}}}){
   require(overlap(pair.first.re,pair.second.re)&&overlap(pair.first.im,pair.second.im),"Finite evaluators disagree");
   std::cout<<"FINITE_COMPARE fast="<<pair.first<<" direct="<<pair.second<<std::endl;
  }
 }
 if(mode=="roots"||mode=="all")roots();
 I winding(0);
 if(mode=="inner"||mode=="all"){for(unsigned j=0;j<m;j++)if(selected<0||selected==(int)j)winding+=boundary(j,false);if(selected<0){I w=winding/(I(2)*pi());require(w.lo>(LD)rootDisks.size()-.5L&&w.hi<(LD)rootDisks.size()+.5L,"Inner winding does not match the disk count");std::cout<<"INNER_WINDING="<<w<<" integer="<<rootDisks.size()<<"\n";}}
 if(mode=="outer"||mode=="all")for(unsigned j=0;j<m;j++)if(selected<0||selected==(int)j)boundary(j,true);
 if(mode=="annulus"||mode=="all")for(unsigned j=0;j<m;j++)if(selected<0||selected==(int)j)annulus(j);
 double secs=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cout<<"PASS mode="<<mode<<" sheet="<<selected<<" seconds="<<secs<<" degree="<<16*N<<" target="<<target<<"\n";
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<"\n";return 1;}return 0;}
