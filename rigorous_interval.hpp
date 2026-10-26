#ifndef SMALE_RIGOROUS_INTERVAL_HPP
#define SMALE_RIGOROUS_INTERVAL_HPP
// Outward-rounded long-double intervals. Transcendental endpoints are
// computed by MPFR at 192 bits and rounded toward the appropriate infinity.
// This is NOT a use of unqualified libm transcendental rounding.
#include <gmp.h>
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cfenv>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// The installed runtime need not include the MPFR development header.
// These are the public MPFR 4.x ABI types and the small API subset used here.
#if __has_include(<mpfr.h>)
#include <mpfr.h>
#else
extern "C" {
 typedef long mpfr_prec_t;
 typedef long mpfr_exp_t;
 typedef int mpfr_sign_t;
 typedef enum { MPFR_RNDN=0,MPFR_RNDZ=1,MPFR_RNDU=2,MPFR_RNDD=3,MPFR_RNDA=4,MPFR_RNDF=5 } mpfr_rnd_t;
 typedef struct { mpfr_prec_t _mpfr_prec; mpfr_sign_t _mpfr_sign; mpfr_exp_t _mpfr_exp; mp_limb_t *_mpfr_d; } __mpfr_struct;
 typedef __mpfr_struct mpfr_t[1];
 typedef __mpfr_struct *mpfr_ptr;
 typedef const __mpfr_struct *mpfr_srcptr;
 void mpfr_init2(mpfr_ptr,mpfr_prec_t);
 void mpfr_clear(mpfr_ptr);
 int mpfr_set_ld(mpfr_ptr,long double,mpfr_rnd_t);
 int mpfr_set_q(mpfr_ptr,mpq_srcptr,mpfr_rnd_t);
 long double mpfr_get_ld(mpfr_srcptr,mpfr_rnd_t);
 int mpfr_sqrt(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
 int mpfr_log(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
 int mpfr_exp(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
 int mpfr_sin(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
 int mpfr_cos(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
 int mpfr_atan2(mpfr_ptr,mpfr_srcptr,mpfr_srcptr,mpfr_rnd_t);
 int mpfr_const_pi(mpfr_ptr,mpfr_rnd_t);
 const char *mpfr_get_version(void);
}
#endif

namespace rigorous {
using LD=long double;
constexpr LD inf=std::numeric_limits<LD>::infinity();
struct Inconclusive:std::exception {};
inline LD down(LD x){return std::nextafterl(x,-inf);}
inline LD up(LD x){return std::nextafterl(x,inf);}
struct Work {
 mpfr_t x,y,z;
 Work(){mpfr_init2(x,192);mpfr_init2(y,192);mpfr_init2(z,192);}
 ~Work(){mpfr_clear(x);mpfr_clear(y);mpfr_clear(z);}
};
inline Work& work(){static thread_local Work w;return w;}
using Fn=int(*)(mpfr_ptr,mpfr_srcptr,mpfr_rnd_t);
inline LD endpoint(Fn f,LD x,mpfr_rnd_t rnd){
 auto&w=work();mpfr_set_ld(w.x,x,MPFR_RNDN);f(w.z,w.x,rnd);return mpfr_get_ld(w.z,rnd);
}
inline LD angle_endpoint(LD y,LD x,mpfr_rnd_t rnd){
 auto&w=work();mpfr_set_ld(w.x,x,MPFR_RNDN);mpfr_set_ld(w.y,y,MPFR_RNDN);
 mpfr_atan2(w.z,w.y,w.x,rnd);return mpfr_get_ld(w.z,rnd);
}
struct I {
 LD lo,hi;
 I():lo(0),hi(0){}
 I(LD x):lo(x),hi(x){if(!std::isfinite(x))throw Inconclusive();}
 I(LD l,LD h):lo(l),hi(h){if(!std::isfinite(l)||!std::isfinite(h)||l>h)throw Inconclusive();}
 explicit I(const std::string&s){
  mpq_class q(s);q.canonicalize();auto&w=work();
  mpfr_set_q(w.z,q.get_mpq_t(),MPFR_RNDD);lo=mpfr_get_ld(w.z,MPFR_RNDD);
  mpfr_set_q(w.z,q.get_mpq_t(),MPFR_RNDU);hi=mpfr_get_ld(w.z,MPFR_RNDU);
 }
 bool finite()const{return std::isfinite(lo)&&std::isfinite(hi);}
 bool zero()const{return lo<=0&&hi>=0;}
 LD mid()const{return lo+(hi-lo)/2;}
 LD width()const{return up(hi-lo);}
};
inline I operator+(const I&a,const I&b){return I(down(a.lo+b.lo),up(a.hi+b.hi));}
inline I operator-(const I&a,const I&b){return I(down(a.lo-b.hi),up(a.hi-b.lo));}
inline I operator-(const I&a){return I(-a.hi,-a.lo);}
inline I operator*(const I&a,const I&b){
 std::array<LD,4>v={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
 return I(down(*std::min_element(v.begin(),v.end())),up(*std::max_element(v.begin(),v.end())));
}
inline I inv(const I&a){if(a.zero())throw Inconclusive();return I(down(1/a.hi),up(1/a.lo));}
inline I operator/(const I&a,const I&b){return a*inv(b);}
inline I& operator+=(I&a,const I&b){a=a+b;return a;}
inline I& operator-=(I&a,const I&b){a=a-b;return a;}
inline I& operator*=(I&a,const I&b){a=a*b;return a;}
inline I sqr(const I&a){
 if(a.zero())return I(0,up(std::max(a.lo*a.lo,a.hi*a.hi)));
 return I(down(std::min(a.lo*a.lo,a.hi*a.hi)),up(std::max(a.lo*a.lo,a.hi*a.hi)));
}
inline I abs(const I&a){if(a.zero())return I(0,std::max(-a.lo,a.hi));return a.hi<0?-a:a;}
inline I min(const I&a,const I&b){return I(std::min(a.lo,b.lo),std::min(a.hi,b.hi));}
inline I max(const I&a,const I&b){return I(std::max(a.lo,b.lo),std::max(a.hi,b.hi));}
inline I sqrt(const I&a){if(a.hi<0)throw Inconclusive();return I(endpoint(mpfr_sqrt,std::max(LD(0),a.lo),MPFR_RNDD),endpoint(mpfr_sqrt,std::max(LD(0),a.hi),MPFR_RNDU));}
inline I log(const I&a){if(a.lo<=0)throw Inconclusive();return I(endpoint(mpfr_log,a.lo,MPFR_RNDD),endpoint(mpfr_log,a.hi,MPFR_RNDU));}
inline I exp(const I&a){return I(endpoint(mpfr_exp,a.lo,MPFR_RNDD),endpoint(mpfr_exp,a.hi,MPFR_RNDU));}
inline I pow(const I&a,unsigned n){
 if(!n)return I(1);
 I v(1),b=a;
 while(n){if(n&1)v*=b;n>>=1;if(n)b*=b;}
 return v;
}
inline I pow(const I&a,const I&b){return exp(b*log(a));}
inline I pi(){auto&w=work();mpfr_const_pi(w.z,MPFR_RNDD);LD l=mpfr_get_ld(w.z,MPFR_RNDD);mpfr_const_pi(w.z,MPFR_RNDU);return I(l,mpfr_get_ld(w.z,MPFR_RNDU));}
inline bool overlap(const I&a,const I&b){return a.lo<=b.hi&&b.lo<=a.hi;}
inline I trig(const I&a,bool cosine){
 if(!a.finite()||a.width()>6.4L)return I(-1,1);
 Fn f=cosine?mpfr_cos:mpfr_sin;
 LD l=std::min(endpoint(f,a.lo,MPFR_RNDD),endpoint(f,a.hi,MPFR_RNDD));
 LD h=std::max(endpoint(f,a.lo,MPFR_RNDU),endpoint(f,a.hi,MPFR_RNDU));
 I p=pi();long k0=(long)std::floor(a.lo/3.14L)-4,k1=(long)std::ceil(a.hi/3.14L)+4;
 // The input ranges in this verifier are small. The conservative extra
 // indices cover all extrema even though 3.14 is used only for enumeration.
 if(k1-k0>1000)return I(-1,1);
 for(long k=k0;k<=k1;k++){
  I t=(I((LD)k)+(cosine?I(0):I(1)/I(2)))*p;
  if(overlap(a,t)){if((k%2)==0)h=1;else l=-1;}
 }
 return I(std::max(LD(-1),l),std::min(LD(1),h));
}
inline I sin(const I&a){return trig(a,false);}
inline I cos(const I&a){return trig(a,true);}
inline std::ostream&operator<<(std::ostream&o,const I&a){return o<<'['<<std::setprecision(22)<<a.lo<<','<<a.hi<<']';}
struct C {
 I re,im;
 C():re(0),im(0){}
 C(LD x):re(x),im(0){}
 C(const I&r):re(r),im(0){}
 C(const I&r,const I&i):re(r),im(i){}
 bool zero()const{return re.zero()&&im.zero();}
 C midpoint()const{return C(I(re.mid()),I(im.mid()));}
};
inline C operator+(const C&a,const C&b){return C(a.re+b.re,a.im+b.im);}
inline C operator-(const C&a,const C&b){return C(a.re-b.re,a.im-b.im);}
inline C operator-(const C&a){return C(-a.re,-a.im);}
inline C operator*(const C&a,const C&b){return C(a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re);}
inline C operator*(const C&a,const I&b){return C(a.re*b,a.im*b);}
inline C operator*(const I&a,const C&b){return b*a;}
inline C conj(const C&a){return C(a.re,-a.im);}
inline I norm(const C&a){return sqr(a.re)+sqr(a.im);}
inline I abs(const C&a){return sqrt(norm(a));}
inline C operator/(const C&a,const C&b){I d=norm(b);if(d.lo<=0)throw Inconclusive();return (a*conj(b))*inv(d);}
inline C operator/(const C&a,const I&b){return a*inv(b);}
inline C&operator+=(C&a,const C&b){a=a+b;return a;}
inline C&operator-=(C&a,const C&b){a=a-b;return a;}
inline C&operator*=(C&a,const C&b){a=a*b;return a;}
inline C&operator*=(C&a,const I&b){a=a*b;return a;}
inline I arg(const C&a){
 if(a.zero())throw Inconclusive();I p=pi();
 if(a.re.lo<0&&a.im.zero())return I(-p.hi,p.hi);
 LD l=inf,h=-inf;
 for(LD y:{a.im.lo,a.im.hi})for(LD x:{a.re.lo,a.re.hi}){
  l=std::min(l,angle_endpoint(y,x,MPFR_RNDD));h=std::max(h,angle_endpoint(y,x,MPFR_RNDU));
 }
 return I(l,h);
}
inline C exp(const C&a){I v=exp(a.re);return C(v*cos(a.im),v*sin(a.im));}
inline C log(const C&a){return C(log(abs(a)),arg(a));}
inline C pow(const C&a,const I&b){return exp(log(a)*b);}
inline C pow(const C&a,unsigned n){C v(1),b=a;while(n){if(n&1)v*=b;n>>=1;if(n)b*=b;}return v;}
inline C inflate(const C&a,LD radius){I e(-radius,radius);return C(a.re+e,a.im+e);}
inline LD box_radius(const C&a){return abs(a-a.midpoint()).hi;}
inline std::ostream&operator<<(std::ostream&o,const C&a){return o<<a.re<<" + i"<<a.im;}
inline void require(bool b,const std::string&why){if(!b)throw std::runtime_error(why);}
inline void init(){
 require(LDBL_MANT_DIG>=64&&LDBL_MANT_DIG<=192,"long double mantissa must have between 64 and 192 bits");
 require(std::numeric_limits<LD>::is_iec559,"IEEE 754 long double arithmetic required");
 require(std::numeric_limits<LD>::radix==2,"Binary floating point required");
 require(std::fesetround(FE_TONEAREST)==0,"Cannot set round-to-nearest");
 require(sizeof(mpfr_prec_t)==sizeof(long),"Unexpected MPFR ABI");
 I p=pi();require(p.lo<3.142L&&p.hi>3.141L&&p.width()<1e-17L,"MPFR ABI smoke test failed");
}
}
#endif
