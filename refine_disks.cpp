// Refine already-counted critical-point disks of the same exact polynomial.
// The base certificate, not this program, supplies the exhaustive root count.
#define main base_verifier_main
#include "verify.cpp"
#undef main
#include <set>
int main(int argc,char**argv){try{
 std::cout.setf(std::ios::unitbuf);
 require(argc==4,"Usage: refine_disks parameters.txt critical_disks.txt indices.txt");
 auto begin=std::chrono::steady_clock::now();
 init_parameters(argv[1]);read_disks(argv[2]);
 require(rootDisks.size()==m*N-1,"Wrong containing-disk count");
 std::ifstream f(argv[3]);require(bool(f),"Cannot read refinement indices");
 std::set<unsigned> seen;unsigned j;
 while(f>>j){
  require(j<rootDisks.size()&&seen.insert(j).second,"Invalid or repeated refinement index");
  I rho,lhs,rhs,lower,upper;
  finite_disk(j,rho,lhs,rhs,lower,upper);
  require(lhs.hi<rhs.lo&&lower.lo>target.hi,"Refinement did not prove its target");
  std::cout<<"REFINED_DISK "<<j<<" radius="<<rho<<" lhs="<<lhs<<" rhs="<<rhs
           <<" quotient_lower="<<lower.lo<<" quotient_upper="<<upper.hi<<std::endl;
 }
 require(f.eof(),"Malformed index input");
 std::cout<<"PASS REFINEMENT degree="<<m*N<<" count="<<seen.size()<<" seconds="
          <<std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<<std::endl;
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<std::endl;return 1;}return 0;}
