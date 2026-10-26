#define main saved_main
#include "verify.cpp"
#undef main
int main(int argc,char**argv){try{
 std::cout.setf(std::ios::unitbuf); require(argc==4,"usage: check_base_comparison parameters disks index");
 init_parameters(argv[1]);read_disks(argv[2]);unsigned j=std::stoul(argv[3]);require(j<rootDisks.size(),"bad index");
 normalizationBase=a;Jets b=finite_P(a);baseValue=b.f;baseDerivative=b.fp;
 require(abs(baseDerivative).lo>0,"anchor derivative inconclusive");target=I(0);
 I rho,lhs,rhs,lo,hi;finite_disk(j,rho,lhs,rhs,lo,hi);
 require(hi.hi<1,"No rigorous upper bound below one");
 std::cout<<"OLD_BASE_FAILURE root="<<j<<" quotient_lower="<<lo<<" quotient_upper="<<hi<<" root_count=1 comparison=STRICT_PASS\n";
 std::cout<<"PASS BASE_COMPARISON degree="<<m*N<<"\n";
}catch(const std::exception&e){std::cerr<<"FAIL "<<e.what()<<std::endl;return 1;}return 0;}
