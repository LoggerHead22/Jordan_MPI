#include "jordan_mpi.h"


int main(int argc, char * argv[]){

    MPI_Init(&argc, &argv);
    int err=0, p, k, n, m;
    string name="";
    MPI_Comm G=MPI_COMM_WORLD;
    MPI_Comm_size(G,&p);
    MPI_Comm_rank(G,&k);

    if((argc>4) || (argc<3) || ((n=atoi(argv[1]))<=0) || ((m=atoi(argv[2]))<=0) ){
            if(!k) printf("Using %s n m <filename>\n",argv[0]);
            MPI_Finalize();
            return-1;
    }

    if(argc==4) name=argv[3];

    if(JordanSolvingSystem(n,m,name,p,k,G)){
            err=-2;
    }
  //  LOG("Fin");

    MPI_Finalize();
    return 0;
}
