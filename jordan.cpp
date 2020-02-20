#include "jordan_mpi.h"



double Random() {
    return rand() / double(RAND_MAX);
}

void generate_file(int n) {
    ofstream file("data.txt");
    for (int i = 0; i < n * n; i++) {
        if (i % n == 0 && i != 0) {
            file << "\n";
        }
        file << 5 + Random() * 10 << "\t";
    }
}

void print_matrix(double* m, int size) {
    for (int y = 0; y < ((size < 10) ? size : 10) ; y++) {
        for (int x = 0; x < ((size < 10) ? size : 10); x++) {
            cout << m[y * size + x] << "\t";
        }
        cout << endl;
    }
}

void new_print_matrix(double* a, int n, int m) {
    int k=n/m;
    int N=10;
    for (int y = 0; y < ((n < N) ? n : N) ; y++) {
        for (int x = 0; x < ((n < N) ? n : N); x++) {
            int k1 = x /m ; // row
            int l = (k1 == k ? n%m :m);

            cout << a[k1*n*m + y*l + x%m] << "\t";
        }
        cout << endl;
    }
}

void print_equa(double* a, double* b, int size) {
    for (int y = 0; y < ((size < 10) ? size : 10); y++) {
        for (int x = 0; x < ((size < 10) ? size : 10); x++) {
            cout << a[y * size + x] << "\t";
        }
        cout << " | " << b[y] << endl;
    }
}

void print_vector(int* m, int size) {
    int N =15;
    for (int y = 0; y  <((size < N) ? size : N); y++) {
        cout << m[y] << "  ";
    }
    cout << endl;
}

void print_vector(double* m, int size) {
    int N = 15;
    for (int y = 0; y < ((size < N) ? size : N); y++) {
        cout << m[y] << "  ";
    }
    cout << endl;
}

void print_matrix_rect(double* m, int w, int h) {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            if (abs(m[y * w + x]) < 1e-7) {
                cout << 0 << "\t";
            }
            else {
                cout << m[y * w + x] << "\t";
            }
        }
        cout << endl;
    }
}

void matr_to_E(double* m, int n,int blockSize) {
    int k = n / blockSize;
    for (int y = 0; y < n; y++) {
        for (int x = 0; x < n; x++) {
            int k1 = x / blockSize ; // row
            int l = (k1 == k ? n % blockSize : blockSize);

            m[k1*n*blockSize + y*l + x%blockSize] = x == y;
        }
    }
}

void matr_to_NULL(double* m, int size) {
    for (int y = 0; y < size * size; y++) {
        m[y] = 0;
    }
}

double formula_matr_local(int i, int j){
  //  return int(-5 + Random() * 10);
    return fabs(i-j);
   // return 1./(1 + i + j);
}

void formula_matr(double* a, int n,int m,int proc_num, int p, double *b, double &Norma, MPI_Comm G) {
    int k=n/m;
    double tempNorma = 0, buf = 0;
    for (int i = 0; i < n; i++) {
        if(proc_num==0 || proc_num==k%p){
            b[i]=0;
        }
        for (int j = 0; j < n; j++) {
          //  buf = fabs(i - j);

            buf = formula_matr_local(i,j);
            tempNorma+=buf;
            int k1 = j /m ; // row
            int l = (k1 == k ? n%m :m);
            if(proc_num ==0 && i < 10 && j < 10){
                cout<<round(buf*1000)/1000<< "\t";
            }
            //a[k1*n*m + i*l + j%m] = fabs(i-j);
            // p[k] = abs(i - k);
            if((j % (p*m)) / m ==proc_num ){
                a[(k1/p)*n*m + i*l + j%m] = buf;
            }
            if(j%2 == 0 && (proc_num == k%p || proc_num==0) ){
                b[i]+= buf;
            }

        //    a[k1*n*m + i*l + j%m] = (double) 1/( i + j + 1);
            //a[k1*n*m + i*l + j%m] = int(-5 + Random() * 10);
            //          p[k] = (double) 1/( i + k + 1);
            //p[k] = n - max(i,k);
        }
        if(proc_num ==0 && i < 10 ){
            cout<< " ... = "<<b[i]<<endl;
        }
        if(tempNorma > Norma){
            Norma = tempNorma;
        }
        tempNorma = 0;
    }
}



void initMatrix(double *a, int n, int m,int p,int thrnum){
    int c=n/m;
    for(int k=0;k<n;k++){
        for(int i=thrnum*m;i < n;i++){
            int k1 = i /m ; // row
            int l = (k1 == c ? n%m :m);
            a[k1*n*m + k*l + i%m]=thrnum;
            if(i%m==m-1) i+=(p-1)*m;
        }
    }

}

int JordanInv(double* a, int size, double* b, double norma) {
    matr_to_E(b, size, size);
    //print_matrix(a,size);
    double tempNorm=norma_matr(a,size,size);
    //LOG(tempNorm);
    if(norma > tempNorm && tempNorm > EPS){
        norma = tempNorm;
    }
    if(norma > 1){
        norma = 1;
    }

    double *p1,*p2;
    int* columns = new int[size];
    for (int i = 0; i < size; i++) {
        columns[i] = i;
    }
    //	LOG("Nachinaen Jordan");
    for (int q = 0; q < size; q++) {
        // ищем элемент с максимальной нормой
        p1 = a + q*size;
        p2 = b + q*size;
        int index = q;
        for (int j = q + 1; j < size; j++) {
            if (fabs(p1[columns[j]]) > fabs(p1[columns[index]])) {
                index = j;
            }
        }
        if (fabs(p1[columns[index]]) < EPS * norma) {
            delete[] columns;
            return ERR_DEG_MATRIX;
        }
        swap(columns[q], columns[index]);

        // умножим строку на 1 / a_qq
        double k = 1 / a[q * size + columns[q]];
        p1[columns[q]]=1;
        for (int j = q+1; j < size; j++) {
            p1[columns[j]] *= k;
        }
        for (int j = 0; j < size; j++) {
            p2[j] *= k;
        }

        // вычтем
        for (int i = 0; i < size; i++) {
            if (i == q) {
                continue;
            }

            double k = a[i * size + columns[q]];
            for (int j = q ; j < size; j++) {
                a[i * size + columns[j]] -= p1[columns[j]] * k;
            }
            for (int j = 0; j < size; j++) {
                b[i * size + j] -=p2[j] * k;
            }
        }
    }

    for (int i = 0; i < size; i++) {
        p1 = a + columns[i] * size;
        p2 = b + i*size;
        for (int j = 0; j < size; j++) {
            p1[j] = p2[j];
        }
    }

    delete[] columns;
    return 0;
}

//int read_matrix(double* a, int n,int m, const string& name, int proc_num, int p) {
//    ifstream file(name);
//    if (!file.is_open()) {
//        printf("read_matrix: can\'t open the file %s \n", name.c_str());
//        return ERR_CANNOT_OPEN;
//    }
//    int k = n / m;
//    for (int i = 0; i < n; i++) {
//        for(int j=0;j<n;j++){
//            int k1 = j /m ; // row
//            int l = (k1 == k ? n%m :m);
//            if((j % (p*m)) % m ==proc_num ){
//                LOG(k1/p + i*l + j%m);
//                if (!(file >> a[k1/p + i*l + j%m]) ) {
//                    return ERR_READ;
//                }
//            }
//        }
//    }
//    return 0;
//}

int read_matrix(double* a, int n,int m,  const string& name, int proc_num, int p, double *b, double &Norma, double *CopyCol, MPI_Comm G) {
    int k = n / m;
    double tempNorma = 0, buf = 0;
    int loc_err = 0, glob_err = 0;
    ifstream file;
    if (proc_num == 0) {
        file.open(name.c_str());
        if (!file.is_open()) {
            printf("read_matrix: can\'t open the file %s \n", name.c_str());
            loc_err = -1;
        }
    }

    MPI_Allreduce(&loc_err,&glob_err, 1, MPI_INT, MPI_SUM, G);
    if(glob_err < 0){
        return ERR_CANNOT_OPEN;
    }

    for(int l = 0; l < k + 1;l++){
        int end_l = (l==k ? n % m : m);
        if(proc_num==0){
            for(int j=0 ; j< end_l; j++){ //y
                   b[l*m + j] = 0;
                for (int i = 0; i < n; i++) { //x
                    if (!(file >> buf) ) {
                        loc_err = -1;
                    }else{
                        if(i < 10 && l*m + j < 10){
                            cout<<round(buf*1000)/1000<<"\t";
                        }

                        if(i % 2 == 0){
                            b[l*m + j]+=buf;
                        }
                        tempNorma+=buf;
                        CopyCol[j*n + i] = buf;
                    }
                }
                if(tempNorma > Norma){
                    Norma = tempNorma;
                }
                tempNorma=0;
                if(l*m + j < 10){ cout<<endl;}
            }
        }

        MPI_Allreduce(&loc_err,&glob_err, 1, MPI_INT, MPI_SUM, G);
        if(glob_err < 0){
            return ERR_READ;
        }


        MPI_Bcast(CopyCol, n*end_l, MPI_DOUBLE, 0, G);
        for(int j=l*m;j<l*m + end_l;j++){
            for (int i = 0; i < n; i++) {

                int k1 =i /m ; // row
                int l2 = (k1 == k ? n%m :m);
                if((i % (p*m)) / m ==proc_num ){
         //           LOG((k1/p)*n*m + j*l2 + i%m);
                    a[(k1/p)*n*m + j*l2 + i%m] = CopyCol[(j - l*m)*n + i];
                }
            }
        }

            //        for( int y = 0; y < end_l; y++){
            //            for( int x = 0; x < n ; x++){
            //                int end2 = (x / m == k ? n % m : m);
            //                int c = (x / m) / p;
            //                if((x / m) % p == proc_num ){
            //                    a[c*n*m + l*m*end_l + x%m + y*end2] = CopyCol[y*n + x];
            //                }
            //            }
            //        }
        }
    MPI_Status s;
    if(proc_num == 0 && k % p != 0){
        MPI_Send(b ,n , MPI_DOUBLE, k%p, 0,G);
    }else if(proc_num != 0 && proc_num == k % p){
        MPI_Recv(b ,n , MPI_DOUBLE, 0, 0,G, &s);
    }

    MPI_Bcast(&Norma, 1, MPI_DOUBLE, 0, G);
    if(proc_num==0){
        printf("Norma: %lf \n",Norma);
    }
    return 0;
}

int init_matrix_file(double* a, int n, int m, const string& name, int proc_num, int p, double *b, double &Norma, double *CopyCol, MPI_Comm G) {
    int res =-1;
    if (name != "") {
        int res = read_matrix(a, n, m, name, proc_num,p,b,Norma,CopyCol,G);
        if (res < 0 && proc_num == 0) {
            switch (res) {
            case ERR_CANNOT_OPEN: {
                printf("Cannot open %s \n", name.c_str());
                break;
            }
            case ERR_READ: {
                printf("Invalid data in %s \n", name.c_str());
                break;
            }
            default: {
                printf("Unknown error %d in %s n", res, name.c_str());
                break;
            }
            }
        }
        return res;
    }
    return res;
}

double norma_matr(double* a, int n, int m) {
    int k = n / m;
    double ret = -1;
    for (int y = 0; y < n; y++) {
        double sum = 0;
        for (int x = 0; x < n; x++) {
            int k1 = x / m ; // row
            int l = (k1 == k ? n%m :m);
            sum += fabs(a[k1*n*m + y*l + x%m]);
        }
        if (sum > ret) {
            ret = sum;
        }
    }
    return ret;
}

void get_block(int proc_num, int p, double* a, int aSize, int blockSize, int x, int y, double* c3, int count, bool isNumaColumn) {
    int k=aSize/blockSize;
    int step=blockSize*blockSize;
    if(x==(k - k%p)/p && proc_num == k % p && !isNumaColumn){
        step = (aSize%blockSize)*blockSize;
    }
    memcpy(c3,a + x*blockSize*aSize + y*step ,count * sizeof(double));
}

void push_block(int proc_num, int p, double* a, int aSize, int blockSize, int x, int y, double* c3, int count, bool isNumaColumn) {
    int k=aSize/blockSize;
    int step=blockSize*blockSize;
    if(x==(k - k%p)/p && proc_num==k%p && !isNumaColumn){
        step = (aSize%blockSize)*blockSize;
    }
    memcpy(a + x*blockSize*aSize + y*step,c3,count * sizeof(double));
}

void get_b_block(double* a, int blockSize, int y, double* c3, int count) {
    memcpy(c3,a  + y*blockSize,count * sizeof(double));
}

void push_b_block(double* a, int blockSize, int y, double* c3, int count) {
    memcpy(a  + y*blockSize,c3,count * sizeof(double));
}

void mult_matrix(double *a, double *b, double *c, int n, int m, int k){
    int i,j,l;
    double *pc=c,*pa=a,*pb=b, sum[9];

    for(i=0;i<n;i++)
        for(j=0;j<k;j++) c[i*n+j]=0.;

    if(n >=3 && k>=3){
        for(i=0;i<n-n%3;i+=3) {
            for(j=0;j<k-k%3;j+=3){
                sum[0]=sum[1]=sum[2]=sum[3]=sum[4]=sum[5]=sum[6]=sum[7]=sum[8]=0.;
                for(l=0;l<m;l++){
                    pa=a+i*m+l;
                    pb=b+l*k+j;
                    sum[0]+=pa[0]*pb[0];
                    sum[1]+=pa[0]*pb[1];
                    sum[2]+=pa[0]*pb[2];
                    sum[3]+=pa[m]*pb[0];
                    sum[4]+=pa[m]*pb[1];
                    sum[5]+=pa[m]*pb[2];
                    sum[6]+=pa[2*m]*pb[0];
                    sum[7]+=pa[2*m]*pb[1];
                    sum[8]+=pa[2*m]*pb[2];
                }
                pc=c+i*k+j;
                pc[0]			=sum[0];
                pc[1]			=sum[1];
                pc[2]			=sum[2];
                pc[k]			=sum[3];
                pc[k+1]			=sum[4];
                pc[k+2]			=sum[5];
                pc[2*k]			=sum[6];
                pc[2*k+1]		=sum[7];
                pc[2*k+2]		=sum[8];
            }
        }
    }
    sum[0]=sum[1]=sum[2]=0;
    for(i=n - n%3 ;i<n;i++){
        for(j=0;j<k;j++){
            pc=c+i*k+j;
            pc[0]=0;
            for(l=0;l<m;l++){
                pa=a+i*m+l;
                pb=b+l*k+j;
                pc[0]+=pa[0]*pb[0];
            }
        }
    }
    for(i=0 ;i<n-n%3;i++){
        for(j=k - k%3;j<k;j++){
            pc=c+i*k+j;
            pc[0]=0;
            for(l=0;l<m;l++){
                pa=a+i*m+l;
                pb=b+l*k+j;
                pc[0]+=pa[0]*pb[0];
            }
        }
    }
}

void matr_sub_matr(double* a, double* b, double* c, int w, int h) {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            c[y * w + x] = a[y * w + x] - b[y * w + x];
        }
    }
}


int JordanSolvingSystem(int n, int m, string name, int p , int proc_num , MPI_Comm G){
    int loc_err = 0, glob_err = 0, len = 0;
    const int blockCount = n/m;
    const int blockSize = m;
    const int size = n;
    const int l = n % m;
    const int k = n /m;

    double *a = nullptr, *b = nullptr, *x = nullptr;
    double* C1 = new double[3 * blockSize * blockSize];
    double* C2 = C1 + m*m;
    double* C3 = C1 + 2*m*m;
    double* numaColumn = nullptr;
    double aNorma = 0;
    try{
        numaColumn = new double[m * size];
        len = k/p * n*m + (proc_num<(k % p)? (n * m ) : 0) + (proc_num == k % p ? l *n :  0);
        a = new double[len];
        x = new double[n];
        if(proc_num == k % p || proc_num == 0){
            b = new double[n];
        }
    }catch(...){
        loc_err = 1;
        delete[] numaColumn;
        delete[] a;
        delete[] x;
        delete[] b;
        delete[] C1;
    }
//    printf("Proc: %d , len %d, col %d, and ost: %d\n",proc_num, len, len/(n*m), len%(n*m));
//    return -1;
    MPI_Allreduce(&loc_err,&glob_err, 1, MPI_INT, MPI_SUM, G);

    if(glob_err){
            if(!proc_num) printf("Memory allocation error\n");
            return -1;
    }

    if(name!=""){
        loc_err = init_matrix_file(a,n,m,name,proc_num,p,b,aNorma , numaColumn,G);
        MPI_Allreduce(&loc_err,&glob_err, 1, MPI_INT, MPI_SUM, G);
        if(!proc_num && glob_err < 0){
            printf("read_matrix: can\'t read the file %s (invalid format)\n", name.c_str());
        }
        if(glob_err < 0){
            delete[] numaColumn;
            delete[] a;
            delete[] x;
            delete[] b;
            delete[] C1;
            return -1;
        }
        // printf("%d: %lf %lf %lf\n", proc_num , a[0] ,a[1] ,a[2]);
    }else{
        formula_matr(a,n,m,proc_num,p,b, aNorma,G);

    }
    if(proc_num==0){
        LN;LN;
    }

//    printf("Proc_num: %d\n",proc_num);
//    print_vector(a,len);

//    if(proc_num==k%p){
//        print_vector(b,n);
//    }

    int ErrorInv = 0;
    double_int di,di_out;
    int *columns = new int[k];

    for(int i = 0;i<k;i++){
        columns[i]=-1;
    }
    matr_to_E(C1,m,m);
    matr_to_E(C2,m,m);
    matr_to_E(C3,m,m);
    fill(numaColumn, numaColumn + m*n, 0);
//    cout<<aNorma<<" "<<endl;

//    if(proc_num==k%p){
//        LOG("DO CIKLA");
//        print_vector(b,n);
//    }
    MPI_Barrier(G);
    double time_proc = MPI_Wtime();
    for (int row = 0; row < blockCount; row++) {
    //    printf("Proc_num: %d\n",proc_num);
        di.minNorma = 0;
        di.minNormaCol = blockCount + 1;
        di_out.minNorma = 0;
        di_out.minNormaCol = blockCount + 1;
        // ###############################################################
        for (int col = proc_num; col < blockCount; col+=p) {
            if(columns[col]!=-1) continue;

            get_block(proc_num, p, a, size, blockSize, int((col - proc_num)/p), row, C2, blockSize*blockSize);
      //      print_matrix(C2 , blockSize);
            ErrorInv = JordanInv(C2, blockSize, C1, aNorma);


            if (ErrorInv == 0) {
                double n = norma_matr(C2, blockSize,blockSize);

                if (1/n >  di.minNorma) {
                    di.minNormaCol = col;
                    di.minNorma = 1/n;
                }
            }
        }


        MPI_Allreduce(&di,&di_out,1,MPI_DOUBLE_INT,MPI_MAXLOC,G);

        if (di_out.minNormaCol == blockCount + 1) {
            delete[] C1;
            delete[] numaColumn;
            delete[] a;
            delete[] x;
            delete[] b;
            if(proc_num==0){
                printf("Degenerate Matrix \n");
            }
            return 0;
        }
          columns[di_out.minNormaCol] = row;

        if( di_out.minNormaCol % p == proc_num){
            get_block(proc_num, p, a,size,blockSize,(di_out.minNormaCol - proc_num)/p,0,numaColumn,m * size);
        }

        MPI_Bcast(numaColumn, n*m, MPI_DOUBLE, di_out.minNormaCol % p, G);
//        if(proc_num==k%p){
//            LOG("numaVector");
//            print_vector(numaColumn, m*size);
//        }
        //ОБРАЩАЕМ ГЛАВНЫЙ  ЭЛЕМЕНТ
        get_block(proc_num, p, numaColumn, size, blockSize, 0, row, C2, blockSize* blockSize,true);
        ErrorInv = JordanInv(C2, blockSize, C1, aNorma);

 //       printf("l = %d\n",l);
        //C2 - obratnaya (row,row)
//        if(proc_num==k%p){
//            LOG("C2");
//            print_matrix(C2,m);
//        }
        //#############################################NORMIRUEM STROKU#############################


//    LOG(proc_num);
        for (int col = proc_num; col < blockCount; col+=p) {
            if(columns[col]!=-1) continue;
            get_block(proc_num, p, a, size, blockSize, (col - proc_num)/p, row, C1 , blockSize* blockSize);   // C1 = A[mxm]{0, i}
            mult_matrix(C2, C1, C3, blockSize, blockSize, blockSize);  // C3 = C2 * C1
            push_block(proc_num, p, a, size, blockSize, (col - proc_num)/p, row, C3, blockSize* blockSize);  // C3 = A[mxm]{0, i}
        }
        //cout<<size<<" "<<blockSize<<" "<<size % blockSize;
        if(proc_num == k % p){
            if(l>0){
                //   matr_to_NULL(C1,blockSize);
                get_block(proc_num, p, a, size, blockSize, (blockCount - proc_num)/p, row, C1,m*l);   // C1 = A[mxp]{0, k}
                mult_matrix(C2, C1, C3, m, m, l);                           // C3[mxp] = C2 * C1
                push_block(proc_num, p, a, size, blockSize, (blockCount - proc_num)/p, row, C3, m*l);
            }
            //  matr_to_NULL(C1,blockSize);
            get_b_block(b, blockSize, row, C1, m);
//            LOG("DO MULTA");
//            print_vector(b,n);
            mult_matrix(C2, C1, C3, m, m, 1);
            push_b_block(b, blockSize, row, C3 , m);
//            print_vector(b,n);
        }

   //         printf("Mult row %d\n",row);
      //  get_block(proc_num, p, a,size,blockSize,minNormaCol,0,numaColumn,m * size);
        //######################OBNULAYEM COLONKI##################################################
        for (int i = 0; i < blockCount; i++) {
            if (i == row) {
                continue;
            }
            get_block(proc_num, p, numaColumn, size, blockSize, 0, i, C1, blockSize * blockSize, true);         // C1 = A[mxm]{i, 0}
            for (int j = proc_num; j < blockCount; j+=p) {
                if(columns[j]!=-1) continue;
                get_block(proc_num, p, a, size, blockSize, (j - proc_num)/p, row, C2, blockSize * blockSize);     // C2 = A[mxm]{0, j}
                mult_matrix(C1, C2, C3, blockSize, blockSize, blockSize);  // C3 = C2 * C1         // TODO: СЃР»РµРІР° РёР»Рё СЃРїСЂР°РІР°?
                get_block(proc_num, p, a, size, blockSize, (j - proc_num)/p, i, C2,blockSize * blockSize );       // C2 = A[mxm]{i, j}
                matr_sub_matr(C2, C3, C2, blockSize, blockSize);        // C2 = C2 - C3
                push_block(proc_num, p, a, size, blockSize, (j - proc_num)/p, i, C2 ,blockSize * blockSize );
            }
            if(proc_num== k%p){
                if(l>0){
      //              printf("Ia  tut and l = %d\n",l);
                    //matr_to_NULL(C2, blockSize);           // C1 = A[mxm]{i, 0}
                    get_block(proc_num, p, a, size, blockSize,(k - proc_num)/p, row, C2, m*l);
                    mult_matrix(C1, C2, C3, m ,m ,l);

                    get_block(proc_num, p, a, size, blockSize,(k - proc_num)/p, i, C2, m*l);     // C1 = A[mxl]{i, k}
                    matr_sub_matr(C2, C3, C2, m, m);
       //             print_vector(C2,l*m);// C2 = C1 - C3
                    push_block(proc_num, p, a, size, blockSize, (k - proc_num)/p, i, C2 ,m*l);
                }

                // matr_to_NULL(C2, blockSize);
                get_b_block(b, blockSize,row, C2 ,m);
                mult_matrix(C1, C2, C3, m, m, 1);
                get_b_block(b, blockSize, i, C2, m);
                matr_sub_matr(C2, C3, C2, m, m);
                push_b_block(b, blockSize, i, C2 , m);
         //       print_vector(b,n);
            }
        }
        if(l>0){
            //  matr_to_NULL(C1, blockSize);
            get_block(proc_num, p, numaColumn, size, blockSize, 0, k, C1, m*l,true); //
     //       print_vector(C1,l*m);
            for (int i = proc_num; i < blockCount; i+=p) {
                if(columns[i]!=-1) continue;
                get_block(proc_num, p, a, size, m, (i - proc_num)/p, row, C2, m*m);
                mult_matrix(C1, C2, C3,l,m,m);
                get_block(proc_num, p, a, size, blockSize, (i - proc_num)/p, k, C2, m*l);
                matr_sub_matr(C2, C3, C2, m, m);
                push_block(proc_num, p, a, size, blockSize, (i - proc_num)/p, k, C2, m*l);
            }

            if(proc_num == k % p){
                // matr_to_NULL(C2, blockSize);
                get_block(proc_num, p, a, size, blockSize, (k - proc_num)/p, row, C2, m*l);
                //      mult_matrix(C2, C1, C3, m);
                mult_matrix(C1, C2, C3,l, m, l);                           // C1 = C2 * C3
                get_block(proc_num, p, a, size, blockSize, (k - proc_num)/p, k, C2, l*l);             // C2 = A[mxm]{k, k}
                matr_sub_matr(C2, C3, C2, m, m);
      //          print_vector(C2,l*l);// C2 = C2 - C1
                push_block(proc_num, p, a, size, blockSize, (k - proc_num)/p, k, C2, l*l);
                /////////////////////////////////////////////////////

                //    matr_to_NULL(C2, blockSize);
                get_b_block(b, blockSize, row, C2,m);
                mult_matrix(C1, C2, C3, l, m , 1 );
                get_b_block(b, blockSize, k, C2 , l);
                matr_sub_matr(C2, C3, C2, m, m);
                push_b_block(b, blockSize, k, C2, l);
                //////////////////////////////////////////////////
            }
        }
        //СТАВИМ ТОЧКУ СИНХРОНИЗАЦИИ
        //        LOG(__FILE__);
        //reduce_max(p);
        //        LOG(__DATE__);
       // MPI_Barrier(G);
  //      printf("End of step %d\n",row);

    }
    if(proc_num == k%p){
        if(l>0){
            //matr_to_E(C2, m);
            //matr_to_E(C1, l);

            get_block(proc_num, p, a, size, blockSize, (k - proc_num)/p, k, C2, l*l);
            ErrorInv = JordanInv(C2, l, C1, aNorma);
        //    LOG(ErrorInv);
            if (ErrorInv < 0) {
                loc_err = ERR_DEG_MATRIX;
            }
            //  matr_to_NULL(C1, m);
            if(loc_err ==0 ){
                get_b_block(b, blockSize, k, C1, l);
                mult_matrix(C2,C1,C3, l, l, 1);
                push_b_block(b, blockSize, k, C3, l);//C3 - b(k+1)
                for (int i = 0; i < blockCount; i++) {
                    get_block(proc_num, p, a, size, blockSize, (k - proc_num)/p, i, C1, l*m);
                    mult_matrix(C1, C3, C2, m , l ,1);
                    get_b_block(b, blockSize, i, C1 , m);
                    matr_sub_matr(C1, C2, C1, m, m);
                    push_b_block(b, m, i, C1, m);
                }
            }
        }

        if(loc_err == 0){
            for (int i = 0; i < blockCount; i++) {
                for (int j = 0; j < blockSize; j++) {
                    x[i* blockSize + j] = b[columns[i] * blockSize + j];
                }
            }
            if(l>0){
                for (int i = 0; i < l; i++) {
                    x[blockSize * blockCount + i] = b[blockSize * blockCount + i];
                }
            }
        }
    }

    MPI_Allreduce(&loc_err,&glob_err, 1, MPI_INT, MPI_SUM, G);
    time_proc = MPI_Wtime() - time_proc;

    if(proc_num == 0){

        printf("Time : %lf\n",time_proc);
    }
    if(glob_err){
            if(!proc_num) printf("Matrix is Degenerate\n");
            delete[] columns;
            delete[] C1;
            delete[] numaColumn;
            delete[] x;
            delete[] a;
            delete[] b;
            return -1;
    }

    MPI_Status s;
    if(proc_num != 0 && proc_num == k % p){
        MPI_Send(x ,n , MPI_DOUBLE, 0, 0,G);
    }else if( proc_num == 0 && k % p != 0){
        MPI_Recv(x ,n , MPI_DOUBLE, k%p, 0,G, &s);
    }

    if(proc_num == 0){
        printf("Answer: ");
        print_vector(x,n);
    }

    //RESIDUAL BLOCK
    if(proc_num==0){
        double *c = new double[n];
        fill(c,c + n, 0);
        fill(b,b + n, 0);
        double buf = 0;
        if(name!=""){
            ifstream file;
            file.open(name.c_str());

            for(int j=0 ; j<n; j++){ //y
                for (int i = 0; i < n; i++) { //x
                    file>>buf;
                    c[j] += x[i] * buf;
                    if(i % 2 == 0){
                        b[j]+=buf;
                    }
                }
                c[j]-=b[j];
            }
        }else
        {
            for(int j=0 ; j<n; j++){ //y
                for (int i = 0; i < n; i++) { //x
                    buf = formula_matr_local(i,j);
                    c[j] += x[i] * buf;
                    if(i % 2 == 0){
                        b[j]+=buf;
                    }
                }
                c[j]-=b[j];
            }
        }
        double resid = norma_vec(c,n);
        delete [] c;
        cout<<"n "<<n<<" m "<<m<<" p "<<p<<" Residual: "<<resid<<endl;
   }





    delete[] columns;
    delete[] C1;
    delete[] numaColumn;
    delete[] x;
    delete[] a;
    delete[] b;

    return 0;
}




double get_time(){
    struct rusage buf;
    getrusage(RUSAGE_THREAD,&buf);
    return (double)buf.ru_utime.tv_sec+(double)buf.ru_utime.tv_usec/1000000.;
}

double get_full_time(){
    struct timeval buf;
    gettimeofday(&buf,0);
    return (double)buf.tv_sec+(double)buf.tv_usec/1000000.;
}

void RightSide(double* a, double* b, int n, int m) {
    int k = n / m;
    for (int i = 0; i < n; i++) {
        b[i]=0;
        for (int j = 0; j < n; j += 2) {
            int k1 = j / m ; // row
            int l = (k1 == k ? n%m :m);
            b[i] += a[k1*n*m + i*l + j%m];
        }
    }
}

double norma_vec(double *b,int n){
    double norm = 0;
    for(int i = 0;i < n;i++){
        if(fabs(b[i]) > norm){
            norm = fabs(b[i]);
        }
    }
    return norm;
}

double Residual(double *a, double *b, double *x,int n, int m){
    double *c = new double[n];
    int k = n / m;
    for (int i = 0; i < n; i++) {
        c[i]=0;
        for (int j = 0; j < n; j ++) {
            int k1 = j / m ; // row
            int l = (k1 == k ? n%m :m);
                c[i] += a[k1*n*m + i*l + j%m]*x[j];
        }
    }
    for(int i = 0;i < n;i++){
        c[i]-=b[i];
    }

    double norm = norma_vec(c,n);
    delete []c;
    return norm;
}

double Error( double *x,int n){
    for(int i = 0;i<n;i+=2 ){
        x[i]-=1;
    }
    return norma_vec(x,n);
}
