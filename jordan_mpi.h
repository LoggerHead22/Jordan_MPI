#pragma once
#include "mpi.h"
#include <cmath>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cstring>
#include <sys/time.h>
#include <sys/sysinfo.h>
#include <sys/resource.h>
#include <float.h>
using namespace std;

#define ERR_CANNOT_OPEN -1
#define ERR_READ -2
#define ERR_DEG_MATRIX -3   // вырожденная матрица
#define EPS 1e-15

#define LOG(...) std::cout<< #__VA_ARGS__<< " = " <<(__VA_ARGS__)<< "\n"
#define LN       std::cout << "\n";
double Random();
void generate_file(int n);
void print_matrix(double* m, int size);
void new_print_matrix(double* m, int size, int blockSize);
void print_equa(double* a, double* b, int size);
void print_vector(int* m, int size);
void print_vector(double* m, int size);
void print_matrix_rect(double* m, int w, int h);
void matr_to_E(double* m, int size, int blockSize);
void matr_to_NULL(double* m, int size);
void formula_matr(double* a, int n,int m, int proc_num, int p, double *b, double &Norma,MPI_Comm G);
double formula_matr_local(int i, int j);
int JordanInv(double* a, int size, double* b, double norma=1);
int read_matrix(double* a, int n, int m, const string& name, int proc_num, int p, double *b, double &Norma, double *CopyCol, MPI_Comm G);
int init_matrix_file(double* a, int n, int m, const string& name, int proc_num, int p, double *b, double &Norma, double *CopyCol, MPI_Comm G);
double norma_matr(double* a, int n, int m);
void get_block(int proc_num, int p, double* a, int aSize, int blockSize, int x, int y, double* c3, int count, bool isNumaColumn = false);
void push_block(int proc_num, int p, double* a, int aSize, int blockSize, int x, int y, double* c3, int count, bool isNumaColumn = false);
void get_b_block(double* b, int blockSize, int row, double* c3, int count);
void push_b_block(double* b, int blockSize, int row, double* c3, int count);
void mult_matrix(double *a, double *b, double *c, int n, int, int);
void matr_sub_matr(double* a, double* b, double* c, int w, int h);
void Matr_mlt(double* a, double* b, double* c, int _i, int _r, int _j);
// get_block - учесть что size % blockCount не всегда != 0
// учесть что blockCount == 0
void* Jordan_Solving_System(void * arg);
void RightSide(double* a, double* b, int n, int blockSize);
double norma_vec(double *b,int n);
double Residual(double *a, double *b, double *x,int n, int m);
double Error( double *x,int n);
double residual_2(double *a, double *x, double *b, int n);
double residual_1(double *a, double *x, double *b, int n);
double residual_inf(double *a, double *x, double *b, int n);
double get_full_time();
double get_time();


int JordanSolvingSystem(int , int , string, int , int , MPI_Comm);

struct double_int{
        double minNorma = 0;
        int minNormaCol = -1;
};
