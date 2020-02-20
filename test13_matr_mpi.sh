#!/bin/bash

for k in 1 2 ; do echo "k=$k a.txt ---------------\n\n"   ; mpirun -np $k ./a.out 4 3 a.txt   ; done
for k in 1 2 ; do echo "k=$k a20.txt ---------------\n\n" ; mpirun -np $k ./a.out 4 3 a20.txt ; done
for k in 1 2 ; do echo "k=$k b.txt ---------------\n\n"   ; mpirun -np $k ./a.out 4 3 b.txt   ; done
for k in 1 2 ; do echo "k=$k c.txt ---------------\n\n"   ; mpirun -np $k ./a.out 6 3 c.txt   ; done
for k in 1 2 ; do echo "k=$k d.txt ---------------\n\n"   ; mpirun -np $k ./a.out 6 3 d.txt   ; done
for k in 1 2 ; do echo "k=$k e.txt ---------------\n\n"   ; mpirun -np $k ./a.out 6 3 e.txt   ; done
for k in 1 2 ; do echo "k=$k f.txt ---------------\n\n"   ; mpirun -np $k ./a.out 6 3 f.txt   ; done
