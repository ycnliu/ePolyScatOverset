#
echo Compiler and load options -O2 -m64 -Wl,-framework -Wl,Accelerate
/opt/openmpi-64-gfortran/bin/mpif90 -o Testmmp.exe -I. -O2 -m64 -Wl,-framework -Wl,Accelerate Testmmp.f90
/opt/openmpi-64-gfortran/bin/mpirun -np 4 ./Testmmp.exe <<eoi
100   1 32768
200   1  4096
400   1   512
eoi

/opt/openmpi-64-gfortran/bin/mpirun -np 2 ./Testmmp.exe <<eoi
100   1 32768
200   1  4096
400   1   512
eoi


echo Compiler and load options -Ofast -fexternal-blas -m64 -Wl,-framework -Wl,Accelerate
/opt/openmpi-64-gfortran/bin/mpif90 -o Testmmp.exe -I. -Ofast -fexternal-blas -m64 -Wl,-framework -Wl,Accelerate Testmmp.f90
/opt/openmpi-64-gfortran/bin/mpirun -np 4 ./Testmmp.exe <<eoi
100   1 32768
200   1  4096
400   1   512
eoi

/opt/openmpi-64-gfortran/bin/mpirun -np 2 ./Testmmp.exe <<eoi
100   1 32768
200   1  4096
400   1   512
eoi


