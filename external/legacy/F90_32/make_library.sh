# How to create call an external routine from IDL with a C wrapper

rm -rf bvls.o bvls_wrapper.o

#1. Create an object out of the routine you want to use
#   IMPORTANT= This example works only with the Intel Fortran compilers
   
   ifort -fPIC -c bvls.f90
   
#2. Create an object out of the C wrapper
#   IMPORTANT= This example works only with the Intel Fortran compilers
   
   icc -fPIC -c bvls_wrapper.c
   
  
#3. Create a shared library

   ld -fPIC -shared -o bvls.so bvls.o bvls_wrapper.o  


#Now you can use this library in your call to CALL_EXTERNAL in IDL
chmod 644 bvls.so
