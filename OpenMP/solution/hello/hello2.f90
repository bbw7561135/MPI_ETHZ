program hello

integer i
!$ integer OMP_GET_THREAD_NUM
i=-1
!$ i=OMP_GET_THREAD_NUM()
print *,'hello world',i

stop
end program hello
