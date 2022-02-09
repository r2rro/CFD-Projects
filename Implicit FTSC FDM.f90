! Alternative implicit FTSC FDm Solver using TDmA algorithm for parabolic 2D heat transfer equation
!-------------------------------------------------------------
!------------------By Arthur Rostami--------------------------
!-------------------------------------------------------------

program First_CFD_Project_implicit
implicit none

!_____________________________Variables_____________________________

    integer:: i,j,k,xF,yF,xG,yG,z,time
    real:: L,He,e1,alpha,dx,dt,w,t_st,t_fin,e2,m,n
    real,allocatable:: T(:,:,:),X(:),Y(:),G(:,:),H(:,:),P(:,:),Q(:,:)
!_____________________________Define Parameters_____________________

    print*, 'This Method is always stable, you can choose any positive value for dx and dt!'
    print*, 'Please enter delta x:'
    read*, dx
    print*, 'Please enter delta t:'
    read*, dt
    L=2.0                 ! Length of the plate
    He=1.5                 ! Height of the plate
    m=(L/dx)+1            ! number of grid points in x direction
    n=(He/dx)+1            ! number of grid points in x direction
    alpha=23.1E-6         ! thermal diffusivity, m2/s
    w=(alpha*dt)/(dx**2)

 !_____________________________Defining F & G locations_____________________
xF=(0.5/dx)+1 !point F location
yF=(0.5/dx)+1
xG=(1/dx)+1 !point G location
yG=(0.75/dx)+1

!_____________________________Mesh Generation_____________________________

allocate(T(m,n,100000),X(m),Y(n),G(m,n),H(m,n),P(m,n),Q(m,n))
do i=1,m
    X(i)=(i-1)*dx
end do
do j=1,n
    Y(j)=(j-1)*dx
end do

!_____________________________Initial Condition____________________________

do i=1,m
    do j=1,n
        T(i,j,1)=25
    end do
end do

!_____________________________________________________________________________
!___________________________FTCS Implicit Solver______________________________
!_____________________________________________________________________________

k=1           !k is a time step counter

call cpu_time(t_st)
do
    do i=2,m
        do j=1,n
            T(i,j,k+1)=60
        end do
    end do

z=1
!_______________________TDMA______________________
do
    do j=2,n-1
        P(2,j)=w/(1+4*w);H(2,j)=w*(T(2,j+1,k+1)+T(2,j-1,k+1))+T(2,j,k)+T(1,j,k+1)*w ;Q(2,j)=H(2,j)/(1+4*w)
            do i=3,m-1
            H(i,j)=w*T(i,j-1,k+1)+w*T(i,j+1,k+1)+T(i,j,k)
            P(i,j)=w/(1+4*w-w*P(i-1,j))
            Q(i,j)=(w*Q(i-1,j)+H(i,j))/(1+4*w-w*P(i-1,j))
            end do
        T(1,j,k+1)=(T(2,j,k+1)+10*dx)/(0.5*dx+1)!BC on wall AB
        do i=m-1,2,-1
            G(i,j)=T(i,j,k+1)
            T(i,j,k+1)=P(i,j)*T(i+1,j,k+1)+Q(i,j)
        end do
    end do
    do i=1,m
        T(i,1,k+1)=T(i,2,k+1) !BC on wall BD
    end do
    do i=1,m
        T(i,n,k+1)=100  !BC on wall AC
    end do
    do j=2,n-1
        T(m,j,k+1)=25 !BC on wall CD
    end do

e2=0              !error calculation
do i=2,m-1
    do j=2,n-1
        e2=ABS(T(i,j,k+1)-G(i,j))+e2
    end do
end do
if(e2<0.0001) exit
z=z+1
end do
e1=0
do i=2,m-1
    do j=2,n-1
        e1=ABS(T(i,j,k+1)-T(i,j,k))+e1
    end do
end do

if(e1<0.0001) exit

k=k+1
end do
call cpu_time(t_fin)
 !_____________________________Reporting solution information______________________

    print*, 'w=                                                 ',w
    print*, 'Number of iterations to reach steady state:    ',k
    print*, 'Maximum difference between the last two time steps:',e1, 'Celsius'
    print*, 'Solution will be Steady at:                        ',k*dt,'seconds'
    print*, 'Calculation time is:                               ',t_fin-t_st, 'seconds'

!_____________________________Saving the results_____________________________

   open(10,file='Temp contour at SS.plt')
    write(10,*)'VARiABLES = "X", "Y", "T"'
    do j=1,n
        do i=1,m
            write(10,*) x(i),y(j),T(i,j,k)
        end do
    end do

    open(1000,file='Temp contour at 5000.plt')
    write(1000,*)'VARiABLES = "X", "Y", "T"'
    do j=1,n
        do i=1,m
            write(1000,*) x(i),y(j),T(i,j,5000)
        end do
    end do

    open(10000,file='Temp contour at 2000.plt')
    write(10000,*)'VARiABLES = "X", "Y", "T"'
    do j=1,n
        do i=1,m
            write(10000,*) x(i),y(j),T(i,j,2000)
        end do
    end do

    open(100,file='Temp contour at 10000.plt')
    write(100,*)'VARiABLES = "X", "Y", "T"'
    do j=1,n
        do i=1,m
            write(100,*) x(i),y(j),T(i,j,10000)
        end do
    end do

    open(20,file='midline temp at 100000.dat')
    write(20,*)"VARiABLES=X,TEmPERATURE"
    do i=1,m
        write(20,*) (i-1)*dx,T(i,yG,100000)
    end do

    open(30,file='midline initial temp.dat')
    write(30,*)"VARiABLES=X,TEmPERATURE"
    do i=1,m
        write(30,*) (i-1)*dx,T(i,yG,1)
    end do

    open(30,file='midline temp at SS.dat')
    write(30,*)"VARiABLES=X,TEmPERATURE"
    do i=1,m
        write(30,*) (i-1)*dx,T(i,yG,k)
    end do

    open(40,file='midline temp at t=2000dt.dat')
    write(40,*)"VARiABLES=X,TEmPERATURE"
    do i=1,m
        write(40,*) (i-1)*dx,T(i,yG,2000)
    end do

    open(50,file='midline temp at t=10000dt.dat')
    write(50,*)"VARiABLES=X,TEmPERATURE"
    do i=1,m
        write(50,*) (i-1)*dx,T(i,yG,10000)
    end do

    open(60,file='midline temp at t=50000dt.dat')
    write(60,*)"VARiABLES=X,TEmPERATURE"
    do i=1,m
        write(60,*) (i-1)*dx,T(i,yG,50000)
    end do

    open(70,file='midline temp at t=100000dt.dat')
    write(70,*)"VARiABLES=X,TEmPERATURE"
    do i=1,m
        write(70,*) (i-1)*dx,T(i,yG,100000)
    end do

    open(80,file='Temp at F.dat')
    write(80,*)"VARiABLES=t,TEmPERATURE"
    do time=1,k
        write(80,*) (time-1)*dt,T(xF,yF,time)
    end do

    open(90,file='Temp at G.dat')
    write(90,*)"VARiABLES=t,TEmPERATURE"
    do time=1,k
        write(90,*) (time-1)*dt,T(xG,yG,time)
    end do

end program First_CFD_Project_implicit
