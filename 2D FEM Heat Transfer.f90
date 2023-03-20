! Implicit FVM solver for parabolic 2D heat transfer equation

program CFD_SECOND_PROJECT_IMPLICIT
    implicit none

!______________________variables________________________

integer:: i,j,K1,K2,He,ymid,L,h_fluid,xmid,rhoC,step,m1,m2,n1,n2,ii,jj
real:: e1,dx,dt,t_st,t_fin,m,n,ap0,alpha,beta,a,b,c,r
real,allocatable:: T_half(:,:),T_new(:,:),T_old(:,:),T3(:,:),T1(:,:),T2(:,:),K(:,:),er(:,:),x(:),y(:),H(:,:),G(:,:)

!___________________Define Parameters____________________

print*, 'Please Enter delta x:'
read*, dx
print*, 'Please Enter delta t:'
read*, dt

L=2                   ! LENGTH OF PLATE
He=1                  ! HEIGHT OF PLATE
m=(2.0/dx)+2          ! number of grid points in x direction (two half steps)
m1=(0.8/dx)+2
m2=(1.6/dx)+1
n=(1.0/dx)+2          ! number of grid points in y direction (two half steps)
n1=(0.4/dx)+2
n2=(1.0/dx)+2
h_fluid=5
rhoC=2000             ! density * c
K1=10                 ! first material
ap0=rhoC*dx*dx/dt

!_________________________mid line positions______________________________

ymid=(0.5/dx)+2
xmid=(1/dx)+2

!_____________________________Mesh generation_____________________________

allocate(T_old(m,n),x(m),y(n),er(m,n),T_new(m,n),T1(m,n),T2(m,n),K(m,n),T_half(m,n),H(m,n),G(m,n),T3(m,n))
x(1)=0
y(1)=0
x(m)=2
y(n)=1
do i=2,m-1
    x(i)=(i-1)*dx-dx/2
end do
do j=2,n-1
    y(j)=(j-1)*dx-dx/2
end do

!_____________________________Initial Condition____________________________

do i=1,m
    do j=1,n
        T_old(i,j)=50
        T_new(i,j)=50
        T_half(i,j)=50
        er(i,j)=0
        K(i,j)=K1
        H(i,j)=0
        G(i,j)=0
    end do
end do

!__________Boundary condition T old______________

do j=2,n
   T_old(m,j)=(h_fluid*100+(2*K1*T_old(m-1,j))/dx)/(2*K1/dx+5)       !BC Q"=H(100-T) & k=10
end do
do i=2,m-1
   T_old(i,n)=T_old(i,n-1)                                           !BC Q"=0
end do
do i=2,m-1
   T_old(i,1)=50*dx/(2*K1)+T_old(i,2)                                !BC Q"=50
end do
do j=1,n
   T_old(1,j)=100*y(j)                                               !BC T=100y
end do

    !________________________________________________________________________
    !_____________________________Implicit Solver____________________________
    !________________________________________________________________________

step=1     !step is time step counter
call cpu_time(t_st)
do while (step<100000)

!__________________________x-sweep________________________________

    T_half(:,:)=T_old(:,:)

    do i=m1,m2
        do j=n1,n2
            K(i,j)=3*((T_old(i,j)/5)+2)
        end do
    end do

!____________________________TDMA_________________________________

    do j=2,n-1
        do i=2,m-1
            if(i==2)then
                a=0
                b=ap0+K(3,j)+K(1,j)+3*dx*dx
                c=-1*K(3,j)
                r=50*dx*dx+K(2,j-1)*T_old(2,j-1)+K(2,j+1)*T_old(2,j+1)+(ap0-K(2,j-1)-K(2,j+1))*T_old(2,j)+K(i-1,j)*100*y(j)
                H(2,j)=c/b
                G(2,j)=r/b
            else if(i==m-1)then
                a=-1*K(m-2,j)
                b=ap0+K(m,j)+K(m-2,j)+3*dx*dx-K(m,j)*((2*K1/dx)/(2*K1/dx+h_fluid))
                c=0
                alpha=K(m,j)*h_fluid*100/(2*K1/dx+h_fluid)
                r=50*dx*dx+K(m-1,j-1)*T_old(m-1,j-1)+K(m-1,j+1)*T_old(m-1,j+1)+(ap0-K(m-1,j-1)-K(m-1,j+1))*T_old(m-1,j)+alpha
                H(m-1,j)=0
                G(m-1,j)=(r-a*G(m-2,j))/(b-a*H(m-2,j))
            else
                a=-1*K(i-1,j)
                b=ap0+K(i+1,j)+K(i-1,j)+3*dx*dx
                c=-1*K(i+1,j)
                r=50*dx*dx+K(i,j-1)*T_old(i,j-1)+K(i,j+1)*T_old(i,j+1)+(ap0-K(i,j-1)-K(i,j+1))*T_old(i,j)
                H(i,j)=c/(b-a*H(i-1,j))
                G(i,j)=(r-a*G(i-1,j))/(b-a*H(i-1,j))
            end if
    end do
    do ii=m-1,2,-1
        T_half(ii,j)=G(ii,j)-H(ii,j)*T_half(ii+1,j)
    end do
    end do

!_____________Boundary condition update: T half______________

    do j=2,n-1
        T_half(m,j)=(h_fluid*100+2*K1/dx*T_half(m-1,j))/(2*K1/dx+h_fluid)       !BC Q"=H(100-T) & k=10
    end do
    do i=2,m
        T_half(i,1)=50*dx/(2*K1)+T_half(i,2)                                    !BC Q"=50
    end do
    do i=2,m
        T_half(i,n)=T_half(i,n-1)                                               !BC Q"=0
    end do

!__________________________y-sweep________________________________

    T_new(:,:)=T_half(:,:)

    do i=m1,m2
        do j=n1,n2
            K(i,j)=3*((T_half(i,j)/5)+2)
        end do
    end do

!____________________________TDMA_________________________________

do i=2,m-1
    do j=2,n-1
      if(j==2)then
        a=0
        b=ap0+K(i,1)+K(i,3)+3*dx*dx-K(i,1)
        c=-1*K(i,3)
        beta=K(i,j-1)*(50*dx/(2*K1))
        r=50*dx*dx+K(i-1,2)*T_half(i-1,2)+K(i+1,2)*T_half(i+1,2)+(ap0-K(i-1,2)-K(i+1,2))*T_half(i,2)+beta
        H(i,2)=c/b
        G(i,2)=r/b
      else if(j==n-1)then
        a=-1*K(i,n-2)
        b=ap0+K(i,n-2)+3*dx*dx
        c=0
        r=50*dx*dx+K(i-1,n-1)*T_half(i-1,n-1)+K(i+1,n-1)*T_half(i+1,n-1)+(ap0-K(i-1,n-1)-K(i+1,n-1))*T_half(i,n-1)
        H(i,n-1)=0
        G(i,n-1)=(r-a*G(i,n-2))/(b-a*H(i,n-2))
      else
        a=-1*K(i,j-1)
        b=ap0+K(i,j-1)+K(i,j+1)+3*dx*dx
        c=-1*K(i,j+1)
        r=50*dx*dx+K(i+1,j)*T_half(i+1,j)+K(i-1,j)*T_half(i-1,j)+(ap0-K(i+1,j)-K(i-1,j))*T_half(i,j)
        H(i,j)=c/(b-a*H(i,j-1))
        G(i,j)=(r-a*G(i,j-1))/(b-a*H(i,j-1))
        end if
    end do
    do jj=n-1,2,-1
        T_new(i,jj)=G(i,jj)-H(i,jj)*T_new(i,jj+1)
    end do
end do

!_____________Boundary condition update: T half______________

    do j=2,n-1
        T_new(m,j)=(h_fluid*100+2*K1/dx*T_new(m-1,j))/(2*K1/dx+h_fluid)           !BC Q"=H(100-T) & k=10
    end do
    do i=2,m
        T_new(i,1)=50*dx/(2*K1)+T_new(i,2)                                        !BC Q"=50
    end do
    do i=2,m
        T_new(i,n)=T_new(i,n-1)                                                   !BC Q"=0
    end do
!______________________Saving Temp at chosen times ________________________

    if(step==2000) T1(:,:)=T_new(:,:)     !saving temperature at t=2000dt
    if(step==10000) T2(:,:)=T_new(:,:)    !saving temperature at t=10000dt
    if(step==20000) T3(:,:)=T_new(:,:)    !saving temperature at t=20000dt

!________________error calculation_____________________

    er(:,:)=ABS(T_new(:,:)-T_old(:,:))
    e1=maxval(er)
    T_old(:,:)=T_new(:,:)
    if(e1<0.00001) exit
    step=step+1
end do
call cpu_time(t_fin)

!_____________________________Reporting solution information_____________________________

    print*, 'Number of iterations to reach steady state:    ',step
    print*, 'Maximum difference between the last two time steps:',e1, 'Celsius'
    print*, 'Solution will be Steady at:                        ',step*dt,'seconds'
    print*, 'Calculation time is:                               ',t_fin-t_st, 'seconds'

!_____________________________Saving the results_____________________________

    open(10,file='Temp contour at SS.plt')
    write(10,*)'VARIABLES = "X", "Y", "T"'
    write(10,*)'ZONE I=',m,'J=',n
    do j=1,n
        do i=1,m
            write(10,*) x(i),y(j),T_old(i,j)
        end do
    end do


    open(20,file='Temp contour at 2000.plt')
    write(20,*)'VARIABLES = "X", "Y", "T"'
    write(20,*)'ZONE I=',m,'J=',n
    do j=1,n
        do i=1,m
            write(20,*) x(i),y(j),T1(i,j)
        end do
    end do

    open(30,file='Temp contour at 10000.plt')
    write(30,*)'VARIABLES = "X", "Y", "T"'
    write(30,*)'ZONE I=',m,'J=',n
    do j=1,n
        do i=1,m
            write(30,*) x(i),y(j),T2(i,j)
        end do
    end do

    open(40,file='Temp contour at 30000.plt')
    write(40,*)'VARIABLES = "X", "Y", "T"'
    write(40,*)'ZONE I=',m,'J=',n
    do j=1,n
        do i=1,m
            write(40,*) x(i),y(j),T3(i,j)
        end do
    end do

    open(60,file='x-Midline temp at 2000.dat')
    write(60,*)"VARIABLES=X,TEMPERATURE"
    write(60,*)'ZONE I=',m,'J=',n
    do i=1,m
        write(60,*) (i-1)*dx,T1(i,ymid)
    end do

    open(70,file='y-Midline temp at 2000.dat')
    write(70,*)"VARIABLES=Y,TEMPERATURE"
    do j=1,n
        write(70,*) (j-1)*dx,T1(xmid,j)
    end do

    open(61,file='x-Midline temp at 10000.dat')
    write(61,*)"VARIABLES=X,TEMPERATURE"
    do i=1,m
        write(61,*) (i-1)*dx,T2(i,ymid)
    end do

    open(71,file='y-Midline temp at 10000.dat')
    write(71,*)"VARIABLES=Y,TEMPERATURE"
    do j=1,n
        write(71,*) (j-1)*dx,T2(xmid,j)
    end do

    open(62,file='x-Midline temp at SS.dat')
    write(62,*)"VARIABLES=X,TEMPERATURE"
    do i=1,m
        write(62,*) (i-1)*dx,T_old(i,ymid)
    end do

    open(72,file='y-Midline temp at SS.dat')
    write(72,*)"VARIABLES=Y,TEMPERATURE"
    do j=1,n
        write(72,*) (j-1)*dx,T_old(xmid,j)
    end do

end program CFD_SECOND_PROJECT_IMPLICIT
