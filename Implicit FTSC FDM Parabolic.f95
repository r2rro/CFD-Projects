        ! Implicit FTSC FDM Solver using TDMA algorithm and ADI method for parabolic 2D heat transfer equation
        !-------------------------------------------------------------
        !------------------By Arthur Rostami -------------------------
        !-------------------------------------------------------------

program CFD_FIRST_PROJECT_IMPLICIT
    implicit none
    !_____________________________Variables_____________________________
    integer:: i,j,k,xF,yF,xG,yG,L,H,time
    real:: e,alpha,dx,dt,w,t_st,t_fin,m,n
    real,allocatable:: T(:,:),x(:),y(:),er(:,:),TF(:,:),TG(:,:),T0(:,:),T1(:,:),T2(:,:),T3(:,:),T4(:,:)
    real,allocatable:: a(:,:),b(:,:),c(:,:),d(:,:),c_new(:,:),d_new(:,:),temp1(:,:),temp2(:,:)

    !_____________________________Define Parameters_____________________
    print*, 'This Method is always stable, you can choose any positive value for dx and dt!'
    print*, 'Please enter delta x:'
    read*, dx
    print*, 'Please enter delta t:'
    read*, dt
    L=2.0                   ! Length of the plate
    H=1.5                 ! Height of the plate
    m=(L/dx)+1            ! number of grid points in x direction
    n=(H/dx)+1            ! number of grid points in x direction
    alpha=23.1E-6         ! thermal diffusivity, m2/s
    w=(alpha*dt)/(dx**2)

    allocate(T(m,n),x(m),y(n),er(m,n),TF(1,100000),TG(1,100000),T0(m,n),T1(m,n),T2(m,n),T3(m,n),T4(m,n),temp1(m,n),temp2(m,n))
    allocate(a(2:m-1,2:n-1),b(2:m-1,2:n-1),c(2:m-1,2:n-1),d(2:m-1,2:n-1),c_new(2:m-1,2:n-1),d_new(2:m-1,2:n-1))
    !_____________________________Defining F & G locations_____________________
    xF=(0.5/dx)+1 !point F location
    yF=(0.5/dx)+1

    xG=(1/dx)+1   !point G location
    yG=(0.75/dx)+1

    !_____________________________Mesh Generation_____________________________
    do i=1,m
        x(i)=(i-1)*dx
    end do
    do j=1,n
        y(j)=(j-1)*dx
    end do

    !_____________________________Initial Condition____________________________
    do i=1,m
        do j=1,n
            T(i,j)=25
            temp1(i,j)=25
            temp2(i,j)=25
        end do
    end do
    !_____________________________________________________________________________
    !____________________FTCS Implicit Solver Using ADI method____________________
    !_____________________________________________________________________________
    k=1 !time step counter
    call cpu_time(t_st)
do
    !_____________________________Boundary condition_____________________________
    do j=2,n-1
        T(1,j)=(T(2,j)+10*dx)/(0.5*dx+1)!BC on wall AB
        temp1(i,j)=(temp1(2,j)+10*dx)/(0.5*dx+1)
        temp2(i,j)=(temp2(2,j)+10*dx)/(0.5*dx+1)
    end do

    do i=1,m
        T(i,1)=T(i,2) !BC on wall BD
        temp1(i,1)=temp1(i,2)
        temp2(i,1)=temp1(i,2)
    end do

    do i=1,m
        T(i,n)=100 !BC on wall AC
        temp1(i,n)=100
        temp2(i,n)=100
    end do

    do j=2,n-1
        T(m,j)=25 !BC on wall CD
        temp1(m,j)=25
        temp2(m,j)=25
    end do

    !_____________________________Saving Temp at chosen times and locations F and G_____________________________
    if(k==1)     T0(:,:)=T(:,:)       !saving temperature at t=0*dt
    if(k==2000)  T1(:,:)=T(:,:)       !saving temperature at t=2000*dt
    if(k==10000) T2(:,:)=T(:,:)       !saving temperature at t=10000*dt
    if(k==50000) T3(:,:)=T(:,:)       !saving temperature at t=50000*dt
    if(k==100000)T4(:,:)=T(:,:)       !saving temperature at t=100000*dt
    TF(1,k)=T(xF,yF)                  !saving temperature of point F for different times
    TG(1,k)=T(xG,yG)                  !saving temperature of point G for different times

    !_____ADI for first half-step_______
    !___________________________________
            do j=2,n-1
                do i=3,m-2
                    a(i,j)=w
                    b(i,j)=-2-2*w
                    c(i,j)=w
                    d(i,j)=-w*T(i,j+1)-w*T(i,j-1)+(-2+2*w)*T(i,j)
                end do

                a(2,j)=0
                b(2,j)=-2-2*w+1/(1+0.5*dx)
                c(2,j)=w
                d(2,j)=-w*T(i,j+1)-w*T(i,j-1)+(-2+2*w)*T(i,j)-10*dx/(1+0.5*dx)

                a(m-1,j)=w
                b(m-1,j)=-2-2*w
                c(m-1,j)=0
                d(m-1,j)=-w*T(i,j+1)-w*T(i,j-1)+(-2+2*w)*T(i,j)-w*T(m,j)

            end do
                !_____TDMA Solver_____
                do j=2,n-1
                    c_new(2,j)=c(2,j)/b(2,j)
                    d_new(2,j)=d(2,j)/b(2,j)
                    do i=3,m-1
                        c_new(i,j)=c(i,j)/(b(i,j)-a(i,j)*c_new(i-1,j))
                        d_new(i,j)=(d(i,j)-a(i,j)*d_new(i-1,j))/(b(i,j)-a(i,j)*c_new(i-1,j))
                    end do
                    temp1(m-1,j)=d_new(m-1,j)
                    do i=m-2,2,-1
                        temp1(i,j)=d_new(i,j)-c_new(i,j)*temp1(i+1,j) ! Saving half-step(k+0.5) values in a temporary variable
                    end do
                end do
    !_____ADI for second half-step_______
    !____________________________________
            do i=2,m-1
                do j=3,n-2
                    a(i,j)=w
                    b(i,j)=-2-2*w
                    c(i,j)=w
                    d(i,j)=-w*temp1(i+1,j)-w*temp1(i-1,j)+(-2+2*w)*temp1(i,j)
                end do

                a(i,2)=0
                b(i,2)=-2-w
                c(i,2)=w
                d(i,2)=-w*temp1(i+1,j)-w*temp1(i-1,j)+(-2+2*w)*temp1(i,j)

                a(i,n-1)=w
                b(i,n-1)=-2-2*w
                c(i,n-1)=0
                d(i,n-1)=-w*temp1(i+1,j)-w*temp1(i-1,j)+(-2+2*w)*temp1(i,j)-w*T(i,n)

            end do

                !_____TDMA Solver_____
                do i=2,m-1
                    c_new(i,2)=c(i,2)/b(i,2)
                    d_new(i,2)=d(i,2)/b(i,2)
                    do j=3,n-1
                        c_new(i,j)=c(i,j)/(b(i,j)-a(i,j)*c_new(i,j-1))
                        d_new(i,j)=(d(i,j)-a(i,j)*d_new(i,j-1))/(b(i,j)-a(i,j)*c_new(i,j-1))
                    end do
                    temp2(i,n-1)=d_new(i,n-1)
                    do j=n-2,2,-1
                        temp2(i,j)=d_new(i,j)-c_new(i,j)*temp2(i,j+1) ! Saving second half-step(k+1)
                    end do
                end do

    er(:,:)=ABS(temp2(:,:)-T(:,:)) ! This variable is defined to determine the time the solution will reach steady state
    e=maxval(er) !Maximum difference
    T(:,:)=temp2(:,:)

    if(e<0.0001) exit
    k=k+1
end do
    call cpu_time(t_fin)
    !_____________________________Reporting solution information_____________________________
    print*, 'w=                                                 ',w
    print*, 'Number of iterations to reach steady state:    ',k
    print*, 'Maximum difference between the last two time steps:',e, 'Celsius'
    print*, 'Solution will be Steady at:                        ',k*dt,'seconds'
    print*, 'Calculation time is:                               ',t_fin-t_st, 'seconds'

    !_____________________________Saving the results_____________________________

    open(10,file='Temp contour.plt')
    write(10,*)'VARIABLES = "X", "Y", "T"'
    do j=1,n
        do i=1,m
            write(10,*) x(i),y(j),T(i,j)
        end do
    end do

    open(20,file='Midline temp at SS.dat')
    write(20,*)"VARIABLES=X,TEMPERATURE"
    do i=1,m
        write(20,*) (i-1)*dx,T(i,yG)
    end do

    open(30,file='Midline initial temp.dat')
    write(30,*)"VARIABLES=X,TEMPERATURE"
    do i=1,m
        write(30,*) (i-1)*dx,T0(i,yG)
    end do

    open(40,file='Midline temp at t=2000dt.dat')
    write(40,*)"VARIABLES=X,TEMPERATURE"
    do i=1,m
        write(40,*) (i-1)*dx,T1(i,yG)
    end do

    open(50,file='Midline temp at t=10000dt.dat')
    write(50,*)"VARIABLES=X,TEMPERATURE"
    do i=1,m
        write(50,*) (i-1)*dx,T2(i,yG)
    end do

    open(60,file='Midline temp at t=50000dt.dat')
    write(60,*)"VARIABLES=X,TEMPERATURE"
    do i=1,m
        write(60,*) (i-1)*dx,T3(i,yG)
    end do

    open(70,file='Midline temp at t=100000dt.dat')
    write(70,*)"VARIABLES=X,TEMPERATURE"
    do i=1,m
        write(70,*) (i-1)*dx,T4(i,yG)
    end do

    open(80,file='Temp at F.dat')
    write(80,*)"VARIABLES=t,TEMPERATURE"
    do time=1,k
        write(80,*) (time-1)*dt,TF(1,time)
    end do

    open(90,file='Temp at G.dat')
    write(90,*)"VARIABLES=t,TEMPERATURE"
    do time=1,k
        write(90,*) (i-1)*dt,TG(1,time)
    end do

end program CFD_FIRST_PROJECT_IMPLICIT
