module ellipse                    
    use phyDomainInit
    implicit none
    real(8),dimension(imax,jmax)::alpha,beta,gamma,jacobi  !the cofficients of ellipse equation
   
    real(8),dimension(imax,jmax)::x_xi,y_xi , &            !the metric of ellipse equation  
                                  x_eta,y_eta , &
                                  x_xieta,y_xieta
    
    real(8),dimension(imax,0:1)::dd,theta                  !the space and the  angle around the boundary              
    
    character(len=1)::judge_source                         !variable used to determine  wheather to possion's equation  
    integer::flag                                          !the return value of iteration used to determine wheather the coordinates converge 
    integer::conv                                          !the return value of Source_conv_judge used to determine wheather the source term converge 
    
contains

subroutine solver         !solver      
    implicit none  
    integer::i,j
    integer::count             !outside iterative variable
    
    write(*,*) "solve"

    call initialize     
    call coeffCal           

    !*********************************
    !wheather to use possion's equation:  y:use possion's equation
                                      ! n:use laplace's equation
    write(*,*) "Wheather to use possion's equation ?"
    write(*,*) "(y¡ªyes,n¡ªno)"
    do while(.true.)
        read(*,*) judge_source
        if(judge_source=='y' .or. judge_source=='n') then    
            exit
        else
            write(*,*)  "Error£¬please input: y or n !"
        end if
    end do 
 
    write(*,*)  "¡­¡­Calculating¡­¡­"

    do count=1,iter_out
        write(*,*)  
        write(*,*)  "Outside loop count:",count
        
        call explicitMethod  
        !print the situation of wheather x and y coordinates converge
        if(flag==1)   then
            write(*,*)  "...Coordinates converge..."
            if(judge_source=='n')  exit
        else
            write(*,*)  "...Coordinates disconverge..."
        end if  
        
        if(judge_source=='y')  then
            call sourceConvJudge 
            if(conv==1) then 
                write(*,*) "...Satisfy the boundary constraints..." 
                exit
            else
                write(*,*)  "...Dissatisfy the boundary constraints..."
                call sourceRenew
            end if  
        end if
    end do
    write(*,*)  "             ...Reach the maximum outside loops..."
    write(*,*)  "             ...Finish..."
 
end subroutine 
  
subroutine coeffCal        !the cofficients of the possion's or the laplace's equation
    implicit none
    integer::i,j
    !*******************************************  
    !the metric of x to xi 
    !excluding the grid points on the left and the riht boundaries 
    !adopt central-difference method , 2 order accuracy
    do j=1,jmax  
    do i=2,imax-1
        x_xi(i,j)=(x(i+1,j)-x(i-1,j))/2
        y_xi(i,j)=(y(i+1,j)-y(i-1,j))/2
    end do
    end do
    !x to xi of the points on the left and the right boundaries
    !adopt the forward-difference and backward-difference  method,1 order accuracy
    do j=1,jmax
        x_xi(1,j)=x(2,j)-x(1,j)
        y_xi(1,j)=y(2,j)-y(1,j)
        x_xi(imax,j)=x(imax,j)-x(imax-1,j)
        y_xi(imax,j)=y(imax,j)-y(imax-1,j) 
    end do
    
    !x to eta
    !excluding the grid points on the bottom and the top boundaries 
    !adopt central-difference method , 2 order accuracy
    do i=1,imax  
    do j=2,jmax-1
        x_eta(i,j)=(x(i,j+1)-x(i,j-1))/2
        y_eta(i,j)=(y(i,j+1)-y(i,j-1))/2
    end do
    end do
    !x to eta of the points on the top and bottom boundaries
    !adopt the forward-difference and backward-difference  method,1 order accuracy
    do i=1,imax
        x_eta(i,1)=x(i,2)-x(i,1)
        y_eta(i,1)=y(i,2)-y(i,1)
        x_eta(i,jmax)=x(i,jmax)-x(i,jmax-1)
        y_eta(i,jmax)=y(i,jmax)-y(i,jmax-1)
    end do
    
    !*********************************************************************
    !x to xieta and y to xieta  of  inner points
    !adopt 4 points method 
    !alternative method: use x_xi to eta or x_eta to xi to calculate
    do j=2,jmax-1  
    do i=2,imax-1 
        x_xieta(i,j)=(x(i-1,j-1)-x(i-1,j+1)-x(i+1,j-1)+x(i+1,j+1))/4.0
        y_xieta(i,j)=(y(i-1,j-1)-y(i-1,j+1)-y(i+1,j-1)+y(i+1,j+1))/4.0
    end do
    end do
    !points on the left and the right boundaries(exclude the points coincide with the top and the bottom boundaries)
    do j=2,jmax-1
        !x_xi to eta, central_difference method
        x_xieta(1,j)=(x_xi(1,j+1)-x_xi(1,j-1))/2
        y_xieta(1,j)=(y_xi(1,j+1)-y_xi(1,j-1))/2
        x_xieta(imax,j)=(x_xi(imax,j+1)-x_xi(imax,j-1))/2
        y_xieta(imax,j)=(y_xi(imax,j+1)-y_xi(imax,j-1))/2
    end do 
    !points on the bottom and the top boundaries(include the points coincide with the left and the fight boundaries)
    do i=1,imax
        x_xieta(i,1)=x_xi(i,2)-x_xi(i,1)
        y_xieta(i,1)=y_xi(i,2)-y_xi(i,1)
        x_xieta(i,jmax)=x_xi(i,jmax)-x_xi(i,jmax-1)
        y_xieta(i,jmax)=y_xi(i,jmax)-y_xi(i,jmax-1)
    end do 
    
    !*******************************************************************
    !the cofficients of all points
    alpha =x_eta**2 + y_eta**2
    beta  =x_xi*x_eta + y_xi*y_eta
    gamma =x_xi**2 + y_xi**2 
    jacobi=x_xi*y_eta- x_eta*y_xi
     
end subroutine

subroutine explicitMethod
    implicit none 
    integer::iter
    integer::i,j
    real(8),dimension(imax,jmax)::x1,y1
    !***************************************************  
   
    do iter=1,itermax
        !inner points
        do j=2,jmax-1
        do i=2,imax-1
            x1(i,j)= 0.5/(alpha(i,j)+gamma(i,j))*(alpha(i,j)*(x(i-1,j)+x(i+1,j))+gamma(i,j)*(x(i,j+1)+x(i,j-1)) &
                   & -2*beta(i,j)*x_xieta(i,j)+jacobi(i,j)**2*(p(i,j)*x_xi(i,j)+q(i,j)*x_eta(i,j)))
            y1(i,j)= 0.5/(alpha(i,j)+gamma(i,j))*(alpha(i,j)*(y(i-1,j)+y(i+1,j))+gamma(i,j)*(y(i,j+1)+y(i,j-1)) &
                   & -2*beta(i,j)*y_xieta(i,j)+jacobi(i,j)**2*(p(i,j)*y_xi(i,j)+q(i,j)*y_eta(i,j)))       
        end do
        end do 
        !points on the left and the right boundaries
        !i=1
        !do j=2,jmax-1
        !    x1(i,j)=x(i,j)
        !    y1(i,j)= 0.5/gamma(i,j)*(alpha(i,j)*(y_xi(i+1,j)-y_xi(i,j))+gamma(i,j)*(y(i,j+1)+y(i,j-1)) &
        !           & -2*beta(i,j)*y_xieta(i,j)+(jacobi(i,j)**2)*(p(i,j)*y_xi(i,j)+q(i,j)*y_eta(i,j)))       
        !end do
        !i=imax
        !do j=2,jmax-1
        !    x1(i,j)=x(i,j)
        !    y1(i,j)= 0.5/gamma(i,j)*(alpha(i,j)*(y_xi(i,j)-y_xi(i-1,j))+gamma(i,j)*(y(i,j+1)+y(i,j-1)) &
        !           & -2*beta(i,j)*y_xieta(i,j)+(jacobi(i,j)**2)*(p(i,j)*y_xi(i,j)+q(i,j)*y_eta(i,j)))       
        !end do 
        !*************************************************************
        !verify  wheather the x and y coordinates  converge
        flag=1
        do j=2,jmax-1
        do i=2,imax-1
            if(abs(x1(i,j)-x(i,j))>eps .or. abs(y1(i,j)-y(i,j))>eps )  then
                flag=0
                exit
            end if
        end do
        end do
       !*********************************************** 
       !renew the coordinates  of the physical domain
       !exclude the points on the inside and the outside boundaries 
       do j=2,jmax-1
       do i=2,imax-1
           x(i,j)=x(i,j)+omg*(x1(i,j)-x(i,j))
           y(i,j)=y(i,j)+omg*(y1(i,j)-y(i,j))
       end do
       end do
       !!adjust the points on the inside boundary except aerofoil part
       !do i=2,(imax-aerofoil_num)/2-1
       !    x(i,1)=2*x(i,2)-x(i,3)
       !    x(imax+1-i,1)=2*x(imax+1-i,2)-x(imax+1-i,3)
       !end do
       !!adjust the points on the outside boundary
       !do i=2,(imax-aerofoil_num)/2
       !    x(i,jmax)=2*x(i,jmax-1)-x(i,jmax-2)
       !    x(imax+1-i,jmax)=2*x(imax+1-i,jmax-1)-x(imax+1-i,jmax-2)
       !end do
       
  
       !*******************************************
       !renew the cofficients 
        call coeffCal  
       !*******************************
       !
        if(flag==1)    exit   
    end do
    write(*,*)  "Inside loop count:",iter-1
  
end subroutine

subroutine  sourceConvJudge    !judge wheather the source term converge
    implicit none
    integer::i,j
    integer::flag1,flag2,flag3,flag4       !four flags used to verify wheather the four constrains are satisfied

    !write(*,*)  "Source_conv_judge"
    !initialize the four flags  with value '1',which represents that they satisfy the constrains
    flag1=1 ; flag2=1 ; flag3=1 ; flag4=1
    !***********************
    !calculate the real space and angle around the boundaries
    do i=1,imax
        dd(i,0)=sqrt(x_eta(i,1)**2+y_eta(i,1)**2)
        theta(i,0)=abs( acos((x_xi(i,1)*x_eta(i,1)+y_xi(i,1)*y_eta(i,1)) &
                  & /(sqrt(x_xi(i,1)**2+y_xi(i,1)**2)*sqrt(x_eta(i,1)**2+y_eta(i,1)**2)))   )
        dd(i,1)=sqrt(x_eta(i,jmax)**2+y_eta(i,jmax)**2)
        theta(i,1)=abs( acos((x_xi(i,jmax)*x_eta(i,jmax)+y_xi(i,jmax)*y_eta(i,jmax)) &
                  & /(sqrt(x_xi(i,jmax)**2+y_xi(i,jmax)**2)*sqrt(x_eta(i,jmax)**2+y_eta(i,jmax)**2)))  )
    end do
    do i=1,imax
        !verify wheather the inside boundary satisfies the boundary's constrains 
        if(abs(dd(i,0)-dr_i)>eps1)                flag1=0
        if(abs(theta(i,0)-thetar_i)>eps2)         flag2=0
    end do
    do i=1,imax
        !verify wheather the outside boundary satisfies the boundary's constrains 
        if(abs(dd(i,1)-dr_o)>eps1)                flag3=0
        if(abs(theta(i,1)-thetar_o)>eps2)         flag4=0
    end do
    !if the four constrains are satisfied, conv is setted to value '1' 
    if(flag1*flag2*flag3*flag4==1) then
        conv=1             
    else
        conv=0              
    end if
    
end subroutine

subroutine sourceRenew   !renew source term
    implicit none 
    integer::i,j
     
    !renew the source term  according to the real space and the angle
    !inside boundary
    do i=1,imax
        p(i,1)=p(i,1) - sigma*atan(thetar_i-theta(i,0))
        q(i,1)=q(i,1) + sigma*atan(dr_i-dd(i,0))
    end do
 
    !outside boundary
    do i=1,imax
        p(i,jmax)=p(i,jmax) + sigma*atan(thetar_o-theta(i,1))
        q(i,jmax)=q(i,jmax) - sigma*atan(dr_o-dd(i,1))
    end do
 
    !renew the source term of all points with exponential-interpolation method
    do j=2,jmax-1
        do i=1,imax
            p(i,j)=p(i,1)*exp(-a*(j-1)) + p(i,jmax)*exp(-b*(jmax-j))
            q(i,j)=q(i,1)*exp(-c*(j-1)) + q(i,jmax)*exp(-d*(jmax-j))
        end do
    end do
 
end subroutine

end module


   
        
   
   
    