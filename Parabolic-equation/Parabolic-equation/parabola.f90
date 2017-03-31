module paraloba
    use controlPara
    use wall
    implicit none
    
    real(8),allocatable:: p(:,:),q(:,:),x1(:,:),y1(:,:)
    real(8),allocatable:: lk(:,:),theta(:,:),delta_xi(:,:),delta_eta(:,:)
    
    
    
contains

subroutine  initalize
    implicit none
    
    allocate( p(nnodes,my))
    allocate( q(nnodes,my))
    
    allocate( lk(nnodes,my+1))
    allocate( x1(nnodes,my+1))
    allocate( y1(nnodes,my+1))
    
    allocate( delta_eta( nnodes,my+1))
    allocate( delta_xi( nnodes,my+1 ))
    
    delta_eta =1.0
    delta_xi=1.0
    
    do i=1,my
        lk(:,i) = thick * growthRate*i
    end do
    
    theta = angle
    
    p=0.0
    q=0.0
    
    xy(:,:,1) = inode
    
    
end subroutine

subroutine  advanceAlgebra(layer)
    implicit none
    
    integer::layer
    integer::i,j,k
    real(8)::delta_r(2,nnodes),dxy(2),ds
    
    do k=layer+1,layer+3
        do i=1,nnodes
            if(i==1)  then
                dxy= xy(,i+1,k-1)- xy(,nnodes,k-1)
            else if(i==n) then
                dxy= xy(,1,k-1)- xy(,i-1,k-1)
            else
                dxy= xy(,i+1,k-1)- xy(,i-1,k-1)
            end if
            
            ds= sqrt(sum(dxy**2))
            
            xy(1,i,k) = xy(1,i,k-1) + lk(i,k)*( -dxy(2)/ds)
            xy(2,i,k) = xy(2,i,k-1) + lk(i,k)*(  dxy(1)/ds)
            
        end do
    end do
    
         
end subroutine


subroutine getCoeff(layer)       !possion's or laplace's equation cofficinets calculating
    implicit none
    integer::i,j
   
    !-------------------------------------------------
    !the metric of x to xi 
    !the grid points excluding the left and the riht boundaries 
    !adopt central-difference method , 2 order accuracy
    do j=1,jmax  
    do i=2,imax-1
        x_xi(i,j)=(x(i+1,j)-x(i-1,j))/2.0
        y_xi(i,j)=(y(i+1,j)-y(i-1,j))/2.0
    end do
    end do

    !the left and the riht boundaries 
    !adopt the central-difference method
    !notice that the left and the right boundaries  are the same boundary
    do j=1,jmax
       ! x_xi(1,j)=x(2,j)-x(1,j)
       ! y_xi(1,j)=y(2,j)-y(1,j)
       ! x_xi(imax,j)=x(imax,j)-x(imax-1,j)
       ! y_xi(imax,j)=y(imax,j)-y(imax-1,j)
        x_xi(1,j)=(x(2,j)-x(imax-1,j))/2.0
        y_xi(1,j)=(y(2,j)-y(imax-1,j))/2.0
        x_xi(imax,j)=x_xi(1,j)
        y_xi(imax,j)=y_xi(1,j)     
    end do
    
    !the metric of x to eta
    !the grid points excluding the bottom and the top boundaries 
    !adopt central-difference method , 2 order accuracy
    do i=1,imax
    do j=2,jmax-1  
        x_eta(i,j)=(x(i,j+1)-x(i,j-1))/2.0
        y_eta(i,j)=(y(i,j+1)-y(i,j-1))/2.0
    end do
    end do
    !the bottom and the top boundaries 
    !adopt the forward-difference and the backward_difference method,1 order accuracy
    do i=1,imax
        x_eta(i,1)=x(i,2)-x(i,1)
        y_eta(i,1)=y(i,2)-y(i,1)
        x_eta(i,jmax)=x(i,jmax)-x(i,jmax-1)
        y_eta(i,jmax)=y(i,jmax)-y(i,jmax-1)
    end do
    
    !-------------------------------------------------
    !x_xieta,y_xietaµÄ¼ÆËã
    !the metric of x to xi.eta and  y to xi.eta
    !inner grid points
    do j=2,jmax-1  
    do i=2,imax-1 
        x_xieta(i,j)=(x(i-1,j-1)-x(i-1,j+1)-x(i+1,j-1)+x(i+1,j+1))/4.0
        y_xieta(i,j)=(y(i-1,j-1)-y(i-1,j+1)-y(i+1,j-1)+y(i+1,j+1))/4.0
    end do
    end do
    !the left and the right boundaries(exclude the points on the topand the bottom boundaries)
    do j=2,jmax-1
       
        !x_xieta(i,j)=x_eta(2,j)-x_eta(1,j)                               !x_eta to xi,forward-difference method
        !y_xieta(i,j)=y_eta(2,j)-y_eta(1,j)
        
        !x_xieta(1,j)=(x_xi(1,j+1)-x_xi(1,j-1))/2                         !x_xi to eta with central-difference method
        !y_xieta(1,j)=(y_xi(1,j+1)-y_xi(1,j-1))/2
        !x_xieta(imax,j)=(x_xi(imax,j+1)-x_xi(imax,j-1))/2
        !y_xieta(imax,j)=(y_xi(imax,j+1)-y_xi(imax,j-1))/2
        
       
        x_xieta(1,j)=(x(imax-1,j-1)-x(imax-1,j+1)-x(2,j-1)+x(2,j+1))/4.0   !adopt 4 points method
        y_xieta(1,j)=(y(imax-1,j-1)-y(imax-1,j+1)-y(2,j-1)+y(2,j+1))/4.0
        x_xieta(imax,j)=x_xieta(1,j)
        y_xieta(imax,j)=y_xieta(1,j)
    end do 
     !the top and the bottom boundaries(includ the points on the left and the right boundaries)
    do i=1,imax
        x_xieta(i,1)=x_xi(i,2)-x_xi(i,1)
        y_xieta(i,1)=y_xi(i,2)-y_xi(i,1)
        x_xieta(i,jmax)=x_xi(i,jmax)-x_xi(i,jmax-1)
        y_xieta(i,jmax)=y_xi(i,jmax)-y_xi(i,jmax-1)
    end do 
    !-------------------------------------------------------
    !the cofficient 
    
    alpha =x_eta**2 + y_eta**2
    beta  =x_xi*x_eta + y_xi*y_eta
    gamma =x_xi**2 + y_xi**2 
    jacobi=x_xi*y_eta- x_eta*y_xi
    
    !----------------------------------------------------
end subroutine
    
subroutine  iter(layer)
    implicit none
    
    integer::layer
    integer::i,j,k
    integer::flag
    real(8)::dxy(2)
    real(8):;allocatable::alpha(:,:),beta(:,:),gamma(:,:),jacobi(:,:)
    
     do i=1,nnodes
         
         if() then
             alpha(i,k) = xy(i,)
    
    
    
    k=0
    
end subroutine

    
end module
