module phyDomainInit     !physical  domain  initialization
    use boundCoor    !inherit module Coor_of_bound
    implicit none
    real(kind=8)::x(imax,jmax),y(imax,jmax)   !x,y coordinates of grid vertex
    real(kind=8)::p(imax,jmax),q(imax,jmax)   !source term  p:control the orthogonality
                                              !             q::conrtrol the space of two layers of grid
  
contains 
subroutine coorInit    !initialize the coordinates of the physical area
    implicit none
    integer::i,j
    
    write(*,*)  "coorInit"
    
    !the coordinates of the boundary initialization
    x(:,1)=x_in
    y(:,1)=y_in
    x(:,jmax)=x_out
    y(:,jmax)=y_out
      
    !the inner coordinates  initialized by linear interpolation method
    do i=1,imax
    do j=2,jmax-1
        x(i,j)=x(i,1)+(x(i,jmax)-x(i,1))/(jmax-1)*(j-1)
        y(i,j)=y(i,1)+(y(i,jmax)-y(i,1))/(jmax-1)*(j-1)
        !the other  means of expression
        !x(i,j)=(jmax-j)/(jmax-1)*x(i,1)+(j-1)/(jmax-1)*x(i,jmax)
        !y(i,j)=(jmax-j)/(jmax-1)*y(i,1)+(j-1)/(jmax-1)*y(i,jmax)
    end do
    end do
    
end subroutine

subroutine sourceInit              !source teram initialization
    implicit none 
    integer::i,j
    write(*,*)  "sourceInit"
    
    !set the source term of inside and outside boundaries to zero
    p(:,1)   =0.0
    q(:,1)   =0.0
    p(:,jmax)=0.0
    q(:,jmax)=0.0
     
    !********************************************
    !initialize the source term of inner grid vertex  by exponential interpolation method
    
    do i=1,imax
    do j=2,jmax-1
        p(i,j)=p(i,1)*exp(-a*(j-1)) + p(i,jmax)*exp(-c*(jmax-j))
        q(i,j)=q(i,1)*exp(-b*(j-1)) + q(i,jmax)*exp(-d*(jmax-j))
    end do
    end do
     
end subroutine

end module

