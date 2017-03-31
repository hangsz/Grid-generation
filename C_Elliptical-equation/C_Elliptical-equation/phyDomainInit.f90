module phyDomainInit
    use boundCoor
    implicit none
    real(8),dimension(imax,jmax)::x,y , &    !coordinates 
                                  p,q        !source term
                                      
contains 

subroutine initialize

    call inside             
    call outside 
    call leftRight
    call coorInit            
    call sourceInit  

end subroutine

subroutine coorInit    !物理域坐标初始化
    implicit none
    integer::i,j
    
     !初始化边界坐标
      
       x(1,:) = x_l  ; y(1,:) = y_l
       x(imax,:) = x_r ; y(imax,:)  = y_r
       x(:,1)=x_inner    ;  y(:,1)=y_inner
       x(:,jmax)=x_outer ;  y(:,jmax)=y_outer
       
       x(1,1) = x_right ;y(1,1) = 0.0
       x(1,jmax) = x_right ;y(1,jmax) = y_top
       x(imax,jmax) = x_right ; y(imax,jmax) = y_bottom
    !通过线性插值，确定内部点的坐标
    do i=1,imax
        do j=2,jmax-1
            x(i,j)=x(i,1)+(x(i,jmax)-x(i,1))/(jmax-1)*(j-1)
            y(i,j)=y(i,1)+(y(i,jmax)-y(i,1))/(jmax-1)*(j-1)
            !x(i,j)=(jmax-j)/(jmax-1)*x(i,1)+(j-1)/(jmax-1)*x(i,jmax)
            !y(i,j)=(jmax-j)/(jmax-1)*y(i,1)+(j-1)/(jmax-1)*y(i,jmax)
        end do
    end do
    
end subroutine

subroutine sourceInit
    implicit none 
    integer::i,j

    p(:,1)   =0.0 ; q(:,1)   =0.0
    p(:,jmax)=0.0 ; q(:,jmax)=0.0

    !********************************************
    !任意一点源项值
    
    do i=1,imax
    do j=2,jmax-1
        p(i,j)=p(i,1)*exp(-a*(j-1)) + p(i,jmax)*exp(-c*(jmax-j))
        q(i,j)=q(i,1)*exp(-b*(j-1)) + q(i,jmax)*exp(-d*(jmax-j))
    end do
    end do
    
end subroutine

end module

