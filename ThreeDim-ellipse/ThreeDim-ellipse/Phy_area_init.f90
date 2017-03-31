 !物理域参数初始化
module Phy_area_init  
    use Control
    use r_oper
   ! use Coor_of_bound
    implicit none
    
    !物面边界坐标向量
    type(vector),save::r(imax,jmax,kmax)
    !物面边界坐标
    real(kind=8),save::x(imax,jmax,kmax),y(imax,jmax,kmax) ,z(imax,jmax,kmax)  
    !源项值，任意一个x,y,z点
    real(kind=8),save::phi_p(imax,jmax,kmax),phi_q(imax,jmax,kmax),phi_r(imax,jmax,kmax)      
    
contains 

subroutine Coor_init    
    !物理域坐标初始化
    implicit none
    integer::i,j,k
    integer,parameter::fileid=10  !文件编号
    character(len=80)::filename  !文件名
    integer::error              !读取数据是否成功判定
    logical::alive              !查询文件是否存在判定
   
    
     !初始化边界坐标
    !x(:,:,1)=x_bound(:,:)
    !y(:,:,1)=y_bound(:,:)
    !z(:,:,1)=z_bound(:,:) 
    
    filename="Coordinate.dat"
    inquire(file=filename,exist=alive)    
    if(.not.alive) then
        write(*,*) trim(filename),"Does't exist!"   
        stop
    end if
    
    open(fileid,file=filename,action='read',status='old')
    do k=1,kmax
    do j=1,jmax
    do i=1,imax
        read(fileid,"(2X,F11.8,2X,F11.8,2X,F11.8)",iostat=error) x(i,j,k),y(i,j,k),z(i,j,k) !按照11位宽度，8位小数的格式读取

        if(error/=0) then
           write(*,*) "error,读取错误！"
           exit
        end if 
    end do
    end do 
    end do
    close(fileid)
    
  
    r(:,:,:)%x=x(:,:,:)
    r(:,:,:)%y=y(:,:,:)
    r(:,:,:)%z=z(:,:,:)
    
end subroutine


subroutine Source_init
    implicit none 
    integer::i,j,k
    
    phi_p(:,:,:)   =0.0
    phi_q(:,:,:)   =0.0
    phi_r(:,:,:)   =0.0
  
end subroutine

end module

