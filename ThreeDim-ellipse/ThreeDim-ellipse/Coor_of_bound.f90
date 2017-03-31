module Coor_of_bound  
    !物面网格划分模块
    use Control
    implicit none
    !边界的网格点物理域坐标
    real(kind=8)::x_bound(imax,jmax,0:1),y_bound(imax,jmax,0:1),z_bound(imax,jmax,0:1) 
    
contains


subroutine Bound               !O型网格，内边界点的坐标的读取
    implicit none
    integer,parameter::fileid=10  !文件编号
    character(len=80)::filename  !文件名
    integer::error              !读取数据是否成功判定
    logical::alive              !查询文件是否存在判定
    integer::i,j,k                  !内边界数据循环量，y=0
    
    !write(*,*)  "Filename:"         
    !read(*,"(A80)") filename       !输入要读取的文件名
    
    !读取内边界网格
    filename="inner.dat"
    inquire(file=filename,exist=alive)    
    if(.not.alive) then
        write(*,*) trim(filename),"Does't exist!"   
        stop
    end if
    
    open(fileid,file=filename,action='read',status='old')
    do k=1,kmax
    do j=1,jmax
    do i=1,imax
        read(fileid,"(F10.8,4X,F9.8)",iostat=error) x_bound(i,j,0),y_bound(i,j,0),z_bound(i,j,0) !按照9位有效数字，8位小数的格式读取

        if(error/=0) then
           write(*,*) "error,读取错误！"
           exit
        end if 
    end do
    end do
    end do
    close(fileid)
    
    !读取外边界网格
     filename="outer.dat"
    inquire(file=filename,exist=alive)    
    if(.not.alive) then
        write(*,*) trim(filename),"Does't exist!"   
        stop
    end if
    
    open(fileid,file=filename,action='read',status='old')
    do k=1,kmax
    do j=1,jmax
    do i=1,imax
        read(fileid,"(F10.8,4X,F9.8)",iostat=error) x_bound(i,j,1),y_bound(i,j,1),z_bound(i,j,1) !按照9位有效数字，8位小数的格式读取

        if(error/=0) then
           write(*,*) "error,读取错误！"
           exit
        end if 
    end do
    end do
    end do
    close(fileid)
      
end subroutine

end module



    
         
        