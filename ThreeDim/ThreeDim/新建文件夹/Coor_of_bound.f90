module Coor_of_bound  
    !物面网格划分模块
    use Control
    implicit none
    !real(kind=8),parameter::pi=3.14159265
    !integer,parameter::imax=133,jmax=60    !计算域的最大网格点数,（o型网格，imax为奇数）
    !integer,parameter::imax=201,jmax=60   !计算域的最大网格点数,（c型网格，imax为奇数）
    !c型网格区域
    !real(kind=8),parameter::x_left=-6,x_right=6    !物理域x向左右范围
    !real(kind=8),parameter::y_bottom=-4,y_top=4    !物理域y向上下范围
    !o型网格区域
    !real(kind=8),parameter::R_outer=6.0
    real(kind=8)::x_bound(imax,jmax),y_bound(imax,jmax),z_bound(imax,jmax)  !边界的网格点物理域坐标(二维)
    
contains


subroutine Inner_O                !O型网格，内边界点的坐标的读取
    implicit none
    integer,parameter::fileid=10  !文件编号
    character(len=80)::filename  !文件名
    integer::error              !读取数据是否成功判定
    logical::alive              !查询文件是否存在判定
    integer::i                  !内边界数据循环量，y=0
    
    !write(*,*)  "Filename:"         
    !read(*,"(A80)") filename       !输入要读取的文件名

    filename="naca0012_O.txt"
    inquire(file=filename,exist=alive)    
    if(.not.alive) then
        write(*,*) trim(filename),"Does't exist!"   
        stop
    end if
    
    open(fileid,file=filename,action='read',status='old')
    do i=1,imax
        read(fileid,"(F10.8,4X,F9.8)",iostat=error) x_inner(i),y_inner(i)  !按照9位有效数字，8位小数的格式读取
        if(error/=0) then
           write(*,*) "error,读取错误！"
           exit
        end if 
     !   x_outer(i)= x_outer(i)-Center
    end do
    close(fileid)
   
      
end subroutine

subroutine Outer_O          !这里不需要读取，需要进行几何划分，采用等角度划分
    implicit none           !笛卡尔坐标在物理域几何中心
    integer::i
    
    do i=1,imax
        x_outer(i)=R_outer*cos(2*pi/(imax-1)*(i-1))
        y_outer(i)=R_outer*sin(2*pi/(imax-1)*(i-1))
    end do
    
end subroutine

subroutine Inner_C               !O型网格，内边界点的坐标的读取
    implicit none
    integer::airfoil_num=131
    integer::traiedge_num
    integer,parameter::fileid=10  !文件编号
    character(len=80)::filename   !文件名
    integer::error                !读取数据是否成功判定
    logical::alive                !查询文件是否存在判定
    integer::i                    !内边界数据循环量，y=0
    
    !内边界后缘到右边界的插值(采用线性插值，有缺点，可以采用指数插值)
    traiedge_num=(imax-airfoil_num)/2
    do i=1,traiedge_num
        x_inner(i)=x_right-(x_right-1.0084)/(traiedge_num-1)*(i-1)!减去1.01是蒋翼型稍微右延长一点
        y_inner(i)=0.0
        x_inner(imax+1-i)=x_inner(i)
        y_inner(imax+1-i)=y_inner(i)
    end do
    
    !读入翼型边界条件
    !write(*,*)  "Filename:"         
    !read(*,"(A80)") filename       !输入要读取的文件名
    filename="naca0012_C.txt"
    
    inquire(file=filename,exist=alive)    
    if(.not.alive) then
        write(*,*) trim(filename),"Does't exist!"   
        stop
    end if
    
    open(fileid,file=filename,action='read',status='old')
    do i=traiedge_num+1,traiedge_num+airfoil_num
        read(fileid,"(F10.8,4X,F11.8)",iostat=error) x_inner(i),y_inner(i)  !按照9为有效数字，8位小数的格式读取
        if(error/=0) then
           write(*,*) "error,读取错误！"
           exit
        end if
    end do
    close(fileid)
    
    
end subroutine

subroutine Outer_C            !这里不需要读取，需要进行几何划分，采用等弧度划分
    implicit none           !笛卡尔坐标在物理域几何中心
    !real,parameter::r=2    !圆弧半径
    integer::arc_num     !单个圆弧部分分配的坐标点个数
    integer::bottom_num  !底部外边界坐标点个数
    integer::left_num  !左下部外边界坐标点个数
    integer::i,i2,i3,i4,i5
    arc_num=(imax+1)/2*(0.5*pi*r/(abs(x_right)+abs(x_left)-r+abs(y_bottom)-r+0.5*pi*r))  !imax+1之后，奇偶都成立
    left_num=anint(((imax+1)/2-arc_num)*((abs(y_bottom)-r)/(abs(x_right)+abs(x_left)-r+abs(y_bottom)-r)))  !向最接近的整数取整
    bottom_num=(imax+1)/2-arc_num-left_num
    i2=bottom_num
    i3=bottom_num+arc_num
    i4=bottom_num+arc_num+left_num
    !write(*,*) i2,i3,i4,bottom_num,arc_num,left_num
    do i=1,i2
        x_outer(i)=x_right-(i-1)*(abs(x_right)+abs(x_left)-r)/bottom_num
        y_outer(i)=y_top
    end do
    do i=i2+1,i3
        x_outer(i)=x_left+r-r*sin(pi/4.0/arc_num*(i-i2))
        y_outer(i)=y_top-r+r*cos(pi/4.0/arc_num*(i-i2))
    end do
    do i=i3+1,i4
        x_outer(i)=x_left
        y_outer(i)=y_top-r-(abs(y_top)-r)/left_num*(i-i3)
    end do
    
    !由于物理域和物理边界对称，所以采取对称赋值
    do i=imax,imax/2,-1            
        x_outer(i)=x_outer(imax-i+1)
        y_outer(i)=-y_outer(imax-i+1)
    end do
    
end subroutine

end module



    
         
        