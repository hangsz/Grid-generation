!对向量r重新加载操作符或定义操作函数
module r_oper   
    implicit none
    type vector
        real(kind=8)::x
        real(kind=8)::y
        real(kind=8)::z
    end type
    
    interface  operator(+)
      module procedure add
    end interface
    
    interface  operator(-)
      module procedure sub
    end interface
    
    interface  operator(*)
      module procedure multi
    end interface
    
  !  interface  operator(**)
  !    module procedure magn
  !  end interface
    
contains

function add(r1,r2)
    implicit none
    type(vector),intent(in)::r1,r2
    type(vector)::add
    real(kind=8)::xx,yy,zz
    xx=r1%x-r2%x
    yy=r1%y-r2%y
    zz=r1%z-r2%z
    add%x=xx
    add%y=yy
    add%z=zz
    return
end function add

function sub(r1,r2)
    implicit none
    type(vector),intent(in)::r1,r2
    type(vector)::sub
    real(kind=8)::xx,yy,zz
    xx=r1%x-r2%x
    yy=r1%y-r2%y
    zz=r1%z-r2%z
    sub%x=xx
    sub%y=yy
    sub%z=zz
    return
end function sub

function multi(r1,r2)
    implicit none
    type(vector),intent(in)::r1,r2
    real(kind=8)::multi
    multi=r1%x*r2%x + r1%y*r2%y + r1%z*r2%z
    return
end function multi

function pow2(r)
    implicit none
    type(vector),intent(in)::r
    real(kind=8)::pow2
    pow2=r%x **2 + r%y **2 + r%z **2
    return
end function pow2
 

function  ex_product(r1,r2)
    implicit none
    type(vector)::r1,r2
    type(vector)::ex_product
    real(kind=8)::xx,yy,zz
    xx=r1%y * r2%z-r1%z * r2%y
    yy=-(r1%x * r2%z-r1%z * r2%x)
    zz=r1%x * r2%y-r1%y * r2%x
    ex_product%x=xx
    ex_product%y=yy
    ex_product%z=zz
    return
end function ex_product

function  time_product(r,time)
    implicit none
    real(kind=8)::time
    type(vector)::r
    type(vector)::time_product
    real(kind=8)::xx,yy,zz
    xx=time* r%x 
    yy=time* r%y
    zz=time* r%z
    time_product%x=xx
    time_product%y=yy
    time_product%z=zz
    return
end function time_product
  
    
end module