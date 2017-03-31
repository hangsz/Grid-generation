
include 'Control.f90'
!include 'Coor_of_bound.f90'
include 'r_oper.f90'
include 'Phy_area_init.f90'
include 'Solution.f90'
include 'Print_data.f90'


program main 
    use Print_data
    implicit none
    
    call Solve  !生成网格
    write(*,*) "输出结果："
    call Output
    
   
end program
              