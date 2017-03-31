

 
program main 
    use Print_data
    implicit none
    
    call solve  !生成网格
    write(*,*) "输出结果："
    call output
    
   
end program
              