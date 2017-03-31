module wall   !读入边界条件，三维的是面网格，二维的是边界节点，，，，处理的是非结构网格的存储文件
    use controlPara
    
    implicit none
    
    integer::nnodes,nedges
    real(8),allocatable::inode(:,:)
    integer,allocatable::iedge(:,:)
    
contains

subroutine readWall
    implicit none
    integer:: i
    open(10,file='airfoil/naca0012.dat')
    read(10,*) nnodes,nedges
    
      
    !allocate space for iedge array
    allocate(inode(2,nodes))
    allocate(iedge(2,nedges))
    
    do i=1,nnodes
         read(10,*) inode(:,i)                                          !inode(1,i)表示第i个点的x坐标，inode(2,i)表示第i个点的y坐标
    end do
    
    do i=1,nedges
        read(10,*)  iedge(:,i)                              !iedge(1,i)表示第i条边的第一个顶点序号，iedge(2,i)表示第i条边的第2个顶点标号
    end do 
    
    close(10)
    
end subroutine

end module