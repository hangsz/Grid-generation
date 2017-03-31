module controlPara
    
    implicit none
    integer,parameter::Mx=10,My=40
    real(8),parameter::eps=10E-8,eps2=2E-20
    real(8),parameter::thick=2E-5,angle=1.57079663,growthRate=1.2!预设的各层层高和角度，此处可以改进,eps为源项确定之后的迭代精度，eps2为源项的精度
    real(8),parameter::a=0.5,b=0.5,c=0.5,d=0.5,sigama=0.30!0.3   !源项的参数赋值
    real,parameter::steps=300   !光顺
    
end module
