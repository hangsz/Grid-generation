
module Solution
    !求解模块
    use Phy_area_init
    implicit none
    
    type(vector),save::r_xi(imax,jmax,kmax),r_eta(imax,jmax,kmax),r_zeta(imax,jmax,kmax)!,r_xieta(imax,jmax,kmax)
    !泊松方程系数
    real(kind=8),save::alpha1(imax,jmax,kmax),alpha2(imax,jmax,kmax),alpha3(imax,jmax,kmax),beta12(imax,jmax,kmax)
    !求解泊松方程需用的差分值
    real(kind=8),save::x_xi(imax,jmax,kmax),y_xi(imax,jmax,kmax),z_xi(imax,jmax,kmax)  
    real(kind=8),save::x_eta(imax,jmax,kmax),y_eta(imax,jmax,kmax),z_eta(imax,jmax,kmax)!
    real(kind=8),save::x_zeta(imax,jmax,kmax),y_zeta(imax,jmax,kmax),z_zeta(imax,jmax,kmax)
    real(kind=8),save::x_xieta(imax,jmax,kmax),y_xieta(imax,jmax,kmax),z_xieta(imax,jmax,kmax)
   
    !n+1次迭代得到的内点坐标
    real(kind=8),save::x1(imax,jmax,kmax),y1(imax,jmax,kmax),z1(imax,jmax,kmax)
    !网格推进步长（实际值和给定值）
    real(kind=8)::Lk(imax,jmax,kmax)
    !有源项时的，边界处的间距和夹角(理想值和真实值)
    real(kind=8)::dr(imax,jmax,0:1)!theta_rxi,theta_reta  !,dd(imax,jmax,0:1),theta_xz(imax,jmax,0:1),theta_ez(imax,jmax,0:1)
    !凸角修正时某点某方向的凸角
    real(kind=8)::theta_xi(imax,jmax,kmax),theta_eta(imax,jmax,jmax)
    !凹角修正的系数
    real(kind=8)::lamda_xi(imax,jmax,kmax),lamda_eta(imax,jmax,kmax)
    !Ellipse收敛程序的返回值，用于控制循环是否终止
    integer,save::conv      
    !生成每层网格的循环次数         
    integer::iter                  
    integer::count
contains

subroutine solve
    !网格迭代计算程序
    implicit none
   integer::i,j,k
   
   !*******************************************************
   !初始化需用数据
   call Coor_Init 
   call Inc_len_Init
   call Source_Init
   
 !  write(*,*)  "44"
 !  pause
   
   !*********************************************************
    !从k=2层网格开始生成
 do k=2,kmax-1!(k-1保证不越界)
     
     
     count=0
    !生成第一层（k层）代数网格，先进行推进步长修正，再生成网格
    call Convex_amend(k)   
   
    call Algebra_Grid(k)
    
   !每次都要更新第二层网格
   iter=0
   do while(iter /=itermax)
       
    iter=iter+1
       
    !生成第二层（k+1层）代数网格，先进行推进步长修正，再生成网格
    call Convex_amend(k+1) 
   
    
    call Algebra_Grid(k+1)
    
     
   
    !计算椭圆方程中的系数
    call Coeff_Cal_t(k)
    
    call Coeff_Cal_n(k)
     
    
    
    !计算椭圆方程中的源项
    call Source_cal(k)

    
    !进行椭圆方程的凹角修正的两个系数计算
    call Concave_amend(k)
    
    !对第一层网格进行光顺,用椭圆形微分方程进行光顺，只计算一次
    call Ellipse(k)
    count=count+1
    write(*,*)  count
    
    !第k层网格收敛性判定，如果收敛，就跳出循环内层循环，开始生成 k+1 层网格
     if(conv==1) then
         write(*,*) "第",k,"层网格满足要求！" 
         exit
     else
         write(*,*) "第",k,"层网格未满足要求，继续！" 
     end if
   end do
   
 end do
 
end subroutine

subroutine Convex_amend(k)
    !凸角修正
    !k表示生成到第n层网格
    
    implicit none
    real(kind=8)::a,b
    integer::k
    integer::i,j
    
    !k层推进步长都设为初始值
    Lk(:,:,k)=lk_init(k)
   ! do j=1,jmax
   ! do i=1,imax
    !     write(*,*)  Lk(i,j,k)
   ! end do
    !end do
   ! stop 
    !对内点进行修正
    do j=2,jmax-1
    do i=2,imax-1
          theta_xi(i,j,k-1) =acos(  (r(i+1,j,k-1)-r(i,j,k-1))*(r(i-1,j,k-1)-r(i,j,k-1))/(sqrt(pow2(r(i+1,j,k-1)-r(i,j,k-1)))*sqrt(pow2(r(i-1,j,k-1)-r(i,j,k-1))))  )
          theta_eta(i,j,k-1)=acos(  (r(i,j+1,k-1)-r(i,j,k-1))*(r(i,j-1,k-1)-r(i,j,k-1))/(sqrt(pow2(r(i,j+1,k-1)-r(i,j,k-1)))*sqrt(pow2(r(i,j-1,k-1)-r(i,j,k-1))))  )
          a=sin(theta_xi(i,j,k-1)/2)
          b=sin(theta_eta(i,j,k-1)/2)
          if(a<b) then
              Lk(i,j,k)=Lk(i,j,k)*a
          else
              Lk(i,j,k)=Lk(i,j,k)*b
          end if  
       ! write(*,*) r(i,j,k-1)%x,r(i,j,k-1)%y,r(i,j,k-1)%z  
       ! write(*,*)  Lk(i,j,k),theta_xi(i,j,k-1),theta_eta(i,j,k-1)!
    end do
    end do
end subroutine

subroutine Algebra_Grid(k)
    !代数法生产两层网格
    implicit none
    integer::i,j,k
   
        !代数法推进生产一层网格，运动凸角修正
        call Convex_amend(k)  !先进行凸角修正 
     
        call Coeff_Cal_t(k-1)    !再进行切向量计算计算 
      
        !生成第一层网格额
  
        do j=1,jmax
        do i=1,imax
      ! do j=2,jmax-1
      ! do i=2,imax-1
         !r(i,j,k)=r(i,j,k-1) + time_product( ex_product(  (r(i+1,j,k-1)-r(i-1,j,k-1)),(r(i,j+1,k-1)-r(i,j-1,k-1))  ) ,Lk(i,j,k) / sqrt(  pow2( ex_product(  (r(i+1,j,k-1)-r(i-1,j,k-1)),(r(i,j+1,k-1)-r(i,j-1,k-1))  ) ) )  )
         r(i,j,k)=r(i,j,k-1) + time_product( ex_product(  r_xi(i,j,k-1),r_eta(i,j,k-1) )  ,Lk(i,j,k) / sqrt(  pow2( ex_product(  r_xi(i,j,k-1),r_eta(i,j,k-1)  ) ) )  )  
        !write(*,*) r(i,j,k)%x,r(i,j,k)%y,r(i,j,k)%z
        !write(*,*) r_xi(i,j,k-1)%x,r_xi(i,j,k-1)%y,r_xi(i,j,k-1)%z
        !write(*,*) r_eta(i,j,k-1)%x,r_eta(i,j,k-1)%y,r_eta(i,j,k-1)%z
        end do
        end do
           
        !上下边界（包含两端点）
        !do i=1,imax
        !r(i,1,k)=r(i,1,k-1)+ Lk(i,1,k)
        !r(i,jmax,k)=r(i,jmax,k-1) + Lk(i,jmax,k)   
        !end do
        !do  j=2,jmax-1
        !r(1,j,k)=r(1,j,k-1) + Lk(1,j,k)
        !r(imax,j,k)=r(imax,j,k-1) + Lk(imax,j,k)
        !end do
        
        !这里要进行反赋值，坐标和向量存储空间不同
       do j=1,jmax
       do i=1,imax
          x(i,j,k)=r(i,j,k)%x
          y(i,j,k)=r(i,j,k)%y
          z(i,j,k)=r(i,j,k)%z
       end do
       end do
       
end subroutine


subroutine Coeff_Cal_t(k)
   !面切向的系数计算
   !欧拉或泊松方程的系数计算
   implicit none
   integer::k
   integer::i,j
   
    !需用的差分计算
    !*******************************************  
    !xi
    do j=1,jmax  
    do i=2,imax-1
        x_xi(i,j,k)=(x(i+1,j,k)-x(i-1,j,k))/2
        y_xi(i,j,k)=(y(i+1,j,k)-y(i-1,j,k))/2
        z_xi(i,j,k)=(z(i+1,j,k)-z(i-1,j,k))/2
    end do
    end do
       !左右边界点对xi的差分
    do j=1,jmax
        x_xi(1,j,k)=x(2,j,k)-x(1,j,k)
        y_xi(1,j,k)=y(2,j,k)-y(1,j,k)
        z_xi(1,j,k)=z(2,j,k)-z(1,j,k)
        x_xi(imax,j,k)=x(imax,j,k)-x(imax-1,j,k)
        y_xi(imax,j,k)=y(imax,j,k)-y(imax-1,j,k)
        z_xi(imax,j,k)=z(imax,j,k)-z(imax-1,j,k)
    end do
    
  
     !*******************************************  
    !eta 
    do i=1,imax 
    do j=2,jmax-1  
        x_eta(i,j,k)=(x(i,j+1,k)-x(i,j-1,k))/2
        y_eta(i,j,k)=(y(i,j+1,k)-y(i,j-1,k))/2
        z_eta(i,j,k)=(z(i,j+1,k)-z(i,j-1,k))/2
    end do
    end do
      !上下边界点对eta的差分
    do i=1,imax
        x_eta(i,1,k)=x(i,2,k)-x(i,1,k)
        y_eta(i,1,k)=y(i,2,k)-y(i,1,k)
        z_eta(i,1,k)=z(i,2,k)-z(i,1,k)
        x_eta(i,jmax,k)=x(i,jmax,k)-x(i,jmax-1,k)
        y_eta(i,jmax,k)=y(i,jmax,k)-y(i,jmax-1,k)
        z_eta(i,jmax,k)=z(i,jmax,k)-z(i,jmax-1,k)
   end do
    r_xi(:,:,k)%x=x_xi(:,:,k)
    r_xi(:,:,k)%y=y_xi(:,:,k)
    r_xi(:,:,k)%z=z_xi(:,:,k)
    
    r_eta(:,:,k)%x=x_eta(:,:,k)
    r_eta(:,:,k)%y=y_eta(:,:,k)
    r_eta(:,:,k)%z=z_eta(:,:,k)
    
end subroutine

subroutine Coeff_Cal_n(k)
    implicit none
    integer::k
    integer::i,j
    !面法向的系数计算
    !*******************************************  
    !zeta
    do j=1,jmax  
    do i=1,imax
        x_zeta(i,j,k)=(x(i,j,k+1)-x(i,j,k-1))/2
        y_zeta(i,j,k)=(y(i,j,k+1)-y(i,j,k-1))/2
        z_zeta(i,j,k)=(z(i,j,k+1)-z(i,j,k-1))/2
    end do
    end do
    !******************************************* 
    !xieta
    do j=2,jmax-1  
    do i=2,imax-1
        x_xieta(i,j,k)=(x(i-1,j-1,k)-x(i-1,j+1,k)-x(i+1,j-1,k)+x(i+1,j+1,k))/4
        y_xieta(i,j,k)=(y(i-1,j-1,k)-y(i-1,j+1,k)-y(i+1,j-1,k)+y(i+1,j+1,k))/4
        z_xieta(i,j,k)=(z(i-1,j-1,k)-z(i-1,j+1,k)-z(i+1,j-1,k)+z(i+1,j+1,k))/4
    end do
    end do
    !上下边界的处理（不包括与左右边界交点）
    do i=2,imax-1
       x_xieta(i,1,k)=(x_eta(i+1,j,k)-x_eta(i-1,j,k))/2
       y_xieta(i,1,k)=(y_eta(i+1,j,k)-y_eta(i-1,j,k))/2
       z_xieta(i,1,k)=(z_eta(i+1,j,k)-z_eta(i-1,j,k))/2
       x_xieta(i,jmax,k)=(x_eta(i+1,j,k)-x_eta(i-1,j,k))/2
       y_xieta(i,jmax,k)=(y_eta(i+1,j,k)-y_eta(i-1,j,k))/2
       z_xieta(i,jmax,k)=(z_eta(i+1,j,k)-z_eta(i-1,j,k))/2
    end do
    !左右边界的处理（包括交点）
    do j=1,jmax
        x_xieta(1,j,k)=x_eta(2,j,k)-x_eta(1,j,k)
        y_xieta(1,j,k)=y_eta(2,j,k)-y_eta(1,j,k)
        z_xieta(1,j,k)=z_eta(2,j,k)-z_eta(1,j,k)
        x_xieta(imax,j,k)=x_eta(imax,j,k)-x_eta(imax-1,j,k)
        y_xieta(imax,j,k)=y_eta(imax,j,k)-y_eta(imax-1,j,k)
        z_xieta(imax,j,k)=z_eta(imax,j,k)-z_eta(imax-1,j,k)
    end do
    !*************************************************
    !差分结果进行向量化，便于运算
   
    
    r_zeta(:,:,k)%x=x_zeta(:,:,k)
    r_zeta(:,:,k)%y=y_zeta(:,:,k)
    r_zeta(:,:,k)%z=z_zeta(:,:,k)
    
 !   r_xieta(:,:,k)%x=x_xieta(:,:,k)
  !  r_xieta(:,:,k)%y=y_xieta(:,:,k)
  !  r_xieta(:,:,k)%z=z_xieta(:,:,k)
    !*************************************************
    !所有点上的几个参数计算
    do j=1,jmax
    do i=1,imax
    alpha1(i,j,k) =pow2(r_eta(i,j,k) )  *  pow2(r_zeta(i,j,k) )
    alpha2(i,j,k) =pow2(r_zeta(i,j,k)) *  pow2(r_xi(i,j,k))
    alpha3(i,j,k) =pow2(r_xi(i,j,k)) *  pow2(r_eta(i,j,k))  - ( r_xi(i,j,k) * r_eta(i,j,k) )**2
    beta12(i,j,k) =( r_xi(i,j,k)*r_zeta(i,j,k) )* ( r_zeta(i,j,k)*r_eta(i,j,k) )
   ! write(*,*)  pow2(r_eta(i,j,k) ) , pow2(r_zeta(i,j,k) )
    end do
    end do
    !************************************************
    
    !do j=1,jmax
    !do i=1,imax  
    !write(*,*) alpha1(i,j,k),alpha2(i,j,k),alpha3(i,j,k),beta12(i,j,k)
    !end do
    !end do
    !write(*,*)  "521"
    !stop
   
end subroutine     

subroutine Concave_amend(k)
   !凹角修正，求出每个面上的修正系数
   implicit none
   integer::k
   integer::i,j
   real(kind=8)::a,b,c,f1_phi,f2_phi
   real(kind=8)::T=0.8 !!!!!!!!!!!!!!!!!!!!!!!!!!!设置的参数
   
   do j=2,jmax-1
   do i=2,imax-1
       a=pow2(r_xi(i,j,k))
       b=pow2(r_eta(i,j,k))
       c=pow2(r_zeta(i,j,k))
       
       if (a<c)  a=c
       if (b<c)  b=c
       
       if(theta_xi(i,j,k)>=0 .and. theta_xi(i,j,k)<pi/2) then
           f1_phi=1
       else 
           if(theta_xi(i,j,k)>=pi/2 .and. theta_xi(i,j,k)<=pi) then
               f1_phi=(sin(theta_xi(i,j,k)))**T
           else
               f1_phi=0.0
            end if
       end if
       
       if(theta_eta(i,j,k)>=0 .and. theta_eta(i,j,k)<pi/2) then
           f2_phi=1
       else 
           if(theta_eta(i,j,k)>=pi/2 .and. theta_eta(i,j,k)<=pi) then
               f2_phi=(sin(theta_eta(i,j,k)))**T
           else
               f2_phi=0.0
            end if
       end if
       
       lamda_xi(i,j,k)=sqrt(a/c)*f1_phi
       lamda_eta(i,j,k)=sqrt(b/c)*f2_phi
       
   end do
   end do
   
        !do j=2,jmax-1
        !do i=2,imax-1  
         !   write(*,*) theta_xi(i,j,k),theta_eta(i,j,k),lamda_xi(i,j,k),lamda_eta(i,j,k)
        !end do
        !end do
        !write(*,*)  "382"
        !stop
   
end subroutine

 subroutine Ellipse(k)
    implicit none 
    integer::flag
    integer::k
    integer::i,j
    !***************************************************  
    !椭圆法求第k层网格的n+1次坐标新值
    
 do j=2,jmax-1
 do i=2,imax-1
     x1(i,j,k)=0.5/(alpha1(i,j,k)*(1+lamda_xi(i,j,k))+alpha2(i,j,k)*(1+lamda_eta(i,j,k))+alpha3(i,j,k))*&
               &( alpha1(i,j,k)*(1+lamda_xi(i,j,k))*(x(i-1,j,k)+x(i+1,j,k))+alpha1(i,j,k)*phi_p(i,j,k)*x_xi(i,j,k) + &
               &  alpha2(i,j,k)*(1+lamda_eta(i,j,k))*(x(i,j-1,k)+x(i,j+1,k))+alpha2(i,j,k)*phi_q(i,j,k)*x_eta(i,j,k)+ &
               &  alpha3(i,j,k)*(x(i,j,k-1)+x(i,j,k+1))+alpha3(i,j,k)*phi_r(i,j,k)*x_zeta(i,j,k)  + 2*beta12(i,j,k)*x_xieta(i,j,k)  )
     y1(i,j,k)=0.5/(alpha1(i,j,k)*(1+lamda_xi(i,j,k))+alpha2(i,j,k)*(1+lamda_eta(i,j,k))+alpha3(i,j,k))*&
               &( alpha1(i,j,k)*(1+lamda_xi(i,j,k))*(y(i-1,j,k)+y(i+1,j,k))+alpha1(i,j,k)*phi_p(i,j,k)*y_xi(i,j,k) + &
               &  alpha2(i,j,k)*(1+lamda_eta(i,j,k))*(y(i,j-1,k)+y(i,j+1,k))+alpha2(i,j,k)*phi_q(i,j,k)*y_eta(i,j,k)+ &
               &  alpha3(i,j,k)*(y(i,j,k-1)+y(i,j,k+1))+alpha3(i,j,k)*phi_r(i,j,k)*y_zeta(i,j,k)  + 2*beta12(i,j,k)*y_xieta(i,j,k)  )
     z1(i,j,k)=0.5/(alpha1(i,j,k)*(1+lamda_xi(i,j,k))+alpha2(i,j,k)*(1+lamda_eta(i,j,k))+alpha3(i,j,k))*&
               &( alpha1(i,j,k)*(1+lamda_xi(i,j,k))*(z(i-1,j,k)+z(i+1,j,k))+alpha1(i,j,k)*phi_p(i,j,k)*z_xi(i,j,k) + &
               &  alpha2(i,j,k)*(1+lamda_eta(i,j,k))*(z(i,j-1,k)+z(i,j+1,k))+alpha2(i,j,k)*phi_q(i,j,k)*z_eta(i,j,k)+ &
               &  alpha3(i,j,k)*(z(i,j,k-1)+z(i,j,k+1))+alpha3(i,j,k)*phi_r(i,j,k)*z_zeta(i,j,k)  + 2*beta12(i,j,k)*z_xieta(i,j,k)  ) 
 end do
 end do 
        do j=2,jmax-1
        do i=2,imax-1  
       !    write(*,*) x1(i,j,k),y1(i,j,k),z1(i,j,k),"***"
       !    write(*,*) x(i,j,k),y(i,j,k),z(i,j,k)  
        end do
        end do
        write(*,*)  "521"
      !  stop
  
    !*************************************************************
    !收敛判断
    !出现没有满足精度要求的坐标点，就把标记设为0
     conv=1
  do j=2,jmax-1
  do i=2,imax-1
      !未达到最大次数,判定是否满足收敛要求
      if(abs(x1(i,j,k)-x(i,j,k))>eps .or. abs(y1(i,j,k)-y(i,j,k))>eps .or.  abs(z1(i,j,k)-z(i,j,k))>eps)  then
         conv=0
         exit
      end if   
  end do
  end do
  
  !边界上新坐标为
  !上下边界（包括端点）
  do i=1,imax
      x1(i,1,k)=x(i,1,k)
      y1(i,1,k)=y(i,1,k)
      z1(i,1,k)=z(i,1,k)
      x1(i,jmax,k)=x(i,jmax,k)
      y1(i,jmax,k)=y(i,jmax,k)
      z1(i,jmax,k)=z(i,jmax,k)
  end do
  !左右边界
  do j=2,jmax-1
      x1(1,j,k)=x(1,j,k)
      y1(1,j,k)=y(1,j,k)
      z1(1,j,k)=z(1,j,k)
      x1(imax,j,k)=x(imax,j,k)
      y1(imax,j,k)=y(imax,j,k)
      z1(imax,j,k)=z(imax,j,k)
  end do
  
  
    
   !*********************************************** 
   !更新第n层的的坐标
     x(:,:,k)=x(:,:,k)+omg*(x1(:,:,k)-x(:,:,k))
     y(:,:,k)=y(:,:,k)+omg*(y1(:,:,k)-y(:,:,k))
     z(:,:,k)=z(:,:,k)+omg*(z1(:,:,k)-z(:,:,k))

end subroutine

subroutine  Source_cal(k)
    implicit none
    integer::k
    integer::i,j
    real(kind=8)::dd(imax,jmax,0:1),theta_xz(imax,jmax,0:1),theta_ez(imax,jmax,0:1)
!   type(vector)::rr_xi(imax,jmax,kmax),rr_eta(imax,jmax,kmax),rr_zeta(imax,jmax,kmax)
   
     !***********************
     !求解边界处约束真实值,角度需先求差分
     
     !xi和eta计算
     call Coeff_Cal_t(k-1)
     call Coeff_Cal_t(k+1)
  !************************************************
  !zeta
  do j=1,jmax
  do i=1,imax
     x_zeta(i,j,k-1)=x(i,j,k)-x(i,j,k-1)
     y_zeta(i,j,k-1)=y(i,j,k)-y(i,j,k-1)
     z_zeta(i,j,k-1)=z(i,j,k)-z(i,j,k-1)
     x_zeta(i,j,k+1)=x(i,j,k+1)-x(i,j,k)
     y_zeta(i,j,k+1)=y(i,j,k+1)-y(i,j,k)
     z_zeta(i,j,k+1)=z(i,j,k+1)-z(i,j,k)
  end do
  end do
  !***************************************************
    !向量化处理
     r_zeta(:,:,k-1)%x=x_zeta(:,:,k-1) 
     r_zeta(:,:,k-1)%y=y_zeta(:,:,k-1) 
     r_zeta(:,:,k-1)%z=z_zeta(:,:,k-1)
     r_zeta(:,:,k+1)%x=x_zeta(:,:,k+1) 
     r_zeta(:,:,k+1)%y=y_zeta(:,:,k+1) 
     r_zeta(:,:,k+1)%z=z_zeta(:,:,k+1)
  !****************************************************  
  !边界距离和夹角实际值进行计算  
  !k-1 内边界
  do j=1,jmax
  do i=1,imax
     !dd(i,j,0)=sqrt((r(i,j,k+1)-r(i,j,k))**2)
     !theta_xz(i,j,0)=acos((x_xi(i,j,k)*x_zeta(i,j,k))+y_xi(i,j,k)*y_zeta(i,j,k)+z_xi(i,j,k)*z_zeta(i,j,k))/(sqrt(x_xi(i,j,k)**2+y_xi(i,j,k)**2+z_xi(i,j,k)**2)+sqrt(x_zeta(i,j,k)**2+y_zeta(i,j,k)**2+z_zeta(i,j,k)**2)))
     !theta_ez(i,j,0)=acos((x_eta(i,j,k)*x_zeta(i,j,k))+y_eta(i,j,k)*y_zeta(i,j,k)+z_eta(i,j,k)*z_zeta(i,j,k))/(sqrt(x_eta(i,j,k)**2+y_eta(i,j,k)**2+z_eta(i,j,k)**2)+sqrt(x_zeta(i,j,k)**2+y_zeta(i,j,k)**2+z_zeta(i,j,k)**2)))
     dd(i,j,0)=sqrt(pow2(r_zeta(i,j,k-1)))
     theta_xz(i,j,0)=acos((r_xi(i,j,k-1)*r_zeta(i,j,k-1))/(sqrt(pow2(r_xi(i,j,k-1)))*sqrt(pow2(r_zeta(i,j,k-1))) ) )
     theta_ez(i,j,0)=acos((r_eta(i,j,k-1)*r_zeta(i,j,k-1))/(sqrt(pow2(r_eta(i,j,k-1)))*sqrt(pow2(r_zeta(i,j,k-1))) ) )
     
     
    
  end do
  end do
 
  !k+1 外边界
  do j=1,jmax
  do i=1,imax  
     !dd(i,j,1)=sqrt((r(i,j,k)-r(i,j,k-1))**2)
     !theta_xz(i,j,1)=acos((x_xi(i,j,k)*x_zeta(i,j,k))+y_xi(i,j,k)*y_zeta(i,j,k)+z_xi(i,j,k)*z_zeta(i,j,k))/(sqrt(x_xi(i,j,k)**2+y_xi(i,j,k)**2+z_xi(i,j,k)**2)+sqrt(x_zeta(i,j,k)**2+y_zeta(i,j,k)**2+z_zeta(i,j,k)**2)))
     !theta_ez(i,j,1)=acos((x_eta(i,j,k)*x_zeta(i,j,k))+y_eta(i,j,k)*y_zeta(i,j,k)+z_eta(i,j,k)*z_zeta(i,j,k))/(sqrt(x_eta(i,j,k)**2+y_eta(i,j,k)**2+z_eta(i,j,k)**2)+sqrt(x_zeta(i,j,k)**2+y_zeta(i,j,k)**2+z_zeta(i,j,k)**2)))
     dd(i,j,1)=sqrt(pow2(r_zeta(i,j,k+1)))
     theta_xz(i,j,1)=acos((r_xi(i,j,k+1)*r_zeta(i,j,k+1))/(sqrt(pow2(r_xi(i,j,k+1)))*sqrt(pow2(r_zeta(i,j,k+1))) )  )
     theta_ez(i,j,1)=acos((r_eta(i,j,k+1)*r_zeta(i,j,k+1))/(sqrt(pow2(r_eta(i,j,k+1)))*sqrt(pow2(r_zeta(i,j,k+1))) ) )
  end do
  end do

    !边界间距理想值，分别为两步推进步长 
     dr(:,:,0)=Lk(:,:,k)
     dr(:,:,1)=Lk(:,:,k+1)
     
 !根据边界处约束情况，修正源项
  do j=1,jmax
  do i=1,imax  
     phi_p(i,j,k-1)=phi_p(i,j,k-1) - sigma*atan(theta_rxz-theta_xz(i,j,0))
     phi_q(i,j,k-1)=phi_q(i,j,k-1) - sigma*atan(theta_rez-theta_ez(i,j,0))
     phi_r(i,j,k-1)=phi_r(i,j,k-1) + sigma*atan(dr(i,j,0)-dd(i,j,0))
     phi_p(i,j,k+1)=phi_p(i,j,k+1) + sigma*atan(theta_rxz-theta_xz(i,j,1))
     phi_q(i,j,k+1)=phi_q(i,j,k+1) + sigma*atan(theta_rez-theta_ez(i,j,1))
     phi_r(i,j,k+1)=phi_r(i,j,k+1) - sigma*atan(dr(i,j,1)-dd(i,j,1))
  end do
  end do
 
    
     !中间层源项计算
     phi_p(:,:,k)=phi_p(:,:,k-1)*exp(-a)-phi_p(:,:,k+1)*exp(-b)
     phi_q(:,:,k)=phi_q(:,:,k-1)*exp(-c)-phi_q(:,:,k+1)*exp(-d)
     phi_r(:,:,k)=phi_r(:,:,k-1)*exp(-e)-phi_r(:,:,k+1)*exp(-f)

  
end subroutine

end module


   
        
   
   
    