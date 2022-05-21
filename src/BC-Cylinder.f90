!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module cyl

  USE decomp_2d
  USE variables
  USE param
  IMPLICIT NONE

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: init_cyl, boundary_conditions_cyl, postprocess_cyl, &
            geomcomplex_cyl, visu_cyl, visu_cyl_init
  real(mytype), dimension (4,2) :: body_centroid, body_vel !======added by Xiantao Fan
      
contains

  subroutine geomcomplex_cyl(epsi,nxi,nxf,ny,nyi,nyf,nzi,nzf,dx,yp,remp)

    use decomp_2d, only : mytype
    use param, only : one, two, ten
    use ibm_param
    use dbg_schemes, only: sqrt_prec
    
    implicit none

    integer                    :: nxi,nxf,ny,nyi,nyf,nzi,nzf
    real(mytype),dimension(nxi:nxf,nyi:nyf,nzi:nzf) :: epsi
    real(mytype),dimension(ny) :: yp
    real(mytype)               :: dx
    real(mytype)               :: remp
    integer                    :: i,j,k,p
    real(mytype)               :: xm,ym,r,rads2,kcon
    real(mytype)               :: zeromach
    real(mytype)               :: cexx,ceyy,dist_axi
    real(mytype)               :: k1, k2, k3, k4,L1, L2, L3, L4
    character(len=30) :: filename

   
    zeromach=one
    do while ((one + zeromach / two) .gt. one)
       zeromach = zeromach/two
    end do
    zeromach = ten*zeromach

    ! Intitialise epsi
    epsi(:,:,:)=zero

    ! Update center of moving Cylinder
    !cexx=cex+ubcx*t
    !ceyy=cey+ubcy*t
    ! Update center of moving Cylinder
    !============= added by Xiantao Fan
    ! (idyn 1: Moving Objects-prescribed velocity, 2: Moving Objects-prescribed sine motion, 3: Moving Objects-1-DOF dynamic, 4: Moving Objects-2-DOF dynamic )
    
    !===============
    !===============
    !===============responses computation 
 do p=1,nvol
  if (imove(p).eq.1) then
   if (idyn(p).eq.3.) then
    ! initial
     if (t.eq.0.) then
     body_centroid(p,2)=0
     body_vel(p,2)= 0.01 !provide initial values to avoid NAN      
    endif
    k1=0.5*u1**2*2*ra(p)*zlz* Lift_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_y(p)*body_vel(p,2)-2*pi*fn_y(p)*2*pi*fn_y(p)*body_centroid(p,2)
    k2=0.5*u1**2*2*ra(p)*zlz* Lift_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_y(p)*(body_vel(p,2)+k1*dt/2)-2*pi*fn_y(p)*2*pi*fn_y(p)*(body_centroid(p,2)+dt/2*body_vel(p,2))
    k3=0.5*u1**2*2*ra(p)*zlz* Lift_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_y(p)*(body_vel(p,2)+k2*dt/2)-2*pi*fn_y(p)*2*pi*fn_y(p)*(body_centroid(p,2)+dt/2*body_vel(p,2)+dt*dt/4*k1)       
    k4=0.5*u1**2*2*ra(p)*zlz* Lift_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_y(p)*(body_vel(p,2)+k3*dt)-2*pi*fn_y(p)*2*pi*fn_y(p)*(body_centroid(p,2)+dt*body_vel(p,2)+dt*dt/2*k2)        
    n_body_vel(p,2)=body_vel(p,2)+dt*(k1+2*k2+2*k3+k4)/6
    n_body_centroid(p,2)=body_centroid(p,2)+dt*body_vel(p,2)+dt*dt*(k1+k2+k3)/6    
    
    body_centroid(p,2)=n_body_centroid(p,2)
    body_vel(p,2)=n_body_vel(p,2)
  endif
  
   if (idyn(p).eq.4.) then
    ! initial
    if (t.eq.0.) then
    body_centroid(p,1)=0
    body_vel(p,1) = 0.0 !provide initial values to avoid NAN     
    body_centroid(p,2)=0
    body_vel(p,2) = 0.01 !provide initial values to avoid NAN      
    endif
    !====y-dof:0.5*u1**2*2*ra(p)*zlz* 
    k1=0.5*u1**2*2*ra(p)*zlz* Lift_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_y(p)*body_vel(p,2)-2*pi*fn_y(p)*2*pi*fn_y(p)*body_centroid(p,2)
    k2=0.5*u1**2*2*ra(p)*zlz* Lift_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_y(p)*(body_vel(p,2)+k1*dt/2)-2*pi*fn_y(p)*2*pi*fn_y(p)*(body_centroid(p,2)+dt/2*body_vel(p,2))
    k3=0.5*u1**2*2*ra(p)*zlz* Lift_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_y(p)*(body_vel(p,2)+k2*dt/2)-2*pi*fn_y(p)*2*pi*fn_y(p)*(body_centroid(p,2)+dt/2*body_vel(p,2)+dt*dt/4*k1)       
    k4=0.5*u1**2*2*ra(p)*zlz* Lift_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_y(p)*(body_vel(p,2)+k3*dt)-2*pi*fn_y(p)*2*pi*fn_y(p)*(body_centroid(p,2)+dt*body_vel(p,2)+dt*dt/2*k2)        
    n_body_vel(p,2)=body_vel(p,2)+dt*(k1+2*k2+2*k3+k4)/6
    n_body_centroid(p,2)=body_centroid(p,2)+dt*body_vel(p,2)+dt*dt*(k1+k2+k3)/6    
    
    body_centroid(p,2)=n_body_centroid(p,2)
    body_vel(p,2)=n_body_vel(p,2)
    
    !=======x-dof
    L1=0.5*u1**2*2*ra(p)*zlz* drag_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_x(p)*body_vel(p,1)-2*pi*fn_x(p)*2*pi*fn_x(p)*body_centroid(p,1)
    L2=0.5*u1**2*2*ra(p)*zlz* drag_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_x(p)*(body_vel(p,1)+L1*dt/2)-2*pi*fn_x(p)*2*pi*fn_x(p)*(body_centroid(p,1)+dt/2*body_vel(p,1))
    L3=0.5*u1**2*2*ra(p)*zlz* drag_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_x(p)*(body_vel(p,1)+L2*dt/2)-2*pi*fn_x(p)*2*pi*fn_x(p)*(body_centroid(p,1)+dt/2*body_vel(p,1)+dt*dt/4*L1)       
    L4=0.5*u1**2*2*ra(p)*zlz* drag_tr(p)/(mass(p))-2*epsion(p)*2*pi*fn_x(p)*(body_vel(p,1)+L3*dt)-2*pi*fn_x(p)*2*pi*fn_x(p)*(body_centroid(p,1)+dt*body_vel(p,1)+dt*dt/2*L2)        
    n_body_vel(p,1)=body_vel(p,1)+dt*(L1+2*L2+2*L3+L4)/6
    n_body_centroid(p,1)=body_centroid(p,1)+dt*body_vel(p,1)+dt*dt*(L1+L2+L3)/6    
    
    body_centroid(p,1)=n_body_centroid(p,1)
    body_vel(p,1)=n_body_vel(p,1)
  endif
     
  !=============update coordinate
    
  if (idyn(p).eq.1.) then
      !=============prescribed linear velocity term
     if (t.ne.0.) then
       cexx=cex(p)+ax(p)*(t-ifirst*dt)
       ceyy=cey(p)+ay(p)*(t-ifirst*dt)
     else
       cexx=cex(p)
       ceyy=cey(p)
     endif
       elseif (idyn(p).eq.2.) then
      !=============prescribed sine motion term
     if (t.ne.0.) then
       cexx=cex(p)+ ax(p) * sin (2* pi * fn_x(p) *(t-ifirst*dt))
       ceyy=cey(p)+ ay(p) * sin (2* pi * fn_y(p) *(t-ifirst*dt))
     else
       cexx=cex(p)
       ceyy=cey(p)
     endif

       !======1-DOF
       elseif (idyn(p).eq.3.) then
     if (t.ne.0.) then
       cexx=cex(p)
       ceyy=cey(p)+ n_body_centroid(p,2)
     else
       cexx=cex(p)
       ceyy=cey(p)
     endif
      !======2-DOF
    elseif (idyn(p).eq.4.) then
     if (t.ne.0.) then
       cexx=cex(p)+ n_body_centroid(p,1)
       ceyy=cey(p)+ n_body_centroid(p,2)
     else
       cexx=cex(p)
       ceyy=cey(p)
     endif
 endif
endif

    !
    ! Define adjusted smoothing constant
!    kcon = log((one-0.0001)/0.0001)/(smoopar*0.5*dx) ! 0.0001 is the y-value, smoopar: desired number of affected points 
!
    do k=nzi,nzf
       do j=nyi,nyf
          ym=yp(j)
          do i=nxi,nxf
             xm=real(i-1,mytype)*dx
             r=sqrt_prec((xm-cexx)**two+(ym-ceyy)**two)
             if (r-ra(p).gt.zeromach) then
                cycle
             endif
             epsi(i,j,k)=remp
          enddo
       enddo
    enddo
enddo
return
  end subroutine geomcomplex_cyl

  !********************************************************************
  subroutine boundary_conditions_cyl (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    call inflow (phi)
    call outflow (ux,uy,uz,phi)

    return
  end subroutine boundary_conditions_cyl
  !********************************************************************
  subroutine inflow (phi)

    USE param
    USE variables
    USE decomp_2d
    USE ibm_param

    implicit none

    integer  :: j,k,is
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi

    !call random_number(bxo)
    !call random_number(byo)
    !call random_number(bzo)
    do k=1,xsize(3)
       do j=1,xsize(2)
          bxx1(j,k)=u1+bxo(j,k)*inflow_noise
          bxy1(j,k)=zero+byo(j,k)*inflow_noise
          bxz1(j,k)=zero+bzo(j,k)*inflow_noise
       enddo
    enddo

    if (iscalar.eq.1) then
       do is=1, numscalar
          do k=1,xsize(3)
             do j=1,xsize(2)
                phi(1,j,k,is)=cp(is)
             enddo
          enddo
       enddo
    endif

    return
  end subroutine inflow
  !********************************************************************
  subroutine outflow (ux,uy,uz,phi)

    USE param
    USE variables
    USE decomp_2d
    USE MPI
    USE ibm_param

    implicit none

    integer :: j,k,code
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype) :: udx,udy,udz,uddx,uddy,uddz,cx,uxmin,uxmax,uxmin1,uxmax1

    udx=one/dx; udy=one/dy; udz=one/dz; uddx=half/dx; uddy=half/dy; uddz=half/dz

    uxmax=-1609._mytype
    uxmin=1609._mytype
    do k=1,xsize(3)
       do j=1,xsize(2)
          if (ux(nx-1,j,k).gt.uxmax) uxmax=ux(nx-1,j,k)
          if (ux(nx-1,j,k).lt.uxmin) uxmin=ux(nx-1,j,k)
       enddo
    enddo

    call MPI_ALLREDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,MPI_COMM_WORLD,code)

    if (u1 == zero) then
       cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
    elseif (u1 == one) then
       cx=uxmax1*gdt(itr)*udx
    elseif (u1 == two) then
       cx=u2*gdt(itr)*udx    !works better
    else
       cx=(half*(u1+u2))*gdt(itr)*udx
    endif

    do k=1,xsize(3)
       do j=1,xsize(2)
          bxxn(j,k)=ux(nx,j,k)-cx*(ux(nx,j,k)-ux(nx-1,j,k))
          bxyn(j,k)=uy(nx,j,k)-cx*(uy(nx,j,k)-uy(nx-1,j,k))
          bxzn(j,k)=uz(nx,j,k)-cx*(uz(nx,j,k)-uz(nx-1,j,k))
       enddo
    enddo

    if (iscalar==1) then
       if (u2==zero) then
          cx=(half*(uxmax1+uxmin1))*gdt(itr)*udx
       elseif (u2==one) then
          cx=uxmax1*gdt(itr)*udx
       elseif (u2==two) then
          cx=u2*gdt(itr)*udx    !works better
       else
          stop
       endif

       do k=1,xsize(3)
          do j=1,xsize(2)
             phi(nx,j,k,:)=phi(nx,j,k,:)-cx*(phi(nx,j,k,:)-phi(nx-1,j,k,:))
          enddo
       enddo
    endif

    if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) &
       write(*,*) "Outflow velocity ux nx=n min max=",real(uxmin1,4),real(uxmax1,4)

    return
  end subroutine outflow
  !********************************************************************
  subroutine init_cyl (ux1,uy1,uz1,phi1)

    USE decomp_2d
    USE decomp_2d_io
    USE variables
    USE param
    USE MPI
    use dbg_schemes, only: exp_prec

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,um
    integer :: k,j,i,ii,is,code

    if (iscalar==1) then

       phi1(:,:,:,:) = zero !change as much as you want

    endif

    ux1=zero; uy1=zero; uz1=zero

    if (iin.ne.0) then
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)

       do k=1,xsize(3)
          do j=1,xsize(2)
             do i=1,xsize(1)
                ux1(i,j,k)=init_noise*(ux1(i,j,k)-0.5)
                uy1(i,j,k)=init_noise*(uy1(i,j,k)-0.5)
                uz1(i,j,k)=init_noise*(uz1(i,j,k)-0.5)
             enddo
          enddo
       enddo

       !modulation of the random noise
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-yly/2.
             if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
             um=exp_prec(-zptwo*y*y)
             do i=1,xsize(1)
                ux1(i,j,k)=um*ux1(i,j,k)
                uy1(i,j,k)=um*uy1(i,j,k)
                uz1(i,j,k)=um*uz1(i,j,k)
             enddo
          enddo
       enddo
    endif

    !INIT FOR G AND U=MEAN FLOW + NOISE
    do k=1,xsize(3)
       do j=1,xsize(2)
          do i=1,xsize(1)
             ux1(i,j,k)=ux1(i,j,k)+u1
             uy1(i,j,k)=uy1(i,j,k)
             uz1(i,j,k)=uz1(i,j,k)
          enddo
       enddo
    enddo

#ifdef DEBG
    if (nrank .eq. 0) write(*,*) '# init end ok'
#endif

    return
  end subroutine init_cyl
  !********************************************************************

  !############################################################################
  subroutine postprocess_cyl(ux1,uy1,uz1,ep1)

    USE MPI
    USE decomp_2d
    USE decomp_2d_io
    USE var, only : uvisu
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    USE ibm_param
    use dbg_schemes, only: sqrt_prec
    
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1

  end subroutine postprocess_cyl

  subroutine visu_cyl_init (visu_initialised)

    use decomp_2d, only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu, only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "vort", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_cyl_init
  !############################################################################
  !!
  !!  SUBROUTINE: visu_cyl
  !!      AUTHOR: FS
  !! DESCRIPTION: Performs cylinder-specific visualization
  !!
  !############################################################################
  subroutine visu_cyl(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var, only : ux2, uy2, uz2, ux3, uy3, uz3
    USE var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    USE var, only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var, ONLY : nxmsize, nymsize, nzmsize
    use visu, only : write_field
    use ibm_param
    

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    character(len=32), intent(in) :: num
    integer :: p


    ! Write vorticity as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    !x-derivatives
    !=====added by Xiantao Fan 
    !==========recompute the fluid velocity, assume the fluid velocity = the moving body velocity
!=======comput the IBM velocity based on different case------important parameter 
!==========recompute the fluid velocity, assume the fluid velocity = the moving body velocity
do p=1,nvol
if (imove(p).eq.1) then
  if (idyn(p).eq.1.) then
         ubcx=ax(p)
         ubcy=ay(p)
    elseif (idyn(p).eq.2.) then
         ubcx = 2* pi * fn_x(p) * ax(p) * cos (2* pi * fn_x(p) *(t-ifirst*dt))
         ubcy = 2* pi * fn_y(p) * ay(p) * cos (2* pi * fn_y(p) *(t-ifirst*dt))
    elseif (idyn(p).eq.3.) then 
         ubcx = 0
         ubcy = n_body_vel(p,2)         
    elseif (idyn(p).eq.4.) then 
         ubcx = n_body_vel(p,1)
         ubcy = n_body_vel(p,2)
    endif
endif    
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    !y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    !!z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
    !!all back to x-pencils
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)
    call transpose_y_to_x(td2,tg1)
    call transpose_y_to_x(te2,th1)
    call transpose_y_to_x(tf2,ti1)
    call transpose_y_to_x(ta2,td1)
    call transpose_y_to_x(tb2,te1)
    call transpose_y_to_x(tc2,tf1)
enddo
    !du/dx=ta1 du/dy=td1 and du/dz=tg1
    !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
    !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
    !VORTICITY FIELD
    di1 = zero
    di1(:,:,:)=sqrt(  (tf1(:,:,:)-th1(:,:,:))**2 &
                    + (tg1(:,:,:)-tc1(:,:,:))**2 &
                    + (tb1(:,:,:)-td1(:,:,:))**2)
    call write_field(di1, ".", "vort", trim(num), flush = .true.) ! Reusing temporary array, force flush

    !Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,: ) = - half*(ta1(:,:,:)**2+te1(:,:,:)**2+ti1(:,:,:)**2) &
                  - td1(:,:,:)*tb1(:,:,:) &
                  - tg1(:,:,:)*tc1(:,:,:) &
                  - th1(:,:,:)*tf1(:,:,:)
    call write_field(di1, ".", "critq", trim(num), flush = .true.) ! Reusing temporary array, force flush

  end subroutine visu_cyl

end module cyl

