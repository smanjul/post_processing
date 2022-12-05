!                                                                      
! This code reads a particular azimuthal mode generated in each z-plane from a matlab code
! and writes a tecplot file to plot iso-surfaces of that mode.
!
      program modes            
      include 'parade.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1,m3md=m3m+2,m3mp=m3m+4) 
      character*4 dummy
      common/dim/n1,n1m,n2,n2m,n3,n3m 
      common/d13/alx3
      common/vperin/vper,epsil,lamb   
      common/njump/n1p,n2p,n3p       
      common/tstep/dt,beta,re       
      common/inslnsw/inslws,inslwn,inslwr
      common/parcoo/r0           
      common/section/nsect                     
      common /ametrz/  g3rc(m3),g3rm(m3),str3
       common/parava/ray,pra,pec,ren
      open(15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) n1,n2,n3,nsst,nwrit,nread
        read(15,301) dummy
        read(15,*) n1p,n2p,n3p
        read(15,301) dummy
        read(15,*) ntst,tprint,tpin,tmax,ireset
        read(15,301) dummy
        read(15,*) alx3,istr3,str3,rmed31,etdp3,strb3
        read(15,301) dummy
        read(15,*)strr,r0,istr,rmed,etdp,strb
        read(15,301) dummy
        read(15,*) ray,ri,ros,pra,dt,resid,cflmax
        read(15,301) dummy
        read(15,*) inslws,inslwn,inslwr,ifugo,icorio,ibuo,isca
        read(15,301) dummy
        read(15,*) epsil,lamb,rper,nsect
        read(15,301) dummy
        read(15,*) istat
        read(15,301) dummy
        read(15,*) nlev,dism,nn
        read(15,301) dummy
        read(15,*) idtv,dtmax,cfllim
        read(15,301) dummy
        read(15,*) nini,nfin,nstri
        read(15,301) dummy
        read(15,*) nson
301     format(a4)                
      close(15)
c       

!       ren=sqrt(ray/pra)
!       pec=sqrt(ray*pra)
       open(333,file='input.in')
         read(333,301)dummy
         read(333,*)ren,pec
       close(333)
c
       print*,1,pec,ren

      n1m=n1-1    
      if(n1.eq.1) n1m=1
      n2m=n2-1   
      n3m=n3-1  
      call gcurv 
      stop      
      end   
c
c************************
c***********************************************************************
      subroutine gcurv    
      include 'parade.f'
      dimension q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)   
      dimension dens(m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
      common/dim/n1,n1m,n2,n2m,n3,n3m       
      common/d1/re,tfin,eps                
      common/tstep/dt,beta,ren            
      common/njump/n1p,n2p,n3p           
      common/d13/alx3                   
      common/corrt/thetac(m1),thetam(m1)
      common/cordvo/rc(m2),rm(m2),zz(m3)       
      common /metrr/  g2rm(m2),g2rc(m2)
      common /ametrz/  g3rc(m3),g3rm(m3),str3
      character*5 nfil
      call meshes(alx3)              
       call cordin_read
      call indic
      write(6,754)n1,n2,n3            
  754 format(/,4x,'Grid :',2x,i3,'x',i3,'x',i3,//)  
      call avgtime
c      call rmstime
      return
      end  
c         
c*********************************************************************************************************************
c***********************************************************************
      subroutine indic  
      include 'parade.f'
c     implicit real*4(a-h,o-z)
      common/dim/n1,n1m,n2,n2m,n3,n3m  
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/inwal/jpc(m2),jup(m2),jmc(m2),jum(m2) 
      n1mm=n1m-1     
      do 1 ic=1,n1m 
        ip=ic+1    
        imv(ic)=ic-1 
        if(ic.eq.1) imv(ic)=n1m   
        ipv(ic)=ic+1             
        if(ic.eq.n1m) ipv(ic)=1 
    1 continue                 
      do 4 kc=1,n3m           
        kmv(kc)=kc-1         
        kpv(kc)=kc+1        
        if(kc.eq.1) kmv(kc)=kc   
        if(kc.eq.n3m) kpv(kc)=kc
    4 continue 
      do 3 jc=1,n2m    
        jp=jc+1       
        jmv(jc)=jc-1 
        jpv(jc)=jc+1
        if(jc.eq.1) jmv(jc)=jc    
        if(jc.eq.n2m) jpv(jc)=jc 
    3 continue                  
      do 15 jc=1,n2m           
        jpc(jc)=jpv(jc)-jc    
        jmc(jc)=jc-jmv(jc)   
        jup(jc)=1-jpc(jc)   
        jum(jc)=1-jmc(jc)  
   15 continue            
      return             
      end               
c                      
c***********************************************************************
      subroutine meshes(alx3)   
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
      common/dim/n1,n1m,n2,n2m,n3,n3m       
      common/vperin/vper,epsil,lamb   
      common/parcoo/r0                     
      common/section/nsect                     
      pi=2.d0*asin(1.d0)                  
      if(nsect.ne.1) then
        dx1=2.d0*pi/float(n1m)
      else
        dx1=2.d0*pi/(float(n1m)*float(lamb))
      endif
      dx2=1./float(n2m)                 
      dx3=1./float(n3m)       
      write(6,99) dx1,dx2,dx3 
 99   format(/,2x,'dx1,dx2,dx3 = ',3f10.4)
      dx1=1.d0/dx1    
      dx2=1.d0/dx2   
      dx3=1.d0/dx3  
      dx1q=dx1*dx1                           
      dx2q=dx2*dx2                          
      dx3q=dx3*dx3                         
      return                              
      end                                
c 
c -------------------------------------------------------------------
C       read coordinates
      subroutine cordin_read 
      include 'parade.f'

       dimension zm(m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
      common/dim/n1,n1m,n2,n2m,n3,n3m       
      common/corrt/thetac(m1),thetam(m1)
      common/cordvo/rc(m2),rm(m2),zz(m3)   
      common /metrr/  g2rm(m2),g2rc(m2)
      common /ametrz/  g3rc(m3),g3rm(m3),str3
       
cc
c     AZIMUTHAL COORDINATE DEFINITION
c
      do i=1,n1
        thetac(i)= float(i-1)/dx1
      end do
      do i=1,n1m
        thetam(i)= (float(i-1)+0.5d0)/dx1
      end do


c     WRITE GRID INFORMATION
c
      open(98,file='radcor.out',status='old')
      do j=1,n2
        read(98,345) nj,rc(j),rm(j),g2rc(j),g2rm(j)
      end do
      close(98)
      open(78,file='axicor.out',status='old')
      do k=1,n3
        read(78,345) nk,zz(k),zm(k),g3rc(k),g3rm(k)
      end do
      close(78)
 345  format(i4,4(2x,e13.5))


      return
       end
C--------------------------------------------------------------------
      subroutine avgtime
      include 'parade.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m 
      dimension temp(m1,m2,m3),disste(m1,m2,m3)
      dimension a1(m1,m2,m3),a3(m1,m2,m3)
      dimension a2(m1,m2,m3),diss(m1,m2,m3),diss2(m2,m3)
      dimension vis(m1,m2,m3),cp(m1,m2,m3)
      dimension dmass(m1,m2,m3),condt(m1,m2,m3)
      dimension expan(m1,m2,m3)
       dimension tempme(m1,m2,m3),dissteme(m1,m2,m3)
      dimension a1me(m1,m2,m3),a3me(m1,m2,m3)
      dimension a2me(m1,m2,m3),dissme(m1,m2,m3)
      dimension a1pert(m1,m2,m3),a3pert(m1,m2,m3)
      dimension a2pert(m1,m2,m3),tempert(m1,m2,m3)
      dimension upert(m1,m2,m3),vpert(m1,m2,m3)
      double precision yp1(m1,m2),yp2(m1,m2),mod_phys(m1,m2)
      integer it
      real pi
      dimension u1(m1,m2,m3),u2(m1,m2,m3)  
      dimension u3(m1,m2,m3), dens(m1,m2,m3)
      dimension voc1(m1,m2,m3),voc2(m1,m2,m3)
      dimension qx(m1,m2,m3),qy(m1,m2,m3)
      dimension densc(m1,m2,m3),voc3(m1,m2,m3)
!       dimension yp1(m1,m2),yp2(m1,m2),rs(m2)
!      dimension qx(m1,m2,m3),qy(m1,m2,m3)
      dimension qx_sec(m1,m2,m3),qy_sec(m1,m2,m3)
 
        character*70 namfi3
        character*70 namfi4
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
      common/corrt/thetac(m1),thetam(m1)
      common/cordvo/rc(m2),rm(m2),zz(m3)   
      common /metrr/  g2rm(m2),g2rc(m2)
      common /ametrz/  g3rc(m3),g3rm(m3),str3
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/parcoo/r0           
      common/d13/alx3
      common/filein/it,itini,incre
       common/asym/isym(m1)
      character*5 ipfi
 

         pi=2.d0*asin(1.)
      if(n1m.ne.1) then
      do i=1,n1m
       isym(i) = i + n1m/2
       if(isym(i).gt.n1m) isym(i) = isym(i) - n1m
      enddo
      end if

!  Read the matlab file
       open(1,file='mode10_re5000a2_5.dat')
        do i=1,n1
          read(1,*)(mod_phys(i,j),j=1,n2)
        enddo

! Transform the grid
! No need to transform the mod_phys as it was already transformed
      do 1 i=1,n1
        do 1 j=1,n2
          yp1(i,j)=rc(j)*cos(thetac(i))
          yp2(i,j)=rc(j)*sin(thetac(i))
    1 continue 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      open(11,file='mode10_re5000a2_5_tec.dat',status='unknown')
      write(11,*) 'VARIABLES = "X","Y","mod"'
      write(11,*) 'ZONE I=',n1,', J=',n2,', F=BLOCK'

c     coordinate x
      write(11,234) ((yp1(i,j),i=1,n1),j=1,n2)
c     coordinate y
      write(11,234) ((yp2(i,j),i=1,n1),j=1,n2)
c     w azimuthal mode
      write(11,234) ((mod_phys(i,j),i=1,n1),j=1,n2)
      close(11)
       write(*,*) 'Tecplot file is written'
234     format(1000(f16.9,2x))


      return
      end
