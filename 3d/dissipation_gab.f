cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c          Calcolo delle medie partendo dai file boum***.dat           c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
      program medie              
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
      dimension yp1(m1,m2),yp2(m1,m2)
      integer it
      real pi
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
 

      if(n1m.ne.1) then
      do i=1,n1m
       isym(i) = i + n1m/2
       if(isym(i).gt.n1m) isym(i) = isym(i) - n1m
      enddo
      end if
c----inizializzazione
       do k=1,n3
              do j=1,n2
                     do i=1,n1
                     tempme(i,j,k)=0.0
                     dissteme(i,j,k)=0.0
                     dissme(i,j,k)=0.0
                     a1me(i,j,k)=0.0
                     a2me(i,j,k)=0.0
                     a3me(i,j,k)=0.0
                     enddo
              enddo
       enddo

c-----lettura dati input 
         write(*,*)'n3m, n3',n3m,n3      
      write(6,*) ' no of time units, starting time, increment'
      read(5,*) it,itini,incre    
c       it=20
c       itini=300
c       incre=5    


ccc ---inizio ciclo esterno (media nel tempo)

      itime=itini

      do iti=1,it
!      iftime=10.*itime+0.5
      iftime=1.*itime+0.5
      write(ipfi,99) iftime
   99 format(i5.5)
c   99 format(i6.6)
       namfi3='dati/boum'//ipfi//'.dat'

      write(6,201) namfi3
  201 format(10x,'files taken ',a70)
c-----lettura file boum***.dat
      open(62,file=namfi3,form='unformatted')
!      open(62,file=namfi3)
      read(62) n1,n2,n3               
      read(62) aaa,nu,aaa,dtime
      read(62)(((a1(i,j,k),i=1,n1),j=1,n2),k=1,n3)
      read(62)(((a2(i,j,k),i=1,n1),j=1,n2),k=1,n3)
      read(62)(((a3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
      close(62)

	write(*,*)'file count= ',iti
      
c-----ricostruzione q2------------------------
!      do k=1,n3m
!        usrnu3=dx3/g3rm(k)
!        kp=k+1
!        do i=1,n1m
!          ip=ipv(i)
!          a2(i,1,k) = 0.
!           a2(i,n2,k) = 0.
!          do j=2,n2m
!            jm=j-1
!            udx2=dx2/g2rm(jm)
!            usurm = 1./rm(jm)
!            a2(i,j,k)=a2(i,jm,k)-(
!     %             (a1(ip,jm,k)-a1(i,jm,k))*dx1*usurm
!     %            +(a3(i,jm,kp)-a3(i,jm,k))*usrnu3
!     %                                       )/(udx2*usurm)
!          end do
!        end do
!      enddo
c-----ricostruzione di u2---------------------------------
!       do k=1,n3m
!              do j=2,n2
!                     do i=1,n1m
!                     a2(i,j,k)=a2(i,j,k)/rc(j)
!                     enddo
!              enddo
!       enddo
c      ricostruzione di u2 sull'asse'
       do k=1,n3m
              do i=1,n1m
            a2(i,1,k)=(a2(i,2,k)-a2(isym(i),2,k))*0.5d0
              enddo
       enddo
       

c-----u1 in n1------------------------------------------
       do k=1,n3m
              do j=1,n2m
              a1(n1,j,k)=a1(1,j,k)
              enddo
       enddo       
c----------------------------------------------------------------------------
       
!       call dissiptemp(disste,condt,temp)
       call dissipation(a1,a2,a3,vis,dmass,diss)

       do k=1,n3
              do j=1,n2
                     do i=1,n1
!                     tempme(i,j,k)=tempme(i,j,k)+temp(i,j,k)
!                     dissteme(i,j,k)=dissteme(i,j,k)+disste(i,j,k)
!                     dissme(i,j,k)=dissme(i,j,k)+diss(i,j,k)
                     a1me(i,j,k)=a1me(i,j,k)+a1(i,j,k)
                     a2me(i,j,k)=a2me(i,j,k)+a2(i,j,k)
                     a3me(i,j,k)=a3me(i,j,k)+a3(i,j,k)
                     enddo
              enddo
       enddo
      

c       do k=1,n3
c	dissteme_k(k)=0.d0
c	dissme_k(k)=0.d0
c              do j=1,n2
c                     do i=1,n1
c                     dissteme_k(k)=dissteme_k(k)+disste(i,j,k)
c                     dissme_k(k)=dissme_k(k)+diss(i,j,k)
c                     enddo
c              enddo
c       enddo
c       do k=1,n3
c              do j=1,n2
c	dissteme_rk(j,k)=0.d0
c	dissme_rk(j,k)=0.d0
c                     do i=1,n1
c                     dissteme_rk(j,k)=dissteme_rk(j,k)+disste(i,j,k)
c                     dissme_rk(j,k)=dissme_rk(j,k)+diss(i,j,k)
c                     enddo
c              enddo
c       enddo
c-----fine ciclo esterno
       itime=itime+incre 
       enddo
       uit=1.0/float(it)
       print*,'indici', it, uit
c      ==========time average=====
       do k=1,n3
              do j=1,n2
                     do i=1,n1
!                     tempme(i,j,k)=tempme(i,j,k)*uit
!                     dissteme(i,j,k)=dissteme(i,j,k)*uit
!                     dissme(i,j,k)=dissme(i,j,k)*uit
                     a1me(i,j,k)=a1me(i,j,k)*uit
                     a2me(i,j,k)=a2me(i,j,k)*uit
                     a3me(i,j,k)=a3me(i,j,k)*uit
                     enddo
              enddo
       enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   Calculations for perturbation
       do k=1,n3
         do j=1,n2
           do i=1,n1
!             tempert(i,j,k)=tempme(i,j,k)-temp(i,j,k)   !temperature perturbation
             a1pert(i,j,k)=a1me(i,j,k)-a1(i,j,k)      !azimuthal velocity perturbation
             a2pert(i,j,k)=a2me(i,j,k)-a2(i,j,k)      !radial velocity perturbation
             a3pert(i,j,k)=a3me(i,j,k)-a3(i,j,k)      !axial velocity perturbation
           enddo
         enddo
       enddo

        do kc=1,n3m
          do jc=1,n2
            a1pert(n1,jc,kc)=a1pert(1,jc,kc)
            a2pert(n1,jc,kc)=a2pert(1,jc,kc)
            a3pert(n1,jc,kc)=a3pert(1,jc,kc)
!            tempert(n1,jc,kc)=tempert(1,jc,kc)
          enddo
        enddo

      do 1 i=1,n1
        do 1 j=1,n2
          yp1(i,j)=rc(j)*cos(thetac(i))
          yp2(i,j)=rc(j)*sin(thetac(i))
    1 continue

!   transformation of velocity from polar to cartesian
      do k=1,n3
        do j=1,n2
          do i=1,n1
!   transformatino of perturbation
            upert(i,j,k)=a2pert(i,j,k)*cos(thetac(i))
     %                  -a1pert(i,j,k)*sin(thetac(i))
            vpert(i,j,k)=a2pert(i,j,k)*sin(thetac(i))
     %                  +a1pert(i,j,k)*cos(thetac(i))
          enddo
        enddo
      enddo
      open(33,file='pertfield_untransformed.dat')
      write(33,*) 'VARIABLES = "X","Y","Z","UP","VP","WP","TP"'
      write(33,*) 'ZONE I=',n1,', J=',n2,', K=',n3,', F=BLOCK'

c     coordinata x
      write(33,235) (((yp1(i,j),i=1,n1),j=1,n2),k=1,n3)
c     coordinata y
      write(33,235) (((yp2(i,j),i=1,n1),j=1,n2),k=1,n3)
c     coordinata z
      write(33,235) (((zz(k),i=1,n1),j=1,n2),k=1,n3)
c     U perturbation
      write(33,235) (((a1pert(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     V perturbation
      write(33,235) (((a2pert(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     W perturbation
      write(33,235) (((a3pert(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     T perturbation
!      write(33,235) (((tempert(i,j,k),i=1,n1),j=1,n2),k=1,n3)

      close(33)
235     format(1000(f16.9,2x))



!       call totaverage(tempme,dissteme,dissme) !volume average
!       call totaverage_k(dissteme,dissme) !r,theta average
!       call totaverage_rk(dissteme,dissme) !theta average
!       call tecplot_diss(dissteme,dissme) !write dissipations fields
       call tecplot_boum(a1me,a2me,a3me,tempme) !write mean fields
       return   
      end   
c-----------------------------------------------------------------------------
       subroutine dissipation(u1,u2,u3,vis,dmass,diss)
       include'parade.f'

c                       2       
c       dissipation:  ----  nu E:E  
c                      Re
c      La dissipazione e` calcolata a centro cella su tutta la velocita`

       common/dim/n1,n1m,n2,n2m,n3,n3m 
       common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
      common/corrt/thetac(m1),thetam(m1)
      common/cordvo/rc(m2),rm(m2),zz(m3)   
      common /metrr/  g2rm(m2),g2rc(m2)
      common /ametrz/  g3rc(m3),g3rm(m3),str3
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/parcoo/r0           
      common/d13/alx3
       common/parava/ray,pra,pec,ren

       dimension u1(m1,m2,m3),u3(m1,m2,m3)
      dimension u2(m1,m2,m3),diss(m1,m2,m3)
      dimension vis(m1,m2,m3),cp(m1,m2,m3)
      dimension dmass(m1,m2,m3),condt(m1,m2,m3)
      dimension expan(m1,m2,m3)
       dimension up1(m1,m2,m3),up2(m1,m2,m3)
       
       do k=1,n3m
              kp=k+1
              udzm=dx3/g3rm(k)
              do j=1,n2m
                     jp=j+1
                     udrm=dx2/g2rm(j)
                     urm=1/rm(j)
                     udth=dx1*urm
                     do i=1,n1m
                     ip=i+1                            


c                          du_r
c                   E_rr=------
c                            dr       
                            
                     dis22=(u2(i,jp,k)-u2(i,j,k))*udrm

c                          du_z
c                   E_zz=------
c                            dz       
                     dis33=(u3(i,j,kp)-u3(i,j,k))*udzm

c                            1    du_th      u_r
c                   E_tt= ---  ------- +  ----
c                             r    dth                r       
                     dis11=(u1(ip,j,k)-u1(i,j,k))*udth+
     %                  (u2(i,jp,k)+u2(i,j,k))*0.5*urm

                     diss(i,j,k)=(dis11*dis11+dis22*dis22+dis33*dis33)
                     enddo
              enddo
       enddo
       
       
c                            1    du_r      du_z
c                   E_rz= --- (----- +   ------)
c                             2    dz               dr       

       call staggu2in3(u2,up1) !up1=u_r
       call staggu3in2(u3,up2) !up2=u_z
       do k=1,n3m
              kp=k+1
              udzm=dx3/g3rm(k)
              udzc=dx3/g3rc(k)
              do j=1,n2m
                     jp=j+1
                     udrm=dx2/g2rm(j)
                     udrc=dx2/g2rc(j)
                     do i=1,n1m       
                     dis23d1=(up1(i,j,kp)-up1(i,j,k))*udzm
                     dis23d2=(up2(i,jp,k)-up2(i,j,k))*udrm
                     dis23=0.5*(dis23d1+dis23d2)
                     diss(i,j,k)=diss(i,j,k)+2.0*(dis23*dis23)
                     enddo
              enddo
       enddo

c                            1     du_th   1   du_z
c                   E_thz= --- (----- +  -- ------)
c                             2     dz             r   dth       

       call staggu1in3(u1,up1) !up1=u_th
       call staggu3in1(u3,up2) !up2=u_z
       do k=1,n3m
              kp=k+1
              udzm=dx3/g3rm(k)
              udzc=dx3/g3rc(k)
              do j=1,n2m
                     urm=1/rm(j)
                     udth=dx1*urm
                     do i=1,n1m
                     ip=i+1
                        dis13d1=(up1(i,j,kp)-up1(i,j,k))*udzm
                        dis13d2=(up2(ip,j,k)-up2(i,j,k))*udth
                        dis13=0.5*(dis13d1+dis13d2)
                        diss(i,j,k)=diss(i,j,k)+2.0*(dis13*dis13)                            
                     enddo
              enddo
       enddo 

c                            1     du_th     u_th       1  du_r
c                   E_thr= --- (-----  -   ----  +   -- ------)
c                             2     dr               r         r  dth       
       
       call staggu1in2(u1,up1) !up1=u_th
       call staggu2in1(u2,up2) !up2=u_r
       do k=1,n3m
              do j=1,n2m
                     jp=j+1
                     udrm=dx2/g2rm(j)
                     udrc=dx2/g2rc(j)
                     urm=1/rm(j)
                     urc=1/rm(j)
                     udth=dx1*urm
                     udthc=dx1*urc
                     do i=1,n1m
                     ip=i+1                            
                        dis12d1=(up1(i,jp,k)-up1(i,j,k))*udrm-
     %                  (u1(ip,j,k)+u1(i,j,k))*0.5*urm
                        dis12d2=(up2(ip,j,k)-up2(i,j,k))*udth
                        dis12=0.5*(dis12d1+dis12d2)
                        diss(i,j,k)=diss(i,j,k)+2.0*(dis12*dis12)
                     enddo
              enddo
       enddo 
       dsren=2.0/ren
       do k=1,n3m
              do j=1,n2m
                     do i=1,n1m       
                 diss(i,j,k)=vis(i,j,k)*diss(i,j,k)*dsren/dmass(i,j,k)
                     enddo
              enddo
       enddo 
       return   
      end   
c***********************************************************************

      subroutine dissiptemp(disste,condt,dens) 
      include'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
      common/corrt/thetac(m1),thetam(m1)
      common/cordvo/rc(m2),rm(m2),zz(m3)   
      common /metrr/  g2rm(m2),g2rc(m2)
      common /ametrz/  g3rc(m3),g3rm(m3),str3
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/parcoo/r0           
      common/d13/alx3
       common/parava/ray,pra,pec,ren
      dimension dens(m1,m2,m3),disste(m1,m2,m3)
      dimension gg2(m2),gg3(m3)
      dimension vis(m1,m2,m3),cp(m1,m2,m3)
      dimension dmass(m1,m2,m3),condt(m1,m2,m3)
      dimension expan(m1,m2,m3)

c       write(6,*) 'dissiptemp '
c
c
c                                1  |                           |  1
c       temperature dissipation ---<| nabla (lambda T) *nabla  T|---- >
c                               Pe  |                           |rho*cp
c
c       T can be either the whole velocity filed or the fluctuating
c       field depending on the value of the flag IMED (0 or 1 respectively)
c
        usdx1 = 1./dx1
        gg2(1)=1./(g2rc(2))
        do j=2,n2m-1
           gg2(j)=1./(g2rc(j+1)+g2rc(j))
        enddo
        gg2(n2m)=1./(g2rc(n2m))

        gg3(1)=1./(g3rc(2))
        do k=2,n3m-1
           gg3(k)=1./(g3rc(k+1)+g3rc(k))
        enddo
        gg3(n3m)=1./(g3rc(n3m))
!$OMP  PARALLEL DO
!$OMP$ DEFAULT(none)
!$OMP$ SHARED(n1m,n2m,n3m,dens)
!$OMP$ SHARED(imv,ipv,jmv,jpv,kmv,kpv)
!$OMP$ SHARED(disste)
!$OMP$ SHARED(rm,dx1,dx2,dx3)
!$OMP$ SHARED(pec)
!$OMP$ FIRSTPRIVATE(usdx1,vol,gg3,gg2)
!$OMP$ PRIVATE(famet)
!$OMP$ PRIVATE(i,j,k,im,ip,jm,jp,jmm,jpp,km,kp,kmm,kpp)
!$OMP$ PRIVATE(udx1,udx2,udx3,h31,h32,h33)
!$OMP$ PRIVATE(dissip2)
        do k=1,n3m
          kpp=kpv(k)
          kmm=kmv(k)
          udx3=dx3*gg3(k)
CG          udx3=dx3/(g3rm(k))*0.5
          do j=1,n2m
          jmm=jmv(j)
          jpp=jpv(j)
          udx1=dx1/rm(j)*0.5
          udx2=dx2*gg2(j)
CG          udx2=dx2/(g2rm(j))*0.5
            do i=1,n1m
              im= imv(i)
              ip= ipv(i)
c
      h31c=((dens(ip,j,k)*condt(ip,j,k))-
     %     (dens(im,j,k)*condt(im,j,k))
     %    )*udx1
c
      h32c=((dens(i,jpp,k)*condt(i,jpp,k))-
     %     (dens(i,jmm,k)*condt(i,jmm,k))
     %    )*udx2
c
      h33c=((dens(i,j,kpp)*condt(i,j,kpp))-
     %     (dens(i,j,kmm)*condt(i,j,kmm))
     %    )*udx3
c       
      h31=((dens(ip,j,k))-
     %     (dens(im,j,k))
     %    )*udx1
c
      h32=((dens(i,jpp,k))-
     %     (dens(i,jmm,k))
     %    )*udx2
c
      h33=((dens(i,j,kpp))-
     %     (dens(i,j,kmm))
     %    )*udx3
c       
       disste(i,j,k)=(h31*h31c + h32*h32c + h33*h33c)/ pec
c	write(*,*)h32,condt(i,j,k),i,j,k
             end do
           end do
         end do
!$OMP  END PARALLEL DO


      return   
      end     

c----------------------------------------------------------------------------------
       subroutine totaverage(tempme,dissteme,dissme)
       include 'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       dimension tempme(m1,m2,m3),dissteme(m1,m2,m3)
       dimension dissme(m1,m2,m3)
       common/parava/ray,pra,pec,ren
       common /metrr/  g2rm(m2),g2rc(m2)
      common /ametrz/  g3rc(m3),g3rm(m3),str3
       common/cordvo/rc(m2),rm(m2),zz(m3)  
       etatemp=0.0
       avertem=0.0
       eta=0.0
       vol=0.0
       do k=1,n3m
              do j=1,n2m
                     fac2=rm(j)*g2rm(j)*g3rm(k)
                     do i=1,n1m
                     etatemp=etatemp+dissteme(i,j,k)*fac2
                     eta=eta+dissme(i,j,k)*fac2
                     avertem=avertem+tempme(i,j,k)*fac2
                     vol=vol+fac2
                     enddo
              enddo
       enddo
       print*,'vol',vol
       print*,'pec',pec
       eta=eta/vol
       etatemp=etatemp/vol
       avertem=avertem/vol
       etnusse=etatemp*pec
       enusse=eta*pec+1.0       
       print*, 'mean temperature',avertem
       print*, 'enegy dissipation rate:',eta,'Nusselt',enusse
       print*, 'temp dissipation rate:',etatemp,'Nusselt',etnusse
	open(42,file='dissipation.out')
	write(42,*)avertem,eta,etatemp,enusse,etnusse
	close(42)
c-----

       return
       end
c-----------------------------------------------------------------
       subroutine totaverage_k(dissteme,dissme)
       include 'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       dimension dissteme(m1,m2,m3),disst_k(m3)
       dimension dissme(m1,m2,m3),diss_k(m3)
       dimension disstec_k(m3),dissc_k(m3)
       common/parava/ray,pra,pec,ren
       common /metrr/  g2rm(m2),g2rc(m2)
      common /ametrz/  g3rc(m3),g3rm(m3),str3
       common/cordvo/rc(m2),rm(m2),zz(m3)  
       common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
       common/parcoo/r0           
       dimension yp1(n1,n2),yp2(n1,n2)
       common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/corrt/thetac(m1),thetam(m1)
	pi = 2.*asin(1.)
       vol = 1./( pi*r0**2*dx1*dx2)
       do k=1,n3m
        etatemp=0.0
        eta=0.0
       vol1=0.0
              do j=1,n2m
                     fac2=rm(j)*g2rm(j)
                     do i=1,n1m
                     etatemp=etatemp+(dissteme(i,j,k)*fac2)
                     eta=eta+(dissme(i,j,k)*fac2)
                     vol1=vol1+fac2
                     enddo
              enddo
	     diss_k(k)=eta*vol
	     disst_k(k)=etatemp*vol
c            print*,'vol',vol,1./vol1
       enddo
c            print*,'vol',vol,1./vol1
	open(42,file='dissipation_k.out')
	do k=1,n3m
	write(42,*)zz(k),diss_k(k),disst_k(k)
	enddo
	close(42)

      do 4 kc=1,n3m 
        km=kmv(kc)                                                    
            dissc_k(kc)=(diss_k(kc)+diss_k(km))*0.5
            disstec_k(kc)=(disst_k(kc)+disst_k(km))*0.5
    4 continue
          dissc_k(1)=dissc_k(2)
          dissc_k(n3)=dissc_k(n3m)
           disstec_k(1)=disstec_k(2)
          disstec_k(n3)=disstec_k(n3m)

CG        calcolo della griglia
      write(*,*) 'scrittura file tec.data'
c     scrittura file per tecplot
      open(11,file='diss_rtheta_avg.data',status='unknown')
c     intestazione: una zona
      write(11,*) 'VARIABLES = "Z","D","DT" '
      write(11,*) 'ZONE K=',n3,', F=BLOCK'

c     coordinata z
      write(11,*) (zz(k),k=1,n3)

c       diss D
      write(11,*) (dissc_k(k),k=1,n3)
c       temperature diss DT
      write(11,*) (disstec_k(k),k=1,n3)

      close(11)


   
       return
       end
c-----------------------------------------------------------------
       subroutine totaverage_rk(dissteme,dissme)
       include 'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       dimension dissteme(m1,m2,m3),disst_rk(m2,m3)
       dimension dissme(m1,m2,m3),diss_rk(m2,m3)
       dimension disstec_rk(m2,m3),dissc_rk(m2,m3)
       common/parava/ray,pra,pec,ren
       common /metrr/  g2rm(m2),g2rc(m2)
       common /ametrz/  g3rc(m3),g3rm(m3),str3
       common/cordvo/rc(m2),rm(m2),zz(m3)  
       common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
       dimension yp1(n1,n2),yp2(n1,n2)
       common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
       common/corrt/thetac(m1),thetam(m1)
	pi = 2.*asin(1.)
        vol = 1./ (2.*pi)
	usdx1 = 1./dx1
        do k=1,n3m
              do j=1,n2m
               etatemp=0.0
               eta=0.0
	       vol1=0.d0
                     fac2=usdx1
                     do i=1,n1m
                     etatemp=etatemp+dissteme(i,j,k)*fac2
                     eta=eta+dissme(i,j,k)*fac2
                     vol1=vol1+fac2
		     enddo	
	     diss_rk(j,k)=eta*vol
	     disst_rk(j,k)=etatemp*vol
c             print*,'vol',vol,1./vol1
              enddo
       enddo
c       print*,'vol',vol,vol1
       print*,'pec',pec
c	open(unit=42,file='dissipation_rk.out')
c	do k=1,n3m
c              do j=1,n2m
c	      write(42,*)zz(k),diss_rk(j,k),disst_rk(j,k)
c              enddo
c	enddo
c	close(42)

      do 4 kc=1,n3m 
        km=kmv(kc)                                                    
        do 4 jc=2,n2m                                                     
          jm=jmv(jc)                                                        
            dissc_rk(jc,kc)=(diss_rk(jc,kc)
     %                      +diss_rk(jm,kc)
     %                      +diss_rk(jc,km)
     %                      +diss_rk(jm,km))*0.25
              disstec_rk(jc,kc)=(disst_rk(jc,kc)
     %                      +disst_rk(jm,kc)
     %                      +disst_rk(jc,km)
     %                      +disst_rk(jm,km))*0.25
    4 continue   
      do 6 kc=1,n3m
        km=kmv(kc)                                                    
          dissc_rk(n2,kc)= (diss_rk(n2m,kc)
     %                      +diss_rk(n2m,km))*0.5d0
           disstec_rk(n2,kc)= (disst_rk(n2m,kc)
     %                      +disst_rk(n2m,km))*0.5d0
           dissc_rk(1,kc)=(diss_rk(1,kc)+diss_rk(1,km)
     %               )*0.5d0
           disstec_rk(1,kc)=(disst_rk(1,kc)+disst_rk(1,km)
     %               )*0.5d0
   6   continue

      do 15 jc=1,n2                                                    
        do 15 ic=1,n1
          dissc_rk(jc,1)=dissc_rk(jc,2)
          dissc_rk(jc,n3)=dissc_rk(jc,n3m)
           disstec_rk(jc,1)=disstec_rk(jc,2)
          disstec_rk(jc,n3)=disstec_rk(jc,n3m)
   15 continue  
c
c      =====for matlab=============
	open(52,file='dissc_matlab.dat')
	open(53,file='disstec_matlab.dat')
      do k=1,n3
            write(52,752)(dissc_rk(j,k),j=1,n2)
            write(53,752)(disstec_rk(j,k),j=1,n2)
      end do
	close(53)
	close(52)
 752  format(400f12.8)

CG    ====== for tecplot ===========
      do 1 i=1,n1
        do 1 j=1,n2
          yp1(i,j)=rc(j)*cos(thetac(i))
          yp2(i,j)=rc(j)*sin(thetac(i))
    1 continue 

      write(*,*) 'scrittura file tec.data'
c     scrittura file per tecplot
      open(11,file='diss_theta_avg.data',status='unknown')
c     intestazione: una zona
      write(11,*) 'VARIABLES = "Y","Z","D","DT" '
      write(11,*) 'ZONE I=',n1,', J=',n2,', K=',n3,', F=BLOCK'

c     coordinata x
      write(11,*) (((yp1(i,j),i=1,n1),j=1,n2),k=1,n3)
c     coordinata y
      write(11,*) (((yp2(i,j),i=1,n1),j=1,n2),k=1,n3)
c     coordinata z
      write(11,*) (((zz(k),i=1,n1),j=1,n2),k=1,n3)

c       diss D
      write(11,*) (((dissc_rk(j,k),i=1,n1),j=1,n2),k=1,n3)
c       diss Dt
      write(11,*) (((disstec_rk(j,k),i=1,n1),j=1,n2),k=1,n3)

      close(11)


       return
       end
c-----------------------------------------------------------------
       subroutine staggu2in3(u2,u23) 
       include'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       dimension u2(m1,m2,m3),up(m1,m2,m3),u23(m1,m2,m3)
       common /ametrz/  g3rc(m3),g3rm(m3),str3

       do k=2,n3m
              km=k-1
              usg3rm=1.0/(g3rm(k)+g3rm(km))
              do j=1,n2
                     do i=1,n1m
                     up(i,j,k)=(u2(i,j,k)*g3rm(k)+
     %                       u2(i,j,km)*g3rm(km))*usg3rm
                     enddo
              enddo
       enddo
       do k=2,n3m
              do j=1,n2m
              jp=j+1
                     do i=1,n1m
                     u23(i,j,k)=(up(i,j,k)+up(i,jp,k))*0.5
                     enddo
              enddo
       enddo
       k=1
       do j=1,n2m
              do i=1,n1m
              u23(i,j,k)=0.0
              enddo
       enddo
       k=n3
       do j=1,n2m
              do i=1,n1m
              u23(i,j,k)=0.0
              enddo
       enddo
       
       return   
      end 

c-----------------------------------------------------------------
       subroutine staggu3in2(u3,u32) 
       include'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
       dimension u3(m1,m2,m3),up(m1,m2,m3),u32(m1,m2,m3)
       common /metrr/  g2rm(m2),g2rc(m2)
       common/asym/isym(m1)
       do k=1,n3
              do j=2,n2m
              jm=j-1
              usg2rm=1.0/(g2rm(j)+g2rm(jm))
                     do i=1,n1m
                     up(i,j,k)=(u3(i,j,k)*g2rm(j)+
     %                    u3(i,jm,k)*g2rm(jm))*usg2rm
                     enddo
              enddo
       enddo
       do k=1,n3m
              kp=k+1
              do j=2,n2m
                     do i=1,n1m
                     u32(i,j,k)=(up(i,j,k)+up(i,j,kp))*0.5
                     enddo
              enddo
       enddo

c      sidewall
       j=n2
       do k=1,n3m
              do i=1,n1m
                     u32(i,j,k)=0.0
              enddo
       enddo
c      sull'asse       
       do k=1,n3m
              do i=1,n1m
              u32(i,1,k)=(u32(i,2,k)+u32(isym(i),2,k))*0.5d0
              enddo
       enddo

       return   
      end

c-----------------------------------------------------------------
       subroutine staggu1in3(u1,u13) 
       include'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       dimension u1(m1,m2,m3),up(m1,m2,m3),u13(m1,m2,m3)
       common /ametrz/  g3rc(m3),g3rm(m3),str3

       do k=2,n3m
              km=k-1
              usg3rm=1.0/(g3rm(k)+g3rm(km))
              do j=1,n2m
                     do i=1,n1
                     up(i,j,k)=(u1(i,j,k)*g3rm(k)+
     %                     u1(i,j,km)*g3rm(km))*usg3rm
                     enddo
              enddo
       enddo
       do k=2,n3m
              do j=1,n2m
                     do i=1,n1m
                     ip=i+1
                     u13(i,j,k)=(up(i,j,k)+up(ip,j,k))*0.5
                     enddo
              enddo
       enddo
       k=1
       do j=1,n2m
              do i=1,n1m
              u13(i,j,k)=0.0
              enddo
       enddo
       k=n3
       do j=1,n2m
              do i=1,n1m
              u13(i,j,k)=0.0
              enddo
       enddo

       return   
      end
c-----------------------------------------------------------------
       subroutine staggu3in1(u3,u31) 
       include'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       dimension u3(m1,m2,m3),up(m1,m2,m3),u31(m1,m2,m3)

       do k=1,n3
              do j=1,n2m
                     do i=2,n1m
                     im=i-1
                     up(i,j,k)=(u3(i,j,k)+u3(im,j,k))*0.5
                     enddo
              enddo
       enddo
       do k=1,n3m
              kp=k+1
              do j=1,n2m
                     do i=2,n1m
                     u31(i,j,k)=(up(i,j,k)+up(i,j,kp))*0.5
                     enddo
              enddo
       enddo
c       interfaccia       
       do k=1,n3
              do j=1,n2m
              up(1,j,k)=(u3(1,j,k)+u3(n1m,j,k))*0.5
              enddo
       enddo
       do k=1,n3m
              kp=k+1
              do j=1,n2m
              u31(1,j,k)=(up(1,j,k)+up(1,j,kp))*0.5
              u31(n1,j,k)=u31(1,j,k)
              enddo
       enddo
       return   
      end
c-----------------------------------------------------------------
       subroutine staggu1in2(u1,u12) 
       include'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
       dimension u1(m1,m2,m3),up(m1,m2,m3),u12(m1,m2,m3)
       common /metrr/  g2rm(m2),g2rc(m2)
       common/asym/isym(m1)
       do k=1,n3m
              do j=2,n2m
              jm=j-1
              usg2rm=1.0/(g2rm(j)+g2rm(jm))
                     do i=1,n1
                     up(i,j,k)=(u1(i,j,k)*g2rm(j)+
     %                  u1(i,jm,k)*g2rm(jm))*usg2rm
                     enddo
              enddo
       enddo
       do k=1,n3m
              do j=2,n2m
                     do i=1,n1m
                     ip=i+1
                     u12(i,j,k)=(up(i,j,k)+up(ip,j,k))*0.5
                     enddo
              enddo
       enddo
c      sidewall
       j=n2
       do k=1,n3m
              do i=1,n1m
                     u12(i,j,k)=0.0
              enddo
       enddo
c      sull'asse       
       do k=1,n3m
              do i=1,n1m
              u12(i,1,k)=(u12(i,2,k)-u12(isym(i),2,k))*0.5d0
              enddo
       enddo

       return   
      end
c-----------------------------------------------------------------
       subroutine staggu2in1(u2,u21) 
       include'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       dimension u2(m1,m2,m3),up(m1,m2,m3),u21(m1,m2,m3)

       do k=1,n3m
              do j=1,n2
                     do i=2,n1m
                     im=i-1
                     up(i,j,k)=(u2(i,j,k)+u2(im,j,k))*0.5
                     enddo
              enddo
       enddo
       do k=1,n3m
              do j=1,n2m
                     jp=j+1
                     do i=2,n1m
                     u21(i,j,k)=(up(i,j,k)+up(i,jp,k))*0.5
                     enddo
              enddo
       enddo
c       interfaccia       
       do k=1,n3m
              do j=1,n2
              up(1,j,k)=(u2(1,j,k)+u2(n1m,j,k))*0.5
              enddo
       enddo
       do k=1,n3m
              do j=1,n2m
              jp=j+1
              u21(1,j,k)=(up(1,j,k)+up(1,jp,k))*0.5
              u21(n1,j,k)=u21(1,j,k)
              enddo
       enddo

       return   
      end
c-----------------------------------------------------------------
       subroutine tecplot_diss(disste,diss)
       include'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
      common/corrt/thetac(m1),thetam(m1)
      common/cordvo/rc(m2),rm(m2),zz(m3)
       common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
       common/asym/isym(m1)  
       dimension disste(m1,m2,m3),diss(m1,m2,m3) 
       dimension disstec(m1,m2,m3),dissc(m1,m2,m3)
       dimension yp1(n1,n2),yp2(n1,n2)
      do 4 kc=1,n3m 
        km=kmv(kc)                                                    
        do 4 jc=2,n2m                                                     
          jm=jmv(jc)                                                        
          do 4 ic=1,n1m                                                     
            im=imv(ic)                                                        
            dissc(ic,jc,kc)=(diss(ic,jc,kc)+diss(im,jc,kc)
     %                      +diss(ic,jm,kc)+diss(im,jm,kc)
     %                      +diss(ic,jc,km)+diss(im,jc,km)
     %                      +diss(ic,jm,km)+diss(im,jm,km))*0.125
              disstec(ic,jc,kc)=(disste(ic,jc,kc)+disste(im,jc,kc)
     %                      +disste(ic,jm,kc)+disste(im,jm,kc)
     %                      +disste(ic,jc,km)+disste(im,jc,km)
     %                      +disste(ic,jm,km)+disste(im,jm,km))*0.125
    4 continue   
c
c
      do 6 kc=1,n3m
        km=kmv(kc)                                                    
        do 61 ic=1,n1m                                                     
          dissc(ic,n2,kc)= (diss(ic,n2m,kc)+diss(im,n2m,kc)
     %                      +diss(ic,n2m,km)+diss(im,n2m,km))*0.25
           disstec(ic,n2,kc)= (disste(ic,n2m,kc)+disste(im,n2m,kc)
     %                      +disste(ic,n2m,km)+disste(im,n2m,km))*0.25
          dissc(ic,1,kc)=(diss(ic,1,kc)+diss(ic,1,km)+
     %               (diss(isym(ic),1,kc)+diss(isym(ic),1,km)))*0.25d0
          disstec(ic,1,kc)=(disste(ic,1,kc)+disste(ic,1,km)+
     %               (disste(isym(ic),1,kc)+disste(isym(ic),1,km)))*0.25d0
   61   continue
    6 continue  
c
       do kc=1,n3m
              do jc=1,n2
                 dissc(n1,jc,kc)= dissc(1,jc,kc)
              disstec(n1,jc,kc)= disstec(1,jc,kc)
              enddo
       enddo              
      do 15 jc=1,n2                                                    
        do 15 ic=1,n1
          dissc(ic,jc,1)=dissc(ic,jc,2)
          dissc(ic,jc,n3)=dissc(ic,jc,n3m)
           disstec(ic,jc,1)=disstec(ic,jc,2)
          disstec(ic,jc,n3)=disstec(ic,jc,n3m)
   15 continue  

c


CG        calcolo della griglia
      do 1 i=1,n1
        do 1 j=1,n2
          yp1(i,j)=rc(j)*cos(thetac(i))
          yp2(i,j)=rc(j)*sin(thetac(i))
    1 continue 

      write(*,*) 'scrittura file tec.data'
c     scrittura file per tecplot
      open(11,file='diss.data',status='unknown')
c     intestazione: una zona
      write(11,*) 'VARIABLES = "X","Y","Z","D","DT" '
      write(11,*) 'ZONE I=',n1,', J=',n2,', K=',n3,', F=BLOCK'

c     coordinata x
      write(11,*) (((yp1(i,j),i=1,n1),j=1,n2),k=1,n3)
c     coordinata y
      write(11,*) (((yp2(i,j),i=1,n1),j=1,n2),k=1,n3)
c     coordinata z
      write(11,*) (((zz(k),i=1,n1),j=1,n2),k=1,n3)

c       diss D
      write(11,*) (((dissc(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c       temperature diss DT
      write(11,*) (((disstec(i,j,k),i=1,n1),j=1,n2),k=1,n3)

      close(11)


       return   
      end
c---------------------------------------------------------------------
       subroutine tecplot_boum(u1,u2,u3,dens)
       include 'parade.f'
       common/dim/n1,n1m,n2,n2m,n3,n3m 
       common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q 
      common/corrt/thetac(m1),thetam(m1)
      common/cordvo/rc(m2),rm(m2),zz(m3)   
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
       common/asym/isym(m1)
      dimension u1(m1,m2,m3),u2(m1,m2,m3)  
      dimension u3(m1,m2,m3), dens(m1,m2,m3)
      dimension voc1(m1,m2,m3),voc2(m1,m2,m3)
      dimension qx(m1,m2,m3),qy(m1,m2,m3)
      dimension densc(m1,m2,m3),voc3(m1,m2,m3)
 
c
      do 4 kc=1,n3m 
        km=kmv(kc)                                                    
        do 4 jc=2,n2m                                                     
          jm=jmv(jc)                                                        
          do 4 ic=1,n1m                                                     
            im=imv(ic)                                                        
            voc2(ic,jc,kc)=(u2(ic,jc,kc)+u2(im,jc,kc)+
     %                      u2(ic,jc,km)+u2(im,jc,km))*0.25
            voc1(ic,jc,kc)=(u1(ic,jc,kc)+u1(ic,jm,kc)+
     %                      u1(ic,jc,km)+u1(ic,jm,km))*0.25
            voc3(ic,jc,kc)=(u3(ic,jc,kc)+u3(im,jc,kc)+
     %                      u3(ic,jm,kc)+u3(im,jm,kc))*0.25 
!            densc(ic,jc,kc)=(dens(ic,jc,kc)+dens(im,jc,kc)
!     %                      +dens(ic,jm,kc)+dens(im,jm,kc)
!     %                      +dens(ic,jc,km)+dens(im,jc,km)
!     %                      +dens(ic,jm,km)+dens(im,jm,km))*0.125
    4 continue   
c
c
      do 6 kc=1,n3m
        km=kmv(kc)                                                    
        do 61 ic=1,n1m                                                     
          voc1(ic,n2,kc)=0.d0
          voc2(ic,n2,kc)=0.d0
          voc3(ic,n2,kc)=0.d0
!          densc(ic,n2,kc)= (dens(ic,n2m,kc)+dens(im,n2m,kc)
!     %                      +dens(ic,n2m,km)+dens(im,n2m,km))*0.25
          voc3(ic,1,kc)=(u3(ic,1,kc)+u3(isym(ic),1,kc))*0.5d0
          voc2(ic,1,kc)=(u2(ic,2,kc)+u2(ic,2,km)-
     %                    (u2(isym(ic),2,km)+u2(isym(ic),2,km)))*0.25d0
          voc1(ic,1,kc)=(u1(ic,1,kc)+u1(ic,1,km)-
     %                    (u1(isym(ic),1,km)+u1(isym(ic),1,km)))*0.25d0
!          densc(ic,1,kc)=(dens(ic,1,kc)+dens(ic,1,km)+
!     %               (dens(isym(ic),1,kc)+dens(isym(ic),1,km)))*0.25d0
   61   continue
    6 continue  
       

       
       
       do kc=1,n3m
              do jc=1,n2
              voc1(n1,jc,kc)=voc1(1,jc,kc)
                 voc2(n1,jc,kc)=voc2(1,jc,kc)
                 voc3(n1,jc,kc)=voc3(1,jc,kc)
!                 densc(n1,jc,kc)= densc(1,jc,kc)
              enddo
       enddo       
c
      do 15 jc=1,n2                                                    
        do 15 ic=1,n1
          voc3(ic,jc,n3)=0.d0           
          voc3(ic,jc,1)=0.d0            
!          densc(ic,jc,1)=1.d0
!          densc(ic,jc,1)=0.d0
!!          densc(ic,jc,n3)=0.d0
!          densc(ic,jc,n3)=1.d0
            voc2(ic,jc,1)=0.d0            
            voc1(ic,jc,1)=0.d0            
            voc2(ic,jc,n3)=0.d0            
!            voc1(ic,jc,n3)=0.d0            
            voc1(ic,jc,n3)=rc(jc)/rc(n2)            
   15 continue  
c
c        
       call tec_vel(voc1,voc2,voc3,densc)
      return
      end
c -----------------------------------------------------
      subroutine tec_vel(voc1,voc2,voc3,densc)
      include 'parade.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m  
       dimension yp1(m1,m2),yp2(m1,m2),rs(m2)
      common/corrt/thetac(m1),thetam(m1)
      common/cordvo/rc(m2),rm(m2),zz(m3) 
      dimension voc1(m1,m2,m3),voc2(m1,m2,m3),voc3(m1,m2,m3)
      dimension qx(m1,m2,m3),qy(m1,m2,m3)
      dimension qx_sec(m1,m2,m3),qy_sec(m1,m2,m3)
      dimension densc(m1,m2,m3)

      character*70 namfi3

         pi=2.*asin(1.)
CG        calcolo della griglia
      do 1 i=1,n1
        do 1 j=1,n2
          yp1(i,j)=rc(j)*cos(thetac(i))
          yp2(i,j)=rc(j)*sin(thetac(i))
    1 continue 

!   transformation of velocity from polar to cartesian
      do k=1,n3
        do j=1,n2
          do i=1,n1
            qx(i,j,k)=voc2(i,j,k)*cos(thetac(i))
     %                  -voc1(i,j,k)*sin(thetac(i))
            qy(i,j,k)=voc2(i,j,k)*sin(thetac(i))
     %                  +voc1(i,j,k)*cos(thetac(i))
!   components of velocity in a plane
            qx_sec(i,j,k)=voc2(i,j,k)*cos(thetac(i))
            qy_sec(i,j,k)=voc2(i,j,k)*sin(thetac(i))
          enddo
        enddo
      enddo

!        do j=1,n2
!          rs(j)=rc(j)/rc(n2)
!        enddo
!
!      open(151,file='profile.dat',form='formatted',status='unknown')
!        do j=1,n2
!          write(151,152)rs(j),(voc3s(12,j,k),k=1,n3)
!        enddo
!      close(151)
!152       format(500(F10.5,2X))
!
!       open(555,file='w_axis.dat')
!         do k=1,n3m
!           write(555,*)zz(k),voc3s(12,1,k)
!         enddo
!       close(555)
       open(666,file='w_azimuth39.dat')
         do i=1,n1
           write(666,*)thetac(i),voc3(i,39,43)
         enddo
       close(666)
       open(666,file='w_azimuth53.dat')
         do i=1,n1
           write(666,*)thetac(i),voc3(i,53,43)
         enddo
       close(666)
       open(666,file='w_azimuth65.dat')
         do i=1,n1
           write(666,*)thetac(i),voc3(i,65,45)
         enddo
       close(666)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     writing the domain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      open(121,file='mean_field.dat',status='unknown')
!c     intestazione: una zona
!c     coordinata x
!      write(121,*) (((yp1(i,j),i=1,n1),j=1,n2),k=1,n3)
!c     coordinata y
!      write(121,*) (((yp2(i,j),i=1,n1),j=1,n2),k=1,n3)
!c     coordinata z
!      write(121,*) (((zz(k),i=1,n1),j=1,n2),k=1,n3)
!!      write(121,*) (((qx(i,j,k),i=1,n1),j=1,n2),k=1,n3)
!      write(121,*) (((voc1(i,j,k),i=1,n1),j=1,n2),k=1,n3)
!c     velocita''V
!!      write(121,*) (((qy(i,j,k),i=1,n1),j=1,n2),k=1,n3)
!      write(121,*) (((voc2(i,j,k),i=1,n1),j=1,n2),k=1,n3)
!c       velocita' W
!      write(121,*) (((voc3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
!      close(121)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'scrittura file tec.data'
c     scrittura file per tecplot
      open(11,file='medie_tec.data',status='unknown')
c     intestazione: una zona
      write(11,*) 'VARIABLES = "X","Y","Z","U","V","W"'
      write(11,*) 'ZONE I=',n1,', J=',n2,', K=',n3,', F=BLOCK'

c     coordinata x
      write(11,234) (((yp1(i,j),i=1,n1),j=1,n2),k=1,n3)
c     coordinata y
      write(11,234) (((yp2(i,j),i=1,n1),j=1,n2),k=1,n3)
c     coordinata z
      write(11,234) (((zz(k),i=1,n1),j=1,n2),k=1,n3)

!     2D section of Q_x
!      write(11,*) (((qx_sec(i,j,k),i=1,n1),j=1,n2),k=1,n3)
!     2D section of Q_y
!      write(11,*) (((qy_sec(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     velocita' U
      write(11,234) (((qx(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     velocita''V
      write(11,234) (((qy(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c       velocita' W
      write(11,234) (((voc3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c       temperatura T
!      write(11,234) (((densc(i,j,k),i=1,n1),j=1,n2),k=1,n3)
      close(11)
       write(*,*) 'scritto file tec.data'
234     format(1000(f16.9,2x))


!  Script for .Vtk format (mayavi)
!
!        open(20,file='vel-boum.vtk')
!         write(20,67)
!67      format('# vtk DataFile Version 2.0')
!         write(20,68)
!68       format('annular example')
!         write(20,69)
!69       format('ASCII')
!!          write(20,71)
!71        format('DATASET STRUCTURED_GRID')
!73       format('DIMENSIONS',I8,1x,I8,1x,I8)
!74       format('POINTS',I8,1x,'float')
!          write(20,73) n1,n2,n3
!          write(20,74) n1*n2*n3
!          do k=1,n3
!          do j=1,n2
!!          do i=1,n1
!          write(20,*) yp1(i,j),yp2(i,j),zz(k)
!!          write(20,*) thetac(i),rc(j),zz(k)
!          enddo
!          enddo
!!          enddo
!          write(20,*) 'POINT_DATA', n1*n2*n3
!          write(20,*)'SCALARS  vtheta  float '
!          write(20,*)'LOOKUP_TABLE  default'
!c75        format('SCALARS  vtheta  float ') 
!!c76        format('LOOKUP_TABLE  default') 
!          do k=1,n3
!          do j=1,n2
!          do i=1,n1
!          write(20,*)qx(i,j,k)
!          enddo
!          enddo
!          enddo
!          write(20,*)'SCALARS  vradial  float '
!          write(20,*)'LOOKUP_TABLE  default'
!c75        format('SCALARS  vradial  float ') 
!c76        format('LOOKUP_TABLE  default') 
!          do k=1,n3
!          do j=1,n2
!          do i=1,n1
!          write(20,*)qy(i,j,k)
!          enddo
!!          enddo
!          enddo
!          write(20,*)'SCALARS  vaxial  float '
!          write(20,*)'LOOKUP_TABLE  default'
!
!c75        format('SCALARS  vaxial  float ') 
!c76        format('LOOKUP_TABLE  default') 
!          do k=1,n3
!          do j=1,n2
!!          do i=1,n1
!!          write(20,*)voc3(i,j,k)
!!          enddo
!!          enddo
!!          enddo
!!          write(20,*)'SCALARS  dens  float '
!!          write(20,*)'LOOKUP_TABLE  default'
!c75        format('SCALARS  dens  float ') 
!c76        format('LOOKUP_TABLE  default') 
!!          do k=1,n3
!!          do j=1,n2
!!          do i=1,n1
!!          write(20,*)densc(i,j,k)
!!          enddo
!!          enddo
!!          enddo
!!      write(20,*)'VECTORS','velocity','float'
!      write(20,*)'VECTORS  velocity  float'
!        do k=1,n3
!          do j=1,n2
!            do i=1,n1
!               write(20,*)qx(i,j,k),qy(i,j,k),voc3(i,j,k)
!            enddo
!          enddo
!        enddo
!
!        close(20)
      return
      end
