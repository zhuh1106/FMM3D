      implicit None
      real *8 ztrg(3),source(3)
      real *8 c1(3),c2(3)
      real *8 xnodes(2000), wts(2000)
      complex *16, allocatable :: pots(:,:),flds(:,:,:)
      complex *16, allocatable :: opots(:,:),oflds(:,:,:)
      complex *16, allocatable :: locexp1(:,:,:),locexp2(:,:,:)
      complex *16, allocatable :: charge(:,:),dipvec(:,:,:)
      real *8 wlege(100000),err_out(2)
      real *8 bsize, rscale1, rscale2, shift
      integer nd,idim,i,ns,nt,nlege,nterms1,nterms2,lw7,lused7,nn
      integer ifinit
      real *8 eps, thresh, radius
      complex *16 zk, eye

c     set scale for center 1 and center 2
      data eye/(0.0d0,1.0d0)/
      zk = 1.02d0
      bsize = 1.0d0
      rscale1 = abs(zk)*bsize/4
      rscale2 = rscale1
      nd = 1

      allocate(charge(nd,1),dipvec(nd,3,1))

c     local exp center 1
      c1(1)=0.248d0
      c1(2)=0.251d0
      c1(3)=0.249d0
c     source position
      source(1)=c1(1)-0.47
      source(2)=c1(2)+0.51
      source(3)=c1(3)+0.49
      do idim=1,nd
        charge(idim,1) = 1.1d0 + eye*0.1
        dipvec(idim,1,1) = 0.31 + eye*1.0
        dipvec(idim,2,1) = 0.33 + eye*1.2
        dipvec(idim,3,1) = 0.43 + eye*0.3
      enddo
      ns = 1

c     local exp center 2
      shift = 0.00001
      c2(1) = c1(1) + shift
      c2(2) = c1(2) + shift*1.01
      c2(3) = c1(3) - shift*0.99

c     target position
      ztrg(1) = c1(1) + 0.01
      ztrg(2) = c1(2) + 0.02
      ztrg(3) = c1(3) - 0.01
      nt = 1

      allocate(opots(nd,nt),pots(nd,nt))
      allocate(oflds(nd,3,nt),flds(nd,3,nt))

c
c     direct calculation:
c
      do i=1,nt
        do idim=1,nd
          opots(idim,i) = 0
          oflds(idim,1,i) = 0
          oflds(idim,2,i) = 0
          oflds(idim,3,i) = 0
        enddo
      enddo

      thresh = 1.0d-15
      call h3ddirectcdg(nd,zk,source,charge,dipvec,ns,ztrg,
     1     nt,opots,oflds,thresh)

c     local exp order for center 1
      nterms1 = 39
c     local exp order for center 2
      nterms2 = 7

c     precompute
      nlege = 100
      lw7 = 100000
      call ylgndrfwini(nlege,wlege,lw7,lused7)
      ifinit = 1
      nn = nterms1
      nn = max(nn,nterms2)
      nn = 2*nn
      call legewhts(nn,xnodes,wts,ifinit)

c     form local exp for center 1
      allocate(locexp1(nd,0:nterms1,-nterms1:nterms1))
      call mpzero(nd,locexp1,nterms1)
      call h3dformtacd(nd, zk, rscale1, source, charge, dipvec,
     1     ns,c1,nterms1,locexp1,wlege,nlege)

      do i=1,nt
        do idim=1,nd
          pots(idim,i) = 0
          flds(idim,1,i) = 0
          flds(idim,2,i) = 0
          flds(idim,3,i) = 0
        enddo
      enddo

c     eval loc exp1 at center 1
      call h3dtaevalg(nd,zk,rscale1,c1,locexp1,nterms1,ztrg,nt,
     1      pots,flds,wlege,nlege)

      print *,"local exp1 eval error:"
      call errprint(pots,opots,flds,oflds,err_out)
      


cc    shift local exp from center 1 to center 2
      radius = sqrt(3.0d0)/4.0d0
      allocate(locexp2(nd,0:nterms2,-nterms2:nterms2))
      call mpzero(nd,locexp2,nterms2)
      call h3dlocloc(nd,zk,rscale1,c1,locexp1,nterms1,
     1      rscale2,c2,locexp2,nterms2,radius,xnodes,wts,nn)
c
      do i=1,nt
        do idim=1,nd
          pots(idim,i) = 0
          flds(idim,1,i) = 0
          flds(idim,2,i) = 0
          flds(idim,3,i) = 0
        enddo
      enddo
c
c     eval loc exp2 at center 2
      call h3dtaevalg(nd,zk,rscale2,c2,locexp2,nterms2,ztrg,nt,
     1      pots,flds,wlege,nlege)
c
      print *,"local exp2 eval error:"
      call errprint(pots,opots,flds,oflds,err_out)

      stop
      end

c
      subroutine errprint(pot,opot,fld,ofld,errs)
      implicit real *8 (a-h,o-z)
      complex *16 pot,opot,fld(3),ofld(3)
      real *8 errs(2)
 1000  format(4D15.5) 
      err = 0
      ddd = 0
      err = err + abs(fld(1)-ofld(1))**2
      err = err + abs(fld(2)-ofld(2))**2
      err = err + abs(fld(3)-ofld(3))**2
      ddd = ddd + abs(ofld(1))**2
      ddd = ddd + abs(ofld(2))**2
      ddd = ddd + abs(ofld(3))**2
      err = sqrt(err)
      ddd = sqrt(ddd)

      err1 = abs(pot-opot)/abs(opot)
      write(*,'(a,e11.4,a,e11.4)') 
     1     'pot error=',err1,'   grad error=',err/ddd

      errs(1) = err1 
      errs(2) = err/ddd

      return
      end
