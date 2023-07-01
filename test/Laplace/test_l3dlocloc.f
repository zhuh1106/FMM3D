      implicit None
      real *8 ztrg(3),source(3)
      real *8 c1(3),c2(3)
      real *8, allocatable :: dc(:,:)
      real *8, allocatable :: pots(:,:),flds(:,:,:),hesss(:,:,:)
      real *8, allocatable :: opots(:,:),oflds(:,:,:),ohesss(:,:,:)
      complex *16, allocatable :: locexp1(:,:,:),locexp2(:,:,:)
      real *8, allocatable :: charge(:,:),dipvec(:,:,:)
      real *8 wlege(100000),err_out(3)
      real *8 scarray_loc(100000)
      real *8 bsize, rscale1, rscale2, shift
      integer nd,idim,i,ns,nt,nlege,nterms1,nterms2,lw7,lused7,nn
      real *8 eps, thresh

c     set scale for center 1 and center 2
      bsize = 1.0d0
      rscale1 = bsize/4
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
        charge(idim,1) = 1.1d0
        dipvec(idim,1,1) = 0.31
        dipvec(idim,2,1) = 0.33
        dipvec(idim,3,1) = 0.43
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
      allocate(ohesss(nd,6,nt),hesss(nd,6,nt))

c
c     direct calculation:
c
      do i=1,nt
        do idim=1,nd
          opots(idim,i) = 0
          oflds(idim,1,i) = 0
          oflds(idim,2,i) = 0
          oflds(idim,3,i) = 0
          ohesss(idim,1,i) = 0
          ohesss(idim,2,i) = 0
          ohesss(idim,3,i) = 0
          ohesss(idim,4,i) = 0
          ohesss(idim,5,i) = 0
          ohesss(idim,6,i) = 0
        enddo
      enddo

      thresh = 1.0d-15
      call l3ddirectch(nd,source,charge,ns,ztrg,
     1    nt,opots,oflds,ohesss,thresh)

c     local exp order for center 1
      nterms1 = 39
c     local exp order for center 2
      nterms2 = 7

c     precompute
      nlege = 100
      lw7 = 100000
      call ylgndrfwini(nlege,wlege,lw7,lused7)
      nn = nterms1
      nn = max(nn,nterms2)
      nn = 4*nn
      allocate(dc(0:nn,0:nn))
      call getsqrtbinomialcoeffs(nn,dc)

c     form local exp for center 1
      allocate(locexp1(nd,0:nterms1,-nterms1:nterms1))
      call mpzero(nd,locexp1,nterms1)
      call l3dformtac(nd, rscale1, source, charge,
     1     ns,c1,nterms1,locexp1,wlege,nlege)

      do i=1,nt
        do idim=1,nd
          pots(idim,i) = 0
          flds(idim,1,i) = 0
          flds(idim,2,i) = 0
          flds(idim,3,i) = 0
          hesss(idim,1,i) = 0
          hesss(idim,2,i) = 0
          hesss(idim,3,i) = 0
          hesss(idim,4,i) = 0
          hesss(idim,5,i) = 0
          hesss(idim,6,i) = 0
        enddo
      enddo

c     eval loc exp1 at center 1
      call l3dtaevalhessdini(nterms1,scarray_loc)
      call l3dtaevalh(nd,rscale1,c1,locexp1,nterms1,ztrg,nt,
     1      pots,flds,hesss,scarray_loc)

      print *,"local exp1 eval error:"
      call errprinth(nd,nt,pots,opots,flds,oflds,hesss,ohesss,
     1   err_out)
      


cc    shift local exp from center 1 to center 2
      allocate(locexp2(nd,0:nterms2,-nterms2:nterms2))
      call mpzero(nd,locexp2,nterms2)
      call l3dlocloc(nd,rscale1,c1,locexp1,nterms1,
     1      rscale2,c2,locexp2,nterms2,dc,nn)

      do i=1,nt
        do idim=1,nd
          pots(idim,i) = 0
          flds(idim,1,i) = 0
          flds(idim,2,i) = 0
          flds(idim,3,i) = 0
          hesss(idim,1,i) = 0
          hesss(idim,2,i) = 0
          hesss(idim,3,i) = 0
          hesss(idim,4,i) = 0
          hesss(idim,5,i) = 0
          hesss(idim,6,i) = 0
        enddo
      enddo

c     eval loc exp2 at center 2
      call l3dtaevalhessdini(nterms2, scarray_loc)
      call l3dtaevalh(nd,rscale2,c2,locexp2,nterms2,ztrg,nt,
     1      pots,flds,hesss,scarray_loc)

      print *,"local exp2 eval error:"
      call errprinth(nd,nt,pots,opots,flds,oflds,hesss,ohesss,
     1   err_out)

      stop
      end

c
      subroutine errprint(pot,opot,fld,ofld,errs)
      implicit real *8 (a-h,o-z)
      real *8 pot,opot,fld(3),ofld(3)
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

c
      subroutine errprinth(nd,nt,pot,opot,fld,ofld,hess,ohess,errs)
      implicit real *8 (a-h,o-z)
      real *8 pot(nd,nt),opot(nd,nt),fld(nd,3,nt)
      real *8 ofld(nd,3,nt),hess(nd,6,nt),ohess(nd,6,nt)
      real *8 errs(3)
 1000  format(4D15.5) 
      
      errp = 0
      dddp = 0
      errf = 0
      dddf = 0
      errh = 0
      dddh = 0

      do i=1,nt
        do idim=1,nd
          errp = errp + abs(pot(idim,i)-opot(idim,i))**2
          errf = errf + abs(fld(idim,1,i)-ofld(idim,1,i))**2
          errf = errf + abs(fld(idim,2,i)-ofld(idim,2,i))**2
          errf = errf + abs(fld(idim,3,i)-ofld(idim,3,i))**2

          do j=1,6
            errh = errh + abs(hess(idim,j,i)-ohess(idim,j,i))**2
            dddh = dddh + abs(ohess(idim,j,i))**2
          enddo

          dddp = dddp + abs(opot(idim,i))**2
          dddf = dddf + abs(ofld(idim,1,i))**2
          dddf = dddf + abs(ofld(idim,2,i))**2
          dddf = dddf + abs(ofld(idim,3,i))**2
        enddo
      enddo

c

      errs(1) = sqrt(errp/dddp)
      errs(2) = sqrt(errf/dddf)
      errs(3) = sqrt(errh/dddh)
      write(*,'(a,e11.4,a,e11.4,a,e11.4)') 
     1     'pot error=',errs(1),'   grad error=',errs(2),
     2     '   hess error=',errs(3)   
      return
      end
c
