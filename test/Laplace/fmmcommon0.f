c
c
c
      subroutine mpalloc0(nd,laddr,iaddr,nlevels,lmptot,nterms,nboxes)
c     This subroutine determines the size of the array
c     to be allocated for the multipole expansions
c     iaddr(1,i) points to the starting location of the multipole
c     expansion of box i and iaddr(2,i) points to the local
c     expansion of box i
c  
c     Input arguments
c     nd          in: Integer
c                 number of expansions per box            
c 
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array provinding access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     nterms      in: Integer(0:nlevels)
c                 Number of terms requried in expansions at each
c                 level
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr       out: Integer *8(2,nboxes)
c                 Points the multipole and local expansions in box i
c 
c     lmptot      out: Integer *8
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,nterms(0:nlevels),nd,nboxes
      integer iaddr(2,nboxes), lmptot 
      integer laddr(2,0:nlevels)
      integer ibox,i,iptr
      integer istart,nn,itmp

      istart = 1
      do i = 0,nlevels

        nn = (2*nterms(i)+1)*2*(nterms(i)+1)*nd 

         do ibox = laddr(1,i),laddr(2,i)
c            Allocate memory for the multipole expansion         

             itmp = ibox-laddr(1,i)
             iaddr(1,ibox) = istart + itmp*2*nn

c            Allocate memory for the local expansion
             iaddr(2,ibox) = istart + itmp*2*nn + nn
         enddo       
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*2*nn
      enddo
      lmptot = istart

      return
      end
c----------------------------------------------------------------     

