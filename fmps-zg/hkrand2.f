c         random number generator
c
        subroutine hkrand2(iseed_hk,hkrand)
        implicit double precision (a-h,o-z)
        integer, intent(in) :: iseed_hk
        real *8, intent(out) :: hkrand
        dimension iseed(4)
        data iseed/0,0,0,1/
        save
c
        if( iseed_hk .ne. 0 ) then
        iseed(1)=0
        iseed(2)=0
        iseed(3)=0
        iseed(4)=mod(2*iseed_hk+1,4096)
        endif
c
        hkrand=dlaran(iseed)
c
        end
