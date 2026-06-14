!  The Fortran example to call  the shared library
!   Windows with MSYS2 install:
!      pacman -S mingw-w64-x86_64-gcc-fortran
!      make -f mf.mk
!

module seuif97
   implicit none
   Interface
      real(c_double) function pt(p,t,wid) bind (C,name="pt")
         use iso_c_binding
         real(c_double), value :: p,t
         integer(c_int), value ::wid
      end function pt

      real(c_double) function pt2s(p,t) bind (C,name="pt2s")
         use iso_c_binding
         real(c_double), value :: p,t
      end function pt2s

   End Interface

end module

program demo
   use iso_c_binding
   use seuif97
   implicit none
   real(c_double) :: p,t,h,s
   p = 16.13;
   t = 535.0;
   ! universal property functions with o_id parameter
   h = pt(p, t, 4);
   ! direct property functions
   s = pt2s(p, t);
   write (*,'(A,F10.2,F10.2,F10.2,F10.4)')  "(p,t) ->h,s",p,t,h,s
end program demo
