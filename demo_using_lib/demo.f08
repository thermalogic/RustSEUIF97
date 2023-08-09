!
!  The Fortran example to call  the if97 library
!     make -f mfortran.mk
!

module if97
   implicit none
   Interface
      real(c_double) function pt(p,t,wid) bind (C,name="pt")
         use iso_c_binding
         real(c_double), value :: p,t
         integer(c_int), value ::wid
      end function pt

   End Interface

end module

program demo
   use iso_c_binding
   use if97
   implicit none
   real(c_double) :: p,t,h,s,v
   p = 16.13;
   t = 535.0;

   h = pt(p, t, 4);
   s = pt(p, t, 5);
   v = pt(p, t, 3);
   write (*,'(A,F10.2,F10.2,F10.2,F10.4,F10.4)') "(p,t),h,s,v",p,t,h,s,v
end program demo
