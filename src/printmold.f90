! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!cccccccccccccccccccccccccccccccc
!    write out Molden input     c
!cccccccccccccccccccccccccccccccc
! ncent  : # atoms
! nmo    : # MOs
! nbf    : # AOs
! nprims : # primitives (in total)
! xyz(4,ncent) : Cartesian coordinates & nuclear charge
! cxip(nprims) : contraction coefficients of primitives
! exip(nprims) : exponents of primitives
! cmo(nbf,nmo) : LCAO-MO coefficients
! eval(nmo)    : orbital eigenvalues
! occ(nmo)     : occupation # of MO
! ipty(nprims) : angular momentum of primitive function
! ipao(nbf)    : # primitives in contracted AO
! ibf(ncent)   : # of contracted AOs on atom


subroutine printmold(ncent,nmo,nbf,xyz,at,cmo,eval,occ,thr,basis)
   use xtb_mctc_symbols, only : toSymbol
   use xtb_type_basisset
   implicit none
   type(TBasisset), intent(in) :: basis
   real*8, intent ( in ) :: xyz(3,ncent)
   real*8, intent ( in ) :: eval(nmo)
   real*8, intent ( in ) :: occ (nmo)
   real*8, intent ( in ) :: cmo(nbf,nmo)
   real*8, intent ( in ) :: thr
   integer, intent( in ) :: at(ncent)
   integer, intent( in ) :: ncent,nmo,nbf
   ! temporary variables
   integer i,j,k,icount,jcount,z,nmomax,nop
   real*8 dum
   character*1 aang
   character*2 atyp
   logical skip
   integer :: iwfn

   iwfn=29
   call open_file(iwfn,'molden.input','w')


   write(iwfn,'(A)',advance='yes')'[Molden Format]'
   write(iwfn,'(A)',advance='yes')'[Title]'

   !cccccccccccccccccccccccccccccccccccccccccccccccc
   ! print out atoms & coordinates                 c
   !cccccccccccccccccccccccccccccccccccccccccccccccc

   write(iwfn,'(A)',advance='yes')'[Atoms] AU'
   ! either coordinates are given in a.u.
   do i = 1,ncent
      atyp = toSymbol(at(i))
      z=at(i)
      ! atom character, running number, nuclear charge, x,y,z coordinates
      write (iwfn,'(a2,2i6,3E22.14)') &
         &   atyp,i,z,xyz(1,i),xyz(2,i),xyz(3,i)
   enddo
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Now print basis set data                                   c
   !                                                            c
   ! go through atoms and assign orbitals (contracted & prims)  c
   ! search for [GTO] statement to read basis set data          c
   ! for each set of atomic orbitals                            c
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   write(iwfn,'(A)',advance='yes') '[GTO]'
   !ccccccccccccccccccc
   ! go through atoms c
   !ccccccccccccccccccc
   icount=1
   do i=1,ncent
      write(iwfn,*) i,'0' ! I don't know what the zero is needed for
      !  now go trough nbfs located on atom i
      do j=basis%fila(1,i),basis%fila(2,i)
         call ang2chr(basis%lao(j),aang,skip)
         if(skip)then
            do k=1,basis%nprim(j)
               icount=icount+1
            enddo
         else
            write(iwfn,*) aang,basis%nprim(j),1.00
            do k=1,basis%nprim(j)
               write(iwfn,*) basis%alp(icount),basis%cont(icount)
               icount=icount+1
            enddo
         endif
      enddo
      write(iwfn,*)
   enddo

   nmomax=0
   do i=1,nmo
      if(eval(i).gt.thr.and.nmomax.eq.0)nmomax=i-1
   enddo
   if(nmomax.eq.0) nmomax=nmo

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! print occupation number of orbitals and orbital energies c
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   write(iwfn,'(A)',advance='yes') '[MO]'
   ! restricted xtb_printout
   do i=1,nmomax
      ! MO info
      write(iwfn,'(A)',advance='no') 'Sym= '
      write(iwfn,'(i5,a1)') i,'a'
      write(iwfn,'(A)',advance='no') 'Ene= '
      dum=0
      write(iwfn,*) eval(i)
      write(iwfn,'(A)',advance='no') 'Spin= '
      write(iwfn,'(A)',advance='yes') 'Alpha' ! for now just consider RHF case
      write(iwfn,'(A)',advance='no') 'Occup= '
      write(iwfn,'(F14.8)') occ(i)
      !now coefficients
      ! a notion: for l>0, the ordering is p: x,y,z ; d: xx,yy,zz,xy,xz,yz ; f: xxx, yyy, zzz, xxy, xxz, yyx, yyz, xzz, yzz, xyz
      do j=1,nbf
         write(iwfn,*) j,cmo(j,i)
      enddo
   enddo
   call close_file(iwfn)
end subroutine

! true U version
subroutine printumold(ncent,nmo,nbf,xyz,at,cmoa,cmob,evala,&
      &                      evalb,occa,occb,thr,basis)
   use xtb_mctc_symbols, only : toSymbol
   use xtb_type_basisset
   implicit none
   type(TBasisset), intent(in) :: basis
   real*8, intent ( in ) :: xyz(3,ncent)
   real*8, intent ( in ) :: evala(nmo)
   real*8, intent ( in ) :: evalb(nmo)
   real*8, intent ( in ) :: occa(nmo)
   real*8, intent ( in ) :: occb(nmo)
   real*8, intent ( in ) :: cmoa(nbf,nmo)
   real*8, intent ( in ) :: cmob(nbf,nmo)
   real*8, intent ( in ) :: thr
   integer, intent( in ) :: at(ncent)
   integer, intent( in ) :: ncent,nmo,nbf
   ! temporary variables
   integer i,j,k,icount,jcount,z,nmomax,nop
   real*8 dum
   character*1 aang
   character*2 atyp
   logical skip
   integer :: iwfn

   iwfn=29
   call open_file(iwfn,'molden.input','w')

   write(iwfn,'(A)',advance='yes')'[Molden Format]'
   write(iwfn,'(A)',advance='yes')'[Title]'

   !cccccccccccccccccccccccccccccccccccccccccccccccc
   ! print out atoms & coordinates                 c
   !cccccccccccccccccccccccccccccccccccccccccccccccc

   write(iwfn,'(A)',advance='yes')'[Atoms] AU'
   ! either coordinates are given in a.u.
   do i = 1,ncent
      atyp = toSymbol(at(i))
      z=at(i)
      ! a notion: for l>0, the ordering is p: x,y,z ; d: xx,yy,zz,xy,xz,yz ; f: xxx, yyy, zzz, xxy, xxz, yyx, yyz, xzz, yzz, xyz
      write (iwfn,'(a2,2i6,3E22.14)') &
         &   atyp,i,z,xyz(1,i),xyz(2,i),xyz(3,i)
   enddo
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Now print basis set data                                   c
   !                                                            c
   ! go through atoms and assign orbitals (contracted & prims)  c
   ! search for [GTO] statement to read basis set data          c
   ! for each set of atomic orbitals                            c
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   write(iwfn,'(A)',advance='yes') '[GTO]'
   !ccccccccccccccccccc
   ! go through atoms c
   !ccccccccccccccccccc
   icount=1
   do i=1,ncent
      write(iwfn,*) i,'0' ! I don't know what the zero is needed for
      !  now go trough nbfs located on atom i
      do j=basis%fila(1,i),basis%fila(2,i)
         call ang2chr(basis%lao(j),aang,skip)
         if(skip)then
            do k=1,basis%nprim(j)
               icount=icount+1
            enddo
         else
            write(iwfn,*) aang,basis%nprim(j),1.00
            do k=1,basis%nprim(j)
               write(iwfn,*) basis%alp(icount),basis%cont(icount)
               icount=icount+1
            enddo
         endif
      enddo
      write(iwfn,*)
   enddo

   nmomax=0
   do i=1,nmo
      if(evala(i).gt.thr.and.nmomax.eq.0)nmomax=i-1
   enddo
   if(nmomax.eq.0) nmomax=nmo

   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! print occupation number of orbitals and orbital energies c
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   write(iwfn,'(A)',advance='yes') '[MO]'
   ! Alpha MOs
   do i=1,nmomax
      write(iwfn,'(A)',advance='no') 'Sym= '
      write(iwfn,'(i5,a1)') i,'a (alpha)'
      write(iwfn,'(A)',advance='no') 'Ene= '
      write(iwfn,*) evala(i)
      write(iwfn,'(A)',advance='no') 'Spin= '
      write(iwfn,'(A)',advance='yes') 'Alpha'
      write(iwfn,'(A)',advance='no') 'Occup= '
      write(iwfn,'(F14.8)') occa(i)
      !now coefficients
      ! a notion: for l>0, the ordering is p: x,y,z ; d: xx,yy,zz,xy,xz,yz ; f: xxx, yyy, zzz, xxy, xxz, yyx, yyz, xzz, yzz, xyz
      do j=1,nbf
         write(iwfn,*) j,cmoa(j,i)
      enddo
   enddo
   ! Beta MOs
   do i=1,nmomax
      write(iwfn,'(A)',advance='no') 'Sym= '
      write(iwfn,'(i5,a1)') i,'a (beta)'
      write(iwfn,'(A)',advance='no') 'Ene= '
      write(iwfn,*) evalb(i)
      write(iwfn,'(A)',advance='no') 'Spin= '
      write(iwfn,'(A)',advance='yes') 'Beta'
      write(iwfn,'(A)',advance='no') 'Occup= '
      write(iwfn,'(F14.8)') occb(i)
      !now coefficients
      ! a notion: for l>0, the ordering is p: x,y,z ; d: xx,yy,zz,xy,xz,yz ; f: xxx, yyy, zzz, xxy, xxz, yyx, yyz, xzz, yzz, xyz
      do j=1,nbf
         write(iwfn,*) j,cmob(j,i)
      enddo
   enddo

   call close_file(iwfn)
end subroutine

! this routine gives character s,p,d etc. if angular momentum is given as input (i.e., 0,1,2 etc.)
subroutine ang2chr(iang,chr,skip)
   implicit none
   character*1, intent( out ) :: chr
   integer, intent( in ) :: iang
   logical, intent( inout ) :: skip

   skip=.false.
   select case(iang)
   case(1)
      chr='s'
   case(2)
      chr='p'
   case(3)
      chr='p'
      skip=.true.
   case(4)
      chr='p'
      skip=.true.
   case(5)
      chr='d'
   case(6)
      chr='d'
      skip=.true.
   case(7)
      chr='d'
      skip=.true.
   case(8)
      chr='d'
      skip=.true.
   case(9)
      chr='d'
      skip=.true.
   case(10)
      chr='d'
      skip=.true.
   end select
end subroutine




!ccccccccccccccccccccccccccccccccccccccccccccc
!
! specifies the length of an integer
!
!ccccccccccccccccccccccccccccccccccccccccccccc
subroutine lenint(iin,iout)
   implicit none
   integer, intent( in ) :: iin
   integer, intent( out) :: iout
   integer mdim,jdim
   jdim=0
   iout=0
   do
      iout=iout+1
      jdim=10*jdim+9
      mdim=iin-jdim
      if(mdim.le.0) exit
   enddo
   return
end subroutine


!     *****************************************************************
#if WITH_TBLITE
!> Module encapsulating the tblite Molden printer
module xtb_tblite_molden
   use mctc_env, only : wp
   use mctc_io_structure, only : structure_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_basis_type, only : basis_type, maxg
   use tblite_basis_cache, only : basis_cache
   implicit none
   private

   public :: print_molden

contains

   !> Writes the MOs, basis set, and geometry to a standard Molden input file
   subroutine print_molden(mol, wfn, basis, bcache, coeff_cart, thr)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Wavefunction structure data
      type(wavefunction_type), intent(in) :: wfn
      !> Basis set container
      class(basis_type), intent(in) :: basis
      !> Basis cache
      type(basis_cache), intent(in) :: bcache
      !> Cartesian MO coefficients, shape [nao_cart, nao, nspin]
      real(wp), intent(in) :: coeff_cart(:,:,:)
      !> Optional energy threshold to truncate the virtual MO output
      real(wp), intent(in), optional :: thr

      integer :: iwfn, iat, ish, iprim, imo, jao, ispin, nmomax
      integer :: l, ncart, icomp, iao_off
      integer :: perm(15)
      character(len=1) :: aang
      real(wp), allocatable :: prim_coeff(:)
      real(wp) :: thr_

      ! Energy threshold to truncate highly virtual states
      if (present(thr)) then
         thr_ = thr
      else
         thr_ = 2.0_wp
      end if

      ! Open Molden output file
      iwfn = 29
      open(unit=iwfn, file='molden.input', status='replace', action='write')

      write(iwfn,'(A)') '[Molden Format]'
      write(iwfn,'(A)') '[Title]'
      write(iwfn,'(A)') ' tblite output'

      ! Print atomic symbols, effective charge and coordinates
      write(iwfn,'(A)') '[Atoms] AU'

      do iat = 1, mol%nat
         ! Notation: symbol, atom, effective nuclear charge, coordinates
         write(iwfn,'(a2,2i6,3E22.14)') &
            & mol%sym(mol%id(iat)), iat, nint(wfn%n0at(iat)), &
            & mol%xyz(1,iat), mol%xyz(2,iat), mol%xyz(3,iat)
      end do

      ! Print basis set data 
      write(iwfn,'(A)') '[GTO]'

      allocate(prim_coeff(maxg))
      do iat = 1, mol%nat
         write(iwfn,'(i6,a)') iat, ' 0'
         do ish = 1, basis%nsh_at(iat)
            associate(p_cgto => basis%cgto(ish, mol%id(iat))%raw)
               select case(p_cgto%ang)
                  case(0); aang = 's'
                  case(1); aang = 'p'
                  case(2); aang = 'd'
                  case(3); aang = 'f'
                  case(4); aang = 'g'
                  case default
                     error stop "[Fatal] Molden writer only supports angular momenta up to g"
               end select

               ! Obtain the primtive contraction coefficients including 
               ! possible charge scaling and normalization factors
               call p_cgto%get_coeffs(bcache%cgto(ish, iat), prim_coeff)

               write(iwfn,'(a,i6,f8.2)') aang, p_cgto%nprim, 1.00_wp
               do iprim = 1, p_cgto%nprim
                  write(iwfn,'(2E22.14)') p_cgto%alpha(iprim), prim_coeff(iprim)
               end do
            end associate
         end do
         write(iwfn,*)
      end do
      deallocate(prim_coeff)

      ! Print occupation number, coefficients and orbital energies
      write(iwfn,'(A)') '[MO]'

      do ispin = 1, wfn%nspin
         ! Truncate the virtual MOs based on the given energy threshhold
         nmomax = 0
         do imo = 1, basis%nao
            if (wfn%emo(imo, ispin) > thr_ .and. nmomax == 0) nmomax = imo - 1
         end do
         if (nmomax == 0) nmomax = basis%nao

         do imo = 1, nmomax
            write(iwfn,'(A)', advance='no') 'Sym= '
            if (wfn%nspin == 2) then
               if (ispin == 1) then
                  write(iwfn,'(i5,a)') imo, 'a (alpha)'
               else
                  write(iwfn,'(i5,a)') imo, 'a (beta)'
               end if
            else
               write(iwfn,'(i5,a)') imo, 'a'
            end if

            write(iwfn,'(A)', advance='no') 'Ene= '
            write(iwfn,*) wfn%emo(imo, ispin)

            write(iwfn,'(A)', advance='no') 'Spin= '
            if (wfn%nspin == 2) then
               if (ispin == 1) then
                  write(iwfn,'(A)') 'Alpha'
               else
                  write(iwfn,'(A)') 'Beta'
               end if
            else
               write(iwfn,'(A)') 'Alpha'
            end if

            write(iwfn,'(A)', advance='no') 'Occup= '
            if (wfn%nspin == 1) then
               write(iwfn,'(F14.8)') wfn%focc(imo, 1) + wfn%focc(imo, 2)
            else
               write(iwfn,'(F14.8)') wfn%focc(imo, ispin)
            end if

            ! Write coefficients shell-by-shell in Molden cartesian order
            iao_off = 0
            jao = 0
            do iat = 1, mol%nat
               do ish = 1, basis%nsh_at(iat)
                  l = basis%cgto(ish, mol%id(iat))%raw%ang
                  ! Reorder the coefficients to match the Molden cartesian ordering
                  call get_molden_cart_perm(l, ncart, perm)

                  do icomp = 1, ncart
                     jao = jao + 1
                     write(iwfn,'(i6,1x,ES24.16)') jao, coeff_cart(iao_off + perm(icomp), imo, ispin)
                  end do

                  iao_off = iao_off + ncart
               end do
            end do

            if (jao /= basis%nao_cart) then
               error stop "[Fatal] Inconsistent AO count while writing Molden [MO] section"
            end if
         end do
      end do

      close(iwfn)
   end subroutine print_molden

   !> Return the cartesian Molden ordering permutation for one shell.
   !> The permutation maps tblite's internal coeff_cart shell order
   !> to the order expected in the Molden [MO] section.
   subroutine get_molden_cart_perm(l, ncart, perm)
      integer, intent(in)  :: l
      integer, intent(out) :: ncart
      integer, intent(out) :: perm(15)

      perm = 0

      select case(l)
      case(0)
         ! s
         ncart = 1
         perm(1) = 1

      case(1)
         ! tblite internal order follows spherical m = -1, 0, +1
         ! which corresponds to (py, pz, px) in the real-harmonic convention
         ! consistent with tblite_integral_trafo.
         ! Molden expects cartesian p order: px, py, pz
         ncart = 3
         perm(1:3) = [3, 1, 2]

      case(2)
         ! tblite internal cartesian d order already matches Molden 6D:
         ! xx, yy, zz, xy, xz, yz
         ncart = 6
         perm(1:6) = [1, 2, 3, 4, 5, 6]

      case(3)
         ! tblite internal 10F from ftrafo columns:
         ! xxx, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz
         !
         ! Molden 10F expects:
         ! xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
         ncart = 10
         perm(1:10) = [1, 2, 3, 6, 4, 5, 8, 9, 7, 10]

      case(4)
         ! tblite internal 15G from gtrafo columns already matches Molden 15G:
         ! xxxx, yyyy, zzzz, xxxy, xxxz, xyyy, yyyz, xzzz, yzzz,
         ! xxyy, xxzz, yyzz, xxyz, xyyz, xyzz
         ncart = 15
         perm(1:15) = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

      case default
         error stop "[Fatal] Molden writer only supports angular momenta up to g"
      end select
   end subroutine get_molden_cart_perm

end module xtb_tblite_molden
#endif


