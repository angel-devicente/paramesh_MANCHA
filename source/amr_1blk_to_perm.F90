!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****f* source/amr_1blk_to_perm
!! NAME
!!
!!   amr_1blk_to_perm
!!
!! SYNOPSIS
!!
!!   Call amr_1blk_to_perm (lcc,lfc,lec,lnc,lb,iopt,idest)
!!   Call amr_1blk_to_perm (logical,logical,logical,logical,
!!                          integer,integer,integer)
!!
!! ARGUMENTS
!!
!!   Logical, Intent(in) :: lcc  copies cell centered data if true
!!   Logical, Intent(in) :: lfc  copies cell face-centered data if true
!!   Logical, Intent(in) :: lec  copies cell edge-centered data if true
!!   Logical, Intent(in) :: lnc  copies cell corner data if true
!!   Integer, Intent(in) :: lb   block into which data is to be copied
!!   Integer, Intent(in) :: iopt data structure to be copied
!!   Integer, Intent(in) :: idest sets value for last dimension index
!!                                 in the 1-blk data arrays
!! 
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!   mpif.h
!!
!! USES
!!
!!    paramesh_dimensions
!!    physicaldata
!!    tree
!!    timings
!!    workspace
!!   
!! CALLS
!!
!!    No other Paramesh routines called.
!!
!! RETURNS
!!
!!    Nothing returned.
!!
!! DESCRIPTION
!!
!!   This routine copies data from the 1-block working arrays with guardcells
!!   to the permanent data arrays, which may or may not have permanent
!!   guardcells, depending on whether NO_PERMANENT_GUARDCELLS is defined 
!!   in physicaldata.fh.
!!
!! AUTHORS
!!
!!   Written :     Peter MacNeice          February 1999
!!  
!!*** 

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"

      subroutine amr_1blk_to_perm( lcc,lfc,lec,lnc,lb,iopt,idest)

!-----Use Statements
      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use workspace

      Implicit None

!-----Include statements
      include 'mpif.h'

!-----Input/Output Arguments
      Integer, Intent(in) :: lb,iopt,idest
      Logical, Intent(in) :: lcc,lfc,lec,lnc

!-----Local arrays and variables
      Integer :: iopt0, ivar, ivar_next
      Integer :: nguard0, nguard_work0
      Double Precision :: time1

!-----Begin Executable Code

      If (timing_mpi) Then
         time1 = mpi_wtime()
      End If

      nguard0 = nguard*(1-npgs)
      nguard_work0 = nguard_work*(1-npgs)

!------cell-centered data
      If (lcc) Then
         If (iopt == 1) Then

           If (ngcell_on_cc < nvar) Then
             Do ivar=1,ngcell_on_cc
               ivar_next = gcell_on_cc_pointer(ivar)
               unk(il_bnd:iu_bnd,jl_bnd:ju_bnd,              & 
                                           kl_bnd:ku_bnd,ivar_next,lb)           & 
               = unk1(il_bnd+nguard0:iu_bnd+nguard0,         & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,            & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,ivar_next,idest)
             End Do
           Else
             unk(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,:,lb)       & 
              = unk1(il_bnd+nguard0:iu_bnd+nguard0,                  & 
                       jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,          & 
                       kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,:,idest)
           End If  ! End If (ngcell_on_cc < nvar)

         ElseIf (iopt >= 2) Then

          iopt0 = iopt-1
           work(ilw:iuw,jlw:juw,klw:kuw,lb,iopt0)                      & 
             = work1(ilw+nguard_work0:iuw+nguard_work0,                & 
                     jlw+nguard_work0*k2d:juw+nguard_work0*k2d,        & 
                     klw+nguard_work0*k3d:kuw+nguard_work0*k3d,        & 
                     idest)

         End If  ! End If (iopt == 1)
      End If  ! End If (lcc)

!------cell face-centered data
      If (lfc) Then
!--------x-face
         If (ngcell_on_fc(1) < nfacevar) Then
           Do ivar=1,ngcell_on_fc(1)
             ivar_next = gcell_on_fc_pointer(1,ivar)
             facevarx(il_bnd:iu_bnd+1,                       & 
                      jl_bnd:ju_bnd,kl_bnd:ku_bnd,ivar_next,lb)                  & 
               = facevarx1(il_bnd+nguard0:iu_bnd+nguard0+1,  & 
                      jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,           & 
                      kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,ivar_next,idest)
           End Do
         Else
           facevarx(il_bnd:iu_bnd+1,                        & 
                             jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nfacevar,lb)           & 
             = facevarx1(il_bnd+nguard0:iu_bnd+nguard0+1,   & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,            & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,1:nfacevar,idest)
         End If  ! End If (ngcell_on_fc(1) < nfacevar)

!--------y-face
         If (ndim > 1) Then
         If (ngcell_on_fc(2) < nfacevar) Then
           Do ivar=1,ngcell_on_fc(2)
             ivar_next = gcell_on_fc_pointer(2,ivar)
             facevary(il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,       & 
                               kl_bnd:ku_bnd,ivar_next,lb)                       & 
             = facevary1(il_bnd+nguard0:iu_bnd+nguard0,      & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d+k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,ivar_next,idest)
           End Do
         Else
           facevary(il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,        & 
                               kl_bnd:ku_bnd,1:nfacevar,lb)                       & 
             = facevary1(il_bnd+nguard0:iu_bnd+nguard0,     & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d+k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,1:nfacevar,idest)
         End If  ! End If (ngcell_on_fc(2) < nfacevar)
         End If  ! End If (ndim > 1)

!--------z-face
         If (ndim == 3) Then
         If (ngcell_on_fc(3) < nfacevar) Then
           Do ivar=1,ngcell_on_fc(3)
             ivar_next = gcell_on_fc_pointer(3,ivar)
             facevarz(il_bnd:iu_bnd,jl_bnd:ju_bnd,           & 
                               kl_bnd:ku_bnd+k3d,ivar_next,lb)                   & 
             = facevarz1(il_bnd+nguard0:iu_bnd+nguard0,      & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,            & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d+k3d,ivar_next,idest)
           End Do
         Else
           facevarz(il_bnd:iu_bnd,jl_bnd:ju_bnd,            & 
                               kl_bnd:ku_bnd+k3d,1:nfacevar,lb)                   & 
             = facevarz1(il_bnd+nguard0:iu_bnd+nguard0,     & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,            & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d+k3d,1:nfacevar,idest)
         End If  ! End If (ngcell_on_fc(3) < nfacevar)

         End If  ! End If (ndim == 3)

      End If  ! End If (lfc)


!------cell edge-centered data

      If (ndim > 1) Then
       If (lec) Then

!--------x-edge
         If (ngcell_on_ec(1) < nvaredge) Then
           Do ivar=1,ngcell_on_ec(1)
             ivar_next = gcell_on_ec_pointer(1,ivar)
             unk_e_x(il_bnd:iu_bnd,                          & 
                     jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,ivar_next,lb)           & 
             = unk_e_x1(il_bnd+nguard0:iu_bnd+nguard0,       & 
                     jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,ivar_next,idest)
           End Do
         Else
           unk_e_x(il_bnd:iu_bnd,                           &  
                     jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,1:nvaredge,lb)           & 
             = unk_e_x1(il_bnd+nguard0:iu_bnd+nguard0,      & 
                     jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,1:nvaredge,idest)
         End If  ! End If (ngcell_on_ec(1) < nvaredge)

!--------y-edge
         If (ngcell_on_ec(2) < nvaredge) Then
           Do ivar=1,ngcell_on_ec(2)
             ivar_next = gcell_on_ec_pointer(2,ivar)
             unk_e_y(il_bnd:iu_bnd+1,                        & 
                     jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,ivar_next,lb)               & 
             = unk_e_y1(il_bnd+nguard0:iu_bnd+1+nguard0,     & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,            & 
                     kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,ivar_next,idest)

           End Do
         Else
           unk_e_y(il_bnd:iu_bnd+1,                         & 
                     jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,1:nvaredge,lb)               & 
             = unk_e_y1(il_bnd+nguard0:iu_bnd+1+nguard0,    & 
                     jl_bnd+nguard0*k2d:ju_bnd+nguard0*k2d,            & 
                     kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,1:nvaredge,idest)
         End If  ! End If (ngcell_on_ec(2) < nvaredge)

!--------z-edge
         If (ndim == 3) Then
         If (ngcell_on_ec(3) < nvaredge) Then
           Do ivar=1,ngcell_on_ec(3)
             ivar_next = gcell_on_ec_pointer(3,ivar)
             unk_e_z(il_bnd:iu_bnd+1,                        & 
                     jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,ivar_next,lb)               & 
             = unk_e_z1(il_bnd+nguard0:iu_bnd+1+nguard0,     & 
                     jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,ivar_next,idest)

           End Do
         Else
           unk_e_z(il_bnd:iu_bnd+1,                         & 
                     jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,1:nvaredge,lb)               & 
             = unk_e_z1(il_bnd+nguard0:iu_bnd+1+nguard0,    & 
                     jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+nguard0*k3d,1:nvaredge,idest)
         End If  ! End If (ngcell_on_ec(3) < nvaredge)
         End If  ! End If (ndim == 3)

        End If  ! End If (lec)
      End If  ! End If (ndim > 1)

!------cell corner data
      If (lnc) Then
         If (ngcell_on_nc < nvarcorn) Then
           Do ivar=1,ngcell_on_nc
             ivar_next = gcell_on_nc_pointer(ivar)
             unk_n(il_bnd:iu_bnd+1,                          & 
                     jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,ivar_next,lb)           & 
             = unk_n1(il_bnd+nguard0:iu_bnd+1+nguard0,       & 
                     jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,ivar_next,idest)
           End Do
         Else
           unk_n(il_bnd:iu_bnd+1,                           & 
                     jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,1:nvarcorn,lb)           & 
             = unk_n1(il_bnd+nguard0:iu_bnd+1+nguard0,      & 
                     jl_bnd+nguard0*k2d:ju_bnd+(nguard0+1)*k2d,        & 
                     kl_bnd+nguard0*k3d:ku_bnd+(nguard0+1)*k3d,1:nvarcorn,idest)
         End If  ! End If (ngcell_on_nc < nvarcorn)
      End If  ! End If (lnc)

      If (timing_mpi) Then
            timer_amr_1blk_to_perm(iopt) =                             & 
                               timer_amr_1blk_to_perm(iopt)            & 
                              + mpi_wtime() - time1
      End If

      Return
      End Subroutine amr_1blk_to_perm
