!=========================================================================
!
! Title         : shm_grid.f90
! Application   : Shape
! Copyright     : Tata Consultancy Services
! Author        : Niranjan Pedanekar, Chetan Malhotra
! Creation Date : 19/07/99
! Requirements  : Module "shm_glob_vars_mod"
! Description   : This file contains the subroutine which divides the 
!                 various inner surfaces of the reheating furnace into
!                 sections and stores their coordinates in structures
! Limitations   : 
! Dependencies  : Fortran 90
! Modifications :		WHO			WHEN		WHY
!				  Chetan Malhotra  20-8-99    Added gas
!				  Chetan Malhotra  9-3-2000   Made conducive to RHF model
!				  Chetan Malhotra  24-3-2000  Split wall into above and 
!											  below gas for gas 
!											  calculations
!=========================================================================
!
! Name              : shm_grid
! Description       : Divides the furnace geometry into arrays of surfaces
!                     and divides the surfaces into sections and stores 
!                     their coordinates into structures
! Parameters        : None
! Limitations		:
! Globals updated   : gv_roof_as,
!					  gv_wtf_as,gv_wtb_as,
!					  gv_wtfag_as,gv_wtbag_as,gv_wtfbg_as,gv_wtbbg_as,
!					  gv_st_as,
!					  gv_gtu_as,gv_gtd_as
!                     gv_floor_as,
!					  gv_wbf_as,gv_wbb_as,
!					  gv_wbfag_as,gv_wbbag_as,gv_wbfbg_as,gv_wbbbg_as,
!					  gv_sb_as,
!					  gv_gbu_as,gv_gbd_as
! Files updated     : None
! Calls             : None
!
!=========================================================================

subroutine shm_grid

use shm_glob_vars_mod
implicit none
integer gr_i_i, &                   ! Surface counter
        gr_j_i, &                   ! Division counter
        gr_k_i, &                   ! Flag
        gr_m_i                      ! Section counter
real*8 gr_sect_w_d, &               ! Section width
       gr_slope_d, &                ! Slope of surface in xy plane
       gr_theta_d                   ! Angle made by surface with x axis

!---------------------------------------------
! Segment Coordinates for Roof Active faces     
!---------------------------------------------

do gr_i_i=1,gv_NoOfTopSurf_i
  gv_roof_as(gr_i_i)%st_x_d=gv_rbpx_ad(gr_i_i)
  gv_roof_as(gr_i_i)%st_y_d=gv_rbpy_ad(gr_i_i)
  gv_roof_as(gr_i_i)%end_x_d=gv_rbpx_ad(gr_i_i+1)
  gv_roof_as(gr_i_i)%end_y_d=gv_rbpy_ad(gr_i_i+1)
  gr_j_i=0
  gr_k_i=1
  do while(gr_k_i.eq.1)
    gr_j_i=gr_j_i+1
    gr_sect_w_d=(((gv_roof_as(gr_i_i)%st_x_d- &
                   gv_roof_as(gr_i_i)%end_x_d)**2.0 + &
                  (gv_roof_as(gr_i_i)%st_y_d- &
                   gv_roof_as(gr_i_i)%end_y_d)**2.0)**0.5)/gr_j_i   
    if(gr_sect_w_d.le.gv_sect_max_d) then
      gv_roof_as(gr_i_i)%n_sections_i=gr_j_i
      gv_roof_as(gr_i_i)%sect_width_d=gr_sect_w_d
      gr_k_i=0
    endif
  end do
  if(dabs(gv_roof_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
    gr_slope_d=((gv_roof_as(gr_i_i)%end_y_d)- &
                (gv_roof_as(gr_i_i)%st_y_d))/ &
               ((gv_roof_as(gr_i_i)%end_x_d)- &
                (gv_roof_as(gr_i_i)%st_x_d))
  else
    gr_slope_d=0.0
  endif
  gr_theta_d=datan(gr_slope_d)
  do gr_m_i=1,gv_roof_as(gr_i_i)%n_sections_i
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_roof_as(gr_i_i)%st_x_d+ &
            (gr_m_i-1)*gv_roof_as(gr_i_i)%sect_width_d*dcos(gr_theta_d)
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_roof_as(gr_i_i)%st_y_d+ &
            (gr_m_i-1)*gv_roof_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_roof_as(gr_i_i)%st_x_d+ &
            gr_m_i*gv_roof_as(gr_i_i)%sect_width_d*dcos(gr_theta_d)
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_roof_as(gr_i_i)%st_y_d+ &
            gr_m_i*gv_roof_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_roof_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0
  end do
enddo

!------------------------------------------------------------
! Segment Coordinates for Wall Top Front / Back Active sfaces
!------------------------------------------------------------

do gr_i_i=1,gv_NoOfTopSurf_i
  gv_wtf_as(gr_i_i)%st_x_d=gv_roof_as(gr_i_i)%st_x_d
  gv_wtf_as(gr_i_i)%st_y_d=gv_roof_as(gr_i_i)%st_y_d
  gv_wtf_as(gr_i_i)%end_x_d=gv_roof_as(gr_i_i)%end_x_d
  gv_wtf_as(gr_i_i)%end_y_d=gv_roof_as(gr_i_i)%end_y_d
  gv_wtf_as(gr_i_i)%n_sections_i=gv_roof_as(gr_i_i)%n_sections_i
  gv_wtf_as(gr_i_i)%sect_width_d=0.0d0
  gv_wtb_as(gr_i_i)%st_x_d=gv_wtf_as(gr_i_i)%st_x_d
  gv_wtb_as(gr_i_i)%st_y_d=gv_wtf_as(gr_i_i)%st_y_d
  gv_wtb_as(gr_i_i)%end_x_d=gv_wtf_as(gr_i_i)%end_x_d
  gv_wtb_as(gr_i_i)%end_y_d=gv_wtf_as(gr_i_i)%end_y_d
  gv_wtb_as(gr_i_i)%n_sections_i=gv_wtf_as(gr_i_i)%n_sections_i
  gv_wtb_as(gr_i_i)%sect_width_d=gv_wtf_as(gr_i_i)%sect_width_d
  do gr_m_i=1,gv_wtf_as(gr_i_i)%n_sections_i
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)=gv_HearthYCoord_d + &
													gv_SlabThick_d
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)=gv_HearthYCoord_d + &
													gv_SlabThick_d
    gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_wtf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_wtb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d
  end do
enddo     

!-------------------------------------------------
! Segment Coordinates for Slab Top Active sfaces
!-------------------------------------------------

do gr_i_i=1,gv_NoOfTopSurf_i
  gv_st_as(gr_i_i)%st_x_d=gv_roof_as(gr_i_i)%st_x_d
  gv_st_as(gr_i_i)%st_y_d=gv_HearthYCoord_d + gv_SlabThick_d
  gv_st_as(gr_i_i)%end_x_d=gv_roof_as(gr_i_i)%end_x_d
  gv_st_as(gr_i_i)%end_y_d=gv_HearthYCoord_d + gv_SlabThick_d
  gv_st_as(gr_i_i)%n_sections_i=gv_roof_as(gr_i_i)%n_sections_i
  gv_st_as(gr_i_i)%sect_width_d=0.0d0
  do gr_m_i=1,gv_st_as(gr_i_i)%n_sections_i
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
			gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)=gv_HearthYCoord_d + &
												   gv_SlabThick_d
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)= &
		-(gv_fwidth_d/2.0d0)*(1.0d0 - gv_SlabCoverageFactor_d)
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
			gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)=gv_HearthYCoord_d + &
												   gv_SlabThick_d
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)= &
		-(gv_fwidth_d/2.0d0)*(1.0d0 - gv_SlabCoverageFactor_d)
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_st_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)=gv_HearthYCoord_d + &
												   gv_SlabThick_d
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)= &
		-(gv_fwidth_d/2.0d0)*(1.0d0 + gv_SlabCoverageFactor_d)
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_st_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)=gv_HearthYCoord_d + &
												   gv_SlabThick_d
    gv_st_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)= &
		-(gv_fwidth_d/2.0d0)*(1.0d0 + gv_SlabCoverageFactor_d)
  enddo
enddo

!-------------------------------------------------
! Segment Coordinates for Gas top-up Active sfaces    
!-------------------------------------------------

do gr_i_i=1,gv_NoOfTopSurf_i
  gv_gtu_as(gr_i_i)%st_x_d=gv_roof_as(gr_i_i)%st_x_d
  gv_gtu_as(gr_i_i)%st_y_d=gv_gtbpy_ad(gr_i_i)
  gv_gtu_as(gr_i_i)%end_x_d=gv_roof_as(gr_i_i)%end_x_d
  gv_gtu_as(gr_i_i)%end_y_d=gv_gtbpy_ad(gr_i_i+1)  
  gv_gtu_as(gr_i_i)%n_sections_i=gv_roof_as(gr_i_i)%n_sections_i
  gv_gtu_as(gr_i_i)%sect_width_d=&
				(((gv_gtu_as(gr_i_i)%st_x_d- &
                   gv_gtu_as(gr_i_i)%end_x_d)**2.0 + &
                  (gv_gtu_as(gr_i_i)%st_y_d- &
                   gv_gtu_as(gr_i_i)%end_y_d)**2.0)**0.5)/ &
				gv_gtu_as(gr_i_i)%n_sections_i
  if(dabs(gv_gtu_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
    gr_slope_d=((gv_gtu_as(gr_i_i)%end_y_d)- &
                (gv_gtu_as(gr_i_i)%st_y_d))/ &
               ((gv_gtu_as(gr_i_i)%end_x_d)- &
                (gv_gtu_as(gr_i_i)%st_x_d))
  else
    gr_slope_d=0.0
  endif
  gr_theta_d=datan(gr_slope_d)
  do gr_m_i=1,gv_gtu_as(gr_i_i)%n_sections_i
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
          gv_gtu_as(gr_i_i)%st_y_d+ &
          (gr_m_i-1)*gv_gtu_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
          gv_gtu_as(gr_i_i)%st_y_d+ &
          gr_m_i*gv_gtu_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d
  end do  !do for gr_m_i no. of sections
enddo   

!---------------------------------------------------
! Segment Coordinates for Gas top-down Active sfaces    
!---------------------------------------------------

do gr_i_i=1,gv_NoOfTopSurf_i
  gv_gtd_as(gr_i_i)%st_x_d=gv_gtu_as(gr_i_i)%st_x_d
  gv_gtd_as(gr_i_i)%st_y_d=gv_gtu_as(gr_i_i)%st_y_d
  gv_gtd_as(gr_i_i)%end_x_d=gv_gtu_as(gr_i_i)%end_x_d
  gv_gtd_as(gr_i_i)%end_y_d=gv_gtu_as(gr_i_i)%end_y_d
  gv_gtd_as(gr_i_i)%n_sections_i=gv_gtu_as(gr_i_i)%n_sections_i
  gv_gtd_as(gr_i_i)%sect_width_d=gv_gtu_as(gr_i_i)%sect_width_d
  do gr_m_i=1,gv_gtd_as(gr_i_i)%n_sections_i
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
		  gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
          gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)= &
		  gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
		  gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
		  gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)= &
		  gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
          gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
          gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)= &
          gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
          gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
          gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_gtd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)= &
		  gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)
  end do  !do for gr_m_i no. of sections
enddo   

!------------------------------------------------------------
! Segment Coordinates for Wall Top Front/Back Above/Below Gas 
! Active sfaces
!------------------------------------------------------------

do gr_i_i=1,gv_NoOfTopSurf_i

  gv_wtfag_as(gr_i_i)%st_x_d=gv_roof_as(gr_i_i)%st_x_d
  gv_wtfag_as(gr_i_i)%st_y_d=gv_roof_as(gr_i_i)%st_y_d
  gv_wtfag_as(gr_i_i)%end_x_d=gv_roof_as(gr_i_i)%end_x_d
  gv_wtfag_as(gr_i_i)%end_y_d=gv_roof_as(gr_i_i)%end_y_d
  gv_wtfag_as(gr_i_i)%n_sections_i=gv_roof_as(gr_i_i)%n_sections_i
  gv_wtfag_as(gr_i_i)%sect_width_d=0.0d0

  gv_wtbag_as(gr_i_i)%st_x_d=gv_wtf_as(gr_i_i)%st_x_d
  gv_wtbag_as(gr_i_i)%st_y_d=gv_wtf_as(gr_i_i)%st_y_d
  gv_wtbag_as(gr_i_i)%end_x_d=gv_wtf_as(gr_i_i)%end_x_d
  gv_wtbag_as(gr_i_i)%end_y_d=gv_wtf_as(gr_i_i)%end_y_d
  gv_wtbag_as(gr_i_i)%n_sections_i=gv_wtf_as(gr_i_i)%n_sections_i
  gv_wtbag_as(gr_i_i)%sect_width_d=gv_wtf_as(gr_i_i)%sect_width_d

  gv_wtfbg_as(gr_i_i)%st_x_d=gv_wtf_as(gr_i_i)%st_x_d
  gv_wtfbg_as(gr_i_i)%st_y_d=gv_wtf_as(gr_i_i)%st_y_d
  gv_wtfbg_as(gr_i_i)%end_x_d=gv_wtf_as(gr_i_i)%end_x_d
  gv_wtfbg_as(gr_i_i)%end_y_d=gv_wtf_as(gr_i_i)%end_y_d
  gv_wtfbg_as(gr_i_i)%n_sections_i=gv_wtf_as(gr_i_i)%n_sections_i
  gv_wtfbg_as(gr_i_i)%sect_width_d=gv_wtf_as(gr_i_i)%sect_width_d

  gv_wtbbg_as(gr_i_i)%st_x_d=gv_wtf_as(gr_i_i)%st_x_d
  gv_wtbbg_as(gr_i_i)%st_y_d=gv_wtf_as(gr_i_i)%st_y_d
  gv_wtbbg_as(gr_i_i)%end_x_d=gv_wtf_as(gr_i_i)%end_x_d
  gv_wtbbg_as(gr_i_i)%end_y_d=gv_wtf_as(gr_i_i)%end_y_d
  gv_wtbbg_as(gr_i_i)%n_sections_i=gv_wtf_as(gr_i_i)%n_sections_i
  gv_wtbbg_as(gr_i_i)%sect_width_d=gv_wtf_as(gr_i_i)%sect_width_d

  do gr_m_i=1,gv_wtf_as(gr_i_i)%n_sections_i


    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_roof_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_gtu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0


    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
			gv_HearthYCoord_d + gv_SlabThick_d
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
			gv_HearthYCoord_d + gv_SlabThick_d
    gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0


    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_wtfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_wtbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d


    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_wtfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_wtbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d


  end do
enddo     

!----------------------------------------------
! Segment Coordinates for Floor Active sfaces    
!----------------------------------------------

do gr_i_i=1,gv_NoOfBotSurf_i
  gv_floor_as(gr_i_i)%st_x_d=gv_fbpx_ad(gr_i_i)
  gv_floor_as(gr_i_i)%st_y_d=gv_fbpy_ad(gr_i_i)
  gv_floor_as(gr_i_i)%end_x_d=gv_fbpx_ad(gr_i_i+1)
  gv_floor_as(gr_i_i)%end_y_d=gv_fbpy_ad(gr_i_i+1)  
  gr_j_i=0
  gr_k_i=1
  do while(gr_k_i.eq.1)
    gr_j_i=gr_j_i+1
    gr_sect_w_d=(((gv_floor_as(gr_i_i)%st_x_d-gv_floor_as(gr_i_i)% &
                   end_x_d)**2.0 + &
                  (gv_floor_as(gr_i_i)%st_y_d-gv_floor_as(gr_i_i)% &
                   end_y_d)**2.0)**0.5)/gr_j_i  
    if (gr_sect_w_d.le.gv_sect_max_d) then
      gv_floor_as(gr_i_i)%n_sections_i=gr_j_i
      gv_floor_as(gr_i_i)%sect_width_d=gr_sect_w_d
      gr_k_i=0
    endif
  end do
  if(dabs(gv_floor_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
    gr_slope_d=((gv_floor_as(gr_i_i)%end_y_d)- &
                (gv_floor_as(gr_i_i)%st_y_d))/ &
               ((gv_floor_as(gr_i_i)%end_x_d)- &
                (gv_floor_as(gr_i_i)%st_x_d))
  else
    gr_slope_d=0.0
  endif
  gr_theta_d=datan(gr_slope_d)
  do gr_m_i=1,gv_floor_as(gr_i_i)%n_sections_i
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
          gv_floor_as(gr_i_i)%st_x_d+ &
          (gr_m_i-1)*gv_floor_as(gr_i_i)%sect_width_d*dcos(gr_theta_d) 
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
          gv_floor_as(gr_i_i)%st_y_d+ &
          (gr_m_i-1)*gv_floor_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
          gv_floor_as(gr_i_i)%st_x_d+ &
          gr_m_i*gv_floor_as(gr_i_i)%sect_width_d*dcos(gr_theta_d)
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
          gv_floor_as(gr_i_i)%st_y_d+ &
          gr_m_i*gv_floor_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_floor_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d
  end do  !do for gr_m_i no. of sections
enddo   

!---------------------------------------------------------------
! Segment Coordinates for Wall Bottom Front / Back Active sfaces
!---------------------------------------------------------------

do gr_i_i=1,gv_NoOfBotSurf_i
  gv_wbf_as(gr_i_i)%st_x_d=gv_floor_as(gr_i_i)%st_x_d
  gv_wbf_as(gr_i_i)%st_y_d=gv_floor_as(gr_i_i)%st_y_d
  gv_wbf_as(gr_i_i)%end_x_d=gv_floor_as(gr_i_i)%end_x_d
  gv_wbf_as(gr_i_i)%end_y_d=gv_floor_as(gr_i_i)%end_y_d
  gv_wbf_as(gr_i_i)%n_sections_i=gv_floor_as(gr_i_i)%n_sections_i
  gv_wbf_as(gr_i_i)%sect_width_d=0.0d0
  gv_wbb_as(gr_i_i)%st_x_d=gv_wbf_as(gr_i_i)%st_x_d
  gv_wbb_as(gr_i_i)%st_y_d=gv_wbf_as(gr_i_i)%st_y_d
  gv_wbb_as(gr_i_i)%end_x_d=gv_wbf_as(gr_i_i)%end_x_d
  gv_wbb_as(gr_i_i)%end_y_d=gv_wbf_as(gr_i_i)%end_y_d
  gv_wbb_as(gr_i_i)%n_sections_i=gv_wbf_as(gr_i_i)%n_sections_i
  gv_wbb_as(gr_i_i)%sect_width_d=gv_wbf_as(gr_i_i)%sect_width_d
  do gr_m_i=1,gv_wbf_as(gr_i_i)%n_sections_i
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)=gv_HearthYCoord_d
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)=gv_HearthYCoord_d
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_wbf_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_wbb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d
  enddo
enddo     
 
!----------------------------------------------------
! Segment Coordinates for Slab Bottom Active sfaces
!----------------------------------------------------

do gr_i_i=1,gv_NoOfBotSurf_i
  gv_sb_as(gr_i_i)%st_x_d=gv_floor_as(gr_i_i)%st_x_d
  gv_sb_as(gr_i_i)%st_y_d=gv_HearthYCoord_d
  gv_sb_as(gr_i_i)%end_x_d=gv_floor_as(gr_i_i)%end_x_d
  gv_sb_as(gr_i_i)%end_y_d=gv_HearthYCoord_d
  gv_sb_as(gr_i_i)%n_sections_i=gv_floor_as(gr_i_i)%n_sections_i
  gv_sb_as(gr_i_i)%sect_width_d=0.0d0
  do gr_m_i=1,gv_sb_as(gr_i_i)%n_sections_i
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)=gv_HearthYCoord_d
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)= &
		-(gv_fwidth_d/2.0d0)*(1.0d0 - gv_SlabCoverageFactor_d)
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)=gv_HearthYCoord_d
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)= &
		-(gv_fwidth_d/2.0d0)*(1.0d0 + gv_SlabCoverageFactor_d)
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)=gv_HearthYCoord_d
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)= &
		-(gv_fwidth_d/2.0d0)*(1.0d0 + gv_SlabCoverageFactor_d)
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_sb_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)=gv_HearthYCoord_d
    gv_sb_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)= &
		-(gv_fwidth_d/2.0d0)*(1.0d0 - gv_SlabCoverageFactor_d)
  enddo
enddo

!----------------------------------------------------
! Segment Coordinates for Gas bottom-up Active sfaces    
!----------------------------------------------------

do gr_i_i=1,gv_NoOfBotSurf_i
  gv_gbu_as(gr_i_i)%st_x_d=gv_floor_as(gr_i_i)%st_x_d
  gv_gbu_as(gr_i_i)%st_y_d=gv_gbbpy_ad(gr_i_i)
  gv_gbu_as(gr_i_i)%end_x_d=gv_floor_as(gr_i_i)%end_x_d
  gv_gbu_as(gr_i_i)%end_y_d=gv_gbbpy_ad(gr_i_i+1)  
  gv_gbu_as(gr_i_i)%n_sections_i=gv_floor_as(gr_i_i)%n_sections_i
  gv_gbu_as(gr_i_i)%sect_width_d=&
				(((gv_gbu_as(gr_i_i)%st_x_d- &
                   gv_gbu_as(gr_i_i)%end_x_d)**2.0 + &
                  (gv_gbu_as(gr_i_i)%st_y_d- &
                   gv_gbu_as(gr_i_i)%end_y_d)**2.0)**0.5)/ &
				gv_gbu_as(gr_i_i)%n_sections_i
  if(dabs(gv_gbu_as(gr_i_i)%sect_width_d).gt.1.0d-7) then
    gr_slope_d=((gv_gbu_as(gr_i_i)%end_y_d)- &
                (gv_gbu_as(gr_i_i)%st_y_d))/ &
               ((gv_gbu_as(gr_i_i)%end_x_d)- &
                (gv_gbu_as(gr_i_i)%st_x_d))
  else
    gr_slope_d=0.0
  endif
  gr_theta_d=datan(gr_slope_d)
  do gr_m_i=1,gv_gbu_as(gr_i_i)%n_sections_i
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
		  gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
          gv_gbu_as(gr_i_i)%st_y_d+ &
          (gr_m_i-1)*gv_gbu_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
	      gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
          gv_gbu_as(gr_i_i)%st_y_d+ &
          gr_m_i*gv_gbu_as(gr_i_i)%sect_width_d*dsin(gr_theta_d)
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d
  end do  !do for gr_m_i no. of sections
enddo   

!------------------------------------------------------
! Segment Coordinates for Gas bottom-down Active sfaces    
!------------------------------------------------------

do gr_i_i=1,gv_NoOfBotSurf_i
  gv_gbd_as(gr_i_i)%st_x_d=gv_gbu_as(gr_i_i)%st_x_d
  gv_gbd_as(gr_i_i)%st_y_d=gv_gbu_as(gr_i_i)%st_y_d
  gv_gbd_as(gr_i_i)%end_x_d=gv_gbu_as(gr_i_i)%end_x_d
  gv_gbd_as(gr_i_i)%end_y_d=gv_gbu_as(gr_i_i)%end_y_d
  gv_gbd_as(gr_i_i)%n_sections_i=gv_gbu_as(gr_i_i)%n_sections_i
  gv_gbd_as(gr_i_i)%sect_width_d=gv_gbu_as(gr_i_i)%sect_width_d
  do gr_m_i=1,gv_gbd_as(gr_i_i)%n_sections_i
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
		  gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
          gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)= &
		  gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
		  gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
		  gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)= &
		  gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
          gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
          gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)= &
          gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
          gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
          gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_gbd_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)= &
		  gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)
  end do  !do for gr_m_i no. of sections
enddo   

!---------------------------------------------------------------
! Segment Coordinates for Wall Bottom Front/Back Above/Below Gas 
! Active sfaces
!---------------------------------------------------------------

do gr_i_i=1,gv_NoOfBotSurf_i

  gv_wbfag_as(gr_i_i)%st_x_d=gv_floor_as(gr_i_i)%st_x_d
  gv_wbfag_as(gr_i_i)%st_y_d=gv_floor_as(gr_i_i)%st_y_d
  gv_wbfag_as(gr_i_i)%end_x_d=gv_floor_as(gr_i_i)%end_x_d
  gv_wbfag_as(gr_i_i)%end_y_d=gv_floor_as(gr_i_i)%end_y_d
  gv_wbfag_as(gr_i_i)%n_sections_i=gv_floor_as(gr_i_i)%n_sections_i
  gv_wbfag_as(gr_i_i)%sect_width_d=0.0d0

  gv_wbbag_as(gr_i_i)%st_x_d=gv_wbf_as(gr_i_i)%st_x_d
  gv_wbbag_as(gr_i_i)%st_y_d=gv_wbf_as(gr_i_i)%st_y_d
  gv_wbbag_as(gr_i_i)%end_x_d=gv_wbf_as(gr_i_i)%end_x_d
  gv_wbbag_as(gr_i_i)%end_y_d=gv_wbf_as(gr_i_i)%end_y_d
  gv_wbbag_as(gr_i_i)%n_sections_i=gv_wbf_as(gr_i_i)%n_sections_i
  gv_wbbag_as(gr_i_i)%sect_width_d=gv_wbf_as(gr_i_i)%sect_width_d

  gv_wbfbg_as(gr_i_i)%st_x_d=gv_wbf_as(gr_i_i)%st_x_d
  gv_wbfbg_as(gr_i_i)%st_y_d=gv_wbf_as(gr_i_i)%st_y_d
  gv_wbfbg_as(gr_i_i)%end_x_d=gv_wbf_as(gr_i_i)%end_x_d
  gv_wbfbg_as(gr_i_i)%end_y_d=gv_wbf_as(gr_i_i)%end_y_d
  gv_wbfbg_as(gr_i_i)%n_sections_i=gv_wbf_as(gr_i_i)%n_sections_i
  gv_wbfbg_as(gr_i_i)%sect_width_d=gv_wbf_as(gr_i_i)%sect_width_d

  gv_wbbbg_as(gr_i_i)%st_x_d=gv_wbf_as(gr_i_i)%st_x_d
  gv_wbbbg_as(gr_i_i)%st_y_d=gv_wbf_as(gr_i_i)%st_y_d
  gv_wbbbg_as(gr_i_i)%end_x_d=gv_wbf_as(gr_i_i)%end_x_d
  gv_wbbbg_as(gr_i_i)%end_y_d=gv_wbf_as(gr_i_i)%end_y_d
  gv_wbbbg_as(gr_i_i)%n_sections_i=gv_wbf_as(gr_i_i)%n_sections_i
  gv_wbbbg_as(gr_i_i)%sect_width_d=gv_wbf_as(gr_i_i)%sect_width_d

  do gr_m_i=1,gv_wbf_as(gr_i_i)%n_sections_i


    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= gv_HearthYCoord_d
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_floor_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= gv_HearthYCoord_d
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_gbu_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0


    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=0.0
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=0.0
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
			gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=0.0
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
			gv_floor_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=0.0


    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_wbfag_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_wbbag_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d


    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)= &
            gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(1)
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)= &
            gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(1)
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(1)=-gv_fwidth_d
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)= &
            gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)= &
            gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(2)=-gv_fwidth_d
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)= &
            gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(3)
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)= &
            gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(3)
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(3)=-gv_fwidth_d
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(4)= &
            gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%x_coord_ad(2)
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(4)= &
            gv_wbfbg_as(gr_i_i)%sect_as(gr_m_i)%y_coord_ad(2)
    gv_wbbbg_as(gr_i_i)%sect_as(gr_m_i)%z_coord_ad(4)=-gv_fwidth_d


  end do
enddo     

! =========================================================
!       Number of sections in top, bottom, gas top, 
!				  gas bottom and slab
! =========================================================

gv_ntop_i=0
do gr_i_i=1,gv_NoOfTopSurf_i
  gv_ntop_i=gv_ntop_i+gv_roof_as(gr_i_i)%n_sections_i
enddo
gv_nbot_i=0
do gr_i_i=1,gv_NoOfBotSurf_i
  gv_nbot_i=gv_nbot_i+gv_floor_as(gr_i_i)%n_sections_i
enddo

! =========================================================
!       Allocate memory for sectionwise shape factors
! =========================================================

allocate(gv_phi_st_r_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_r_st_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_wtf_r_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_r_wtf_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_wtb_r_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_r_wtb_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_wtf_st_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_st_wtf_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_wtb_st_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_st_wtb_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_wtf_wtb_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_wtb_wtf_ad(gv_ntop_i,gv_ntop_i), &
         gv_phi_r_r_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_gtu_r_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_r_gtu_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_gtd_r_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_r_gtd_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_gtd_st_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_st_gtd_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_gtu_wtfag_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_wtfag_gtu_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_gtu_wtbag_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_wtbag_gtu_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_gtd_wtfbg_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_wtfbg_gtd_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_gtd_wtbbg_ad(gv_ntop_i,gv_ntop_i), &
		 gv_phi_wtbbg_gtd_ad(gv_ntop_i,gv_ntop_i))

allocate(gv_phi_sb_f_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_f_sb_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_wbf_f_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_f_wbf_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_wbb_f_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_f_wbb_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_wbf_sb_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_sb_wbf_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_wbb_sb_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_sb_wbb_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_wbf_wbb_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_wbb_wbf_ad(gv_nbot_i,gv_nbot_i), &
         gv_phi_f_f_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_gbu_f_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_f_gbu_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_gbd_f_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_f_gbd_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_gbu_sb_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_sb_gbu_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_gbu_wbfag_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_wbfag_gbu_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_gbu_wbbag_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_wbbag_gbu_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_gbd_wbfbg_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_wbfbg_gbd_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_gbd_wbbbg_ad(gv_nbot_i,gv_nbot_i), &
		 gv_phi_wbbbg_gbd_ad(gv_nbot_i,gv_nbot_i))

return
end subroutine shm_grid
