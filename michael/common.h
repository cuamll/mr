!-----------------------------------------------------------!
!-----------------------------------------------------------!
!                                                           !
!          Global variables common to both the main         !
!                   and to all subroutines                  !
!                                                           !
!-----------------------------------------------------------!
!-----------------------------------------------------------!

implicit none
real*8 pi,twopi,pibytwo,rootpi,gamma,root2
  parameter (pi=3.141592653589793)
  parameter (twopi=6.283185307179590)
  parameter (pibytwo=1.570796326794897)
  parameter (rootpi=1.77245385091)
  parameter (gamma=0.5772156649)
  parameter (root2=1.41421356237309504880)
  integer max_sides,max_sites,max_sweeps
  parameter (max_sides=200)
  parameter (max_sites=max_sides*max_sides)
  integer pos(max_sides),neg(max_sides)
  integer side,sites,accept,acceptCE,Tsteps,Esteps,sidesteps,sidemin,sidemax,sideincr,acceptg,accepth,ratio,ratiorot,ratiowind
  integer therm_sweeps,sweeps,quench_sweeps,bin_size,bin_size_curr,num_bins,measurements,ibin,ibin_curr,twist,MCT,glob,globy,printU

  real*8  theta(max_sides,max_sides),v(max_sides,max_sides)
  real*8  top_x(max_sides,max_sides),top_y(max_sides,max_sides)
  real*8  n_x(max_sides,max_sides),n_y(max_sides,max_sides)
  real*8  Ehat_x(max_sides,max_sides),Ehat_y(max_sides,max_sides),wind_x(max_sides,max_sides),wind_y(max_sides,max_sides)
  real*8  phi(max_sides,max_sides),Ehatsum_x,Ehatsum_y,nx,ny,Jx,Jy
  real*8  SW(max_sides,max_sides),topDfree_x(max_sides,max_sides),topDfree_y(max_sides,max_sides)
  real*8  Eg_x,Eg_y,quench_fluct,delta_fluct,delta_phi,delta_delta_phi
  real*8  volume,length,T,Tstore,Tquench,beta,Tmin,Tmax,Tlow,Tincr,vort,abs_vort,vort_probe,Uself_factor,abel
  real*8  E,Emin,Emax,Elow,Eincr,Eapp,Efluct,abswinding,dx,dy,Ebar_x,Ebar_y,Ebardip_x,Ebardip_y,Ebarwind_x,Ebarwind_y,hel_real
  real*8  sum_energy,sum_enersq,sum_ghost_energy,sum_discrete_ghost_energy,sum_enerf
  real*8  sum_Utot,sum_Utotsq,sum_Uvort,sum_Uvortsq,sum_Uharm,sum_Uharmsq
  real*8  sum_accept,sum_accepth,sum_acceptg,sum_hop_proposal,sum_spins,sum_spinsq,sum_mx,sum_my,sum_spinf,sum_mag
  real*8  sum_vort,sum_vortsq,sum_abs_vort,sum_curr_x,sum_curr_xsq,sum_curr_y,sum_curr_ysq,sum_pr,sum_prsq,sum_core,sum_winding
  real*8  sum_nx,sum_ny,sum_nxsq,sum_nysq,sum_dx,sum_dy,sum_dxsq,sum_dysq
  real*8  sum_Ehatsum_x,sum_Ehatsum_xsq,sum_Ehatsum_y,sum_Ehatsum_ysq,sum_e,sum_s
  real*8  sum_Ebar_x,sum_Ebar_y,sum_Ebar_xsq,sum_Ebar_ysq
  real*8  sum_Ebardip_x,sum_Ebardip_y,sum_Ebardip_xsq,sum_Ebardip_ysq
  real*8  sum_Ebarwind_x,sum_Ebarwind_y,sum_Ebarwind_xsq,sum_Ebarwind_ysq,sum_hel_real
  real*8  sum_senergy,sum_senersq,sum_cenergy,sum_cenersq,sum_fieldsq
  real*8  sum_delta_x,sum_delta_y,av_delta_x,av_delta_y,sum_sizex,sum_sizey,av_sizex,av_sizey
  real*8  sum_factorA,sum_factorB,av_factorA,av_factorB
  real*8  sum_chi_p,sum_chi_w,sum_chiCG,sum_chiCG_p,sum_chiCG_w,sum_chiCG_cross,sum_e_p,sum_e_w
  real*8  bin_energy,bin_enersq,bin_ghost_energy,bin_discrete_ghost_energy,bin_spins,bin_spinsq,bin_mx,bin_my,bin_spinf,bin_enerf
  real*8  bin_vort,bin_vortsq,bin_curr_x,bin_curr_xsq,bin_curr_y,bin_curr_ysq,bin_pr,bin_prsq,bin_core,bin_winding
  real*8  bin_Ebar_x,bin_Ebar_xsq,bin_Ebar_y,bin_Ebar_ysq,bin_Ebardip_x,bin_Ebardip_xsq,bin_Ebardip_y,bin_Ebardip_ysq
  real*8  bin_Ebarwind_x,bin_Ebarwind_xsq,bin_Ebarwind_y,bin_Ebarwind_ysq,bin_totsus,bin_dipsus,bin_windsus
  real*8  bin_Ebarsq,bin_Ebardipsq,bin_Ebarwindsq
  real*8  bin_nx,bin_ny,bin_nxsq,bin_nysq,bin_e,bin_s
  real*8  avbin_spheat,avbin_spheatsq,avbin_susp,avbin_suspsq
  real*8  avbin_enersq,avbin_ghost_enersq,avbin_discrete_ghost_enersq,avbin_spinsq,avbin_spinf,avbin_enerf
  real*8  avbin_curr_xsq,avbin_curr_ysq,avbin_prsq,avbin_usq,avbin_vortsq,avbin_windingsq,avbin_nxsq,avbin_nysq
  real*8  avbin_totsus,avbin_totsussq,avbin_dipsus,avbin_dipsussq,avbin_windsus,avbin_windsussq
  real*8  avbin_Ebarf,avbin_Ebardipf,avbin_Ebarwindf
  real*8  av_energy,av_enersq,av_ghost_energy,av_discrete_ghost_energy,av_spins,av_spinsq,specific_heat
  real*8  av_accept,av_accepth,av_acceptg
  real*8  av_vort,av_vortsq,av_abs_vort,av_curr_x,av_curr_xsq,av_curr_y,av_curr_ysq,av_curr_x_old,av_pr,av_prsq,av_core,av_winding
  real*8  av_nx,av_ny,av_nxsq,av_nysq,av_dx,av_dy,av_dxsq,av_dysq,av_Ehatsum_x,av_Ehatsum_xsq,av_Ehatsum_y,av_Ehatsum_ysq
  real*8  av_Ebar_x,av_Ebar_y,av_Ebar_xsq,av_Ebar_ysq,av_Ebardip_x,av_Ebardip_y,av_Ebardip_xsq,av_Ebardip_ysq
  real*8  av_Ebarwind_x,av_Ebarwind_y,av_Ebarwind_xsq,av_Ebarwind_ysq
  real*8  susceptibility,av_mx,av_my,av_spinf,av_enerf,ul,vl,av_mag
  
  character*20 spinin,spinout,vortout,xtopout,xtopin,ytopout,ytopin
  character*50 outfile1,outfile2,outfile3,outfile4,outfile5,plot_spin_data,plot_vort_data
  character*50 top_snap_shot,spin_snap_shot
  
  common  T,Tstore,Tquench,beta,Tmin,Tmax,Tlow,vort,abs_vort,vort_probe,Uself_factor,abel
  common  E,Emin,Emax,Elow,Eincr,Eapp,Efluct,abswinding,dx,dy,Ebar_x,Ebar_y,Ebardip_x,Ebardip_y,Ebarwind_x,Ebarwind_y
  common  sum_accept,sum_accepth,sum_acceptg,sum_hop_proposal
  common  sum_energy,sum_enersq,sum_ghost_energy,sum_discrete_ghost_energy,sum_enerf
  common  sum_Utot,sum_Utotsq,sum_Uvort,sum_Uvortsq,sum_Uharm,sum_Uharmsq
  common  sum_spins,sum_spinsq,sum_mx,sum_my,sum_spinf,sum_mag
  common  sum_vort,sum_vortsq,sum_abs_vort,sum_curr_x,sum_curr_xsq,sum_curr_y,sum_curr_ysq,sum_pr,sum_prsq,sum_core,sum_winding
  common  sum_nx,sum_ny,sum_nxsq,sum_nysq,sum_dx,sum_dy,sum_dxsq,sum_dysq
  common  sum_Ehatsum_x,sum_Ehatsum_xsq,sum_Ehatsum_y,sum_Ehatsum_ysq,sum_e,sum_s
  common  sum_Ebar_x,sum_Ebar_y,sum_Ebar_xsq,sum_Ebar_ysq
  common  sum_Ebardip_x,sum_Ebardip_y,sum_Ebardip_xsq,sum_Ebardip_ysq
  common  sum_Ebarwind_x,sum_Ebarwind_y,sum_Ebarwind_xsq,sum_Ebarwind_ysq,sum_hel_real
  common  sum_senergy,sum_senersq,sum_cenergy,sum_cenersq,sum_fieldsq
  common  sum_delta_x,sum_delta_y,av_delta_x,av_delta_y,sum_sizex,sum_sizey,av_sizex,av_sizey
  common  sum_factorA,sum_factorB,av_factorA,av_factorB
  common  sum_chi_p,sum_chi_w,sum_chiCG,sum_chiCG_p,sum_chiCG_w,sum_chiCG_cross,sum_e_p,sum_e_w
  common  bin_energy,bin_enersq,bin_ghost_energy,bin_discrete_ghost_energy,bin_spins,bin_spinsq,bin_mx,bin_my,bin_spinf,bin_enerf
  common  bin_vort,bin_vortsq,bin_curr_x,bin_curr_xsq,bin_curr_y,bin_curr_ysq,bin_pr,bin_prsq,bin_core,bin_winding
  common  bin_Ebar_x,bin_Ebar_xsq,bin_Ebar_y,bin_Ebar_ysq,bin_Ebardip_x,bin_Ebardip_xsq,bin_Ebardip_y,bin_Ebardip_ysq
  common  bin_Ebarwind_x,bin_Ebarwind_xsq,bin_Ebarwind_y,bin_Ebarwind_ysq,bin_totsus,bin_dipsus,bin_windsus
  common  bin_Ebarsq,bin_Ebardipsq,bin_Ebarwindsq
  common  bin_nx,bin_ny,bin_nxsq,bin_nysq,bin_e,bin_s
  common  avbin_spheat,avbin_spheatsq,avbin_susp,avbin_suspsq
  common  avbin_enersq,avbin_ghost_enersq,avbin_discrete_ghost_enersq,avbin_spinsq,avbin_spinf,avbin_enerf
  common  avbin_curr_xsq,avbin_curr_ysq,avbin_prsq,avbin_usq,avbin_vortsq,avbin_windingsq,avbin_nxsq,avbin_nysq
  common  avbin_totsus,avbin_totsussq,avbin_dipsus,avbin_dipsussq,avbin_windsus,avbin_windsussq
  common  avbin_Ebarf,avbin_Ebardipf,avbin_Ebarwindf
  common  av_energy,av_enersq,av_ghost_energy,av_discrete_ghost_energy,av_spins,av_spinsq,specific_heat
  common  av_accept,av_accepth,av_acceptg
  common  av_vort,av_vortsq,av_abs_vort,av_curr_x,av_curr_xsq,av_curr_y,av_curr_ysq,av_curr_x_old,av_pr,av_prsq,av_core,av_winding
  common  av_nx,av_ny,av_nxsq,av_nysq,av_dx,av_dy,av_dxsq,av_dysq,av_Ehatsum_x,av_Ehatsum_xsq,av_Ehatsum_y,av_Ehatsum_ysq
  common  av_Ebar_x,av_Ebar_y,av_Ebar_xsq,av_Ebar_ysq,av_Ebardip_x,av_Ebardip_y,av_Ebardip_xsq,av_Ebardip_ysq
  common  av_Ebarwind_x,av_Ebarwind_y,av_Ebarwind_xsq,av_Ebarwind_ysq
  common  susceptibility,av_mx,av_my,av_spinf,av_enerf,ul,vl,av_mag
  common  theta,v
  common  top_x,top_y,SW,topDfree_x,topDfree_y
  common  n_x,n_y
  common  Ehat_x,Ehat_y,wind_x,wind_y
  common  phi,Ehatsum_x,Ehatsum_y,nx,ny,Jx,Jy,Eg_x,Eg_y,quench_fluct,delta_fluct,delta_phi,delta_delta_phi
  common  volume,length,pos,neg,side,sites
  common  therm_sweeps,sweeps,quench_sweeps,bin_size,bin_size_curr,num_bins,measurements,ibin,ibin_curr,twist,MCT,glob,globy,printU
  common  accept,accepth,acceptg,acceptCE,Tsteps,Esteps,sidesteps,sidemin,sidemax,sideincr,ratio,ratiorot,ratiowind
  common  outfile1,outfile2,outfile3,outfile4,outfile5,spinin,spinout,vortout,plot_spin_data,plot_vort_data
  common  top_snap_shot,spin_snap_shot,xtopout,xtopin,ytopout,ytopin
