#ifndef WTURBCONF_H
#define WTURBCONF_H

/*-------------------------------------------------------------------------*/
/* 									   */
/*-------------------------------------------------------------------------*/

typedef struct
{
  int ps_cg_use_old_pressure,
      ps_cg_use_old_pressure_opt;
  int ps_mg_sxyzmin,
      ps_mg_sm_iter,
      ps_mg_sm_bound_iter,
      ps_mg_sm_block_iter;
  //float ps_mg_sm_block_iter_thr;
  int ps_mg_bz_width,
      ps_mg_coarsest_sm_iter;
  int ps_mg_ar_on;
  float ps_mg_ar_theta,
	ps_mg_ar_theta_exp;
  float ps_mg_ar_testval;
  int ps_mg_ar_options,
      ps_mg_ar_use_adaptive_thr;
  int ps_tol_norm,
      ps_tol_opt;
  float ps_tol_factor_air;
  int ps_options;
}
WTURB_PSolverConf;

typedef struct
{
  int damp_on,
      slip_typ,
      slip_new_vers;
  float damp_alpha,
	damp_width;
}
WTURB_ObstConf;

typedef struct
{

  int		do_use_SI_units,
		decouple_air_pressure,
		do_physical_air_velocity,
		do_separate_air_velocity,
		have_liquid;
  double        domain_size_m,
		cell_size_m,
		max_simulation_time;
  int		domain_res,
		domain_res_ampl;
  float 	liq_density,
	        air_density,
		small_density_cutoff,
		rho_min_factor,
		air_part_bump_dist;
  char		init_part_file_name[UT_MAXPATHLEN+1];
  int           do_swap_particles;
  float	        max_particle_array_fragmentation;
  int		particle_sort;
  int	        test_quant,
    		test_quant_pvel,
	        test_quant_ppos;
  
  int		heat_on,
		heat_opt;
  float		heat_temp,
		heat_temp_amb,
		heat_diffusion;
  int		heat_dissolve_typ;
  float		heat_dissolve_rate;

  int		soot_on,
		soot_opt;
  float		soot_density,
		soot_diffusion;
  int		soot_dissolve_typ;
  float		soot_dissolve_rate;
}
WTURB_MiscConf;

typedef struct
{
  int	   on,
	   blocksize0,
	   n_levels;
  int	   refine_db,
	   refine_at_obstacles,
	   refine_heat;
  int	   thr_typ;
  double   thr_val,
	   thr_val_dens,
	   thr_span_dens,
	   thr_val_heat,
	   thr_span_heat;

  int	   thr_liq_nb;
  double   thr_liq_nb_width,
	   thr_liq_nb_width_2,
	   thr_liq_nb_width_3;
  int 	   thr_liq_drop;
  

  int	   thr_air_nb;
  double   thr_air_nb_width,
	   thr_air_nb_width_2,
	   thr_air_nb_width_3;

  int	   thr_soot_on;
  double   thr_soot_eps;

  int	   thr_vort_on;
  int	   thr_vort_opt,
	   thr_vort_min_level_open;
  float    thr_vort_levels[4],
  	   thr_vort_min_level;

  int	   rf_roi_typ;
  double   rf_roi_center[4],
	   rf_roi_radius[4];

  int	   sootsrc_typ,
	   sootsrc_opt;
  double   sootsrc_center[4],
	   sootsrc_radius[4];
  float    sootsrc_duration,
	   sootsrc_strength,
	   sootsrc_randomize;
  
  int	   heatsrc_typ;
  double   heatsrc_center[4],
	   heatsrc_radius[4];
  double   heatsrc_temp,
	   heatsrc_duration,
	   heatsrc_randomize;

  int	   inflow_typ,
	   inflow_opt;
  double   inflow_center[4],
	   inflow_radius[4],
	   inflow_velocity[4],
	   inflow_force[4],
	   inflow_duration;

  double   rbz_width_blk;
}
WTURB_MSBGConf;

typedef struct
{
  int      dva_on,
	   dva_positive_only;
  float	   dva_strength,
	   dva_min_depth;
  int      dva_ui_on,
	   dva_ui_opt;
  float	   dva_ui_divdiscount,
	   dva_ui_maxdist;
}
WTURB_DivAdjConf;

typedef struct
{
  int vc_on,
      vc_options;
  float vc_strength;  

  // D E P R E C A T E D
  int vort_vers,
      vort_on;
  double vort_eps,
    	 vort_thr,
	 vort_bias,
	 vort_gain;
  double hel_factor;
  double ada_factor;
  int use_spline_gradient; 
  int turb_on,
      turb_use_prod_only,
      use_spline_derivs,
      edge_centered_vort;
  double turb_tm_exp,
    	 turb_scale,
	 turb_offset,
	 turb_gamma;
  int hfe_on,
      hfe_bands;
  double hfe_tm_exp,
    	 hfe_scale,
	 hfe_offset,
	 hfe_gamma;

  int	 ghfe_adaptive;

  int	 ta_n_keys;
  struct 
  {
    double t,
	   strength;
  }		
  ta_keys[100];
}
WTURB_VortConfineConf;

typedef struct
{
  float maxCFL,
	maxAirCFL;
  int  rol_active;		// Remove Outliers for vmax 
  float rol_numSigma;
}
WTURB_CFLConf;

typedef struct
{
  int   typ,
	use_pres_grad_mag;
  int   iter;
  float thr,
	thr_tol1,
	thr_tol2;
  float dt;
  float velmag_thr;
  float lambda;
  float unused;
  float taubin;
  float lambda_sdf;
  float redist;
}
WTURB_LapSmConf;

typedef struct
{
  float r_large,
	r_small_mult,
	mu_sep,
	rho_kr,
	rho_kr_supp,
	rho_droplet_mult,
	rho_min;
  unsigned options;
}
WTURB_RsurfSbConf;

typedef struct
{
  int type;
  int ls_typ,
      ls_ext_nb_sdf;
  float dens_thr,
	mass2dens_clip;
  
  float lsr_twk_strength,
        lsr_twk_offset,
	lsr_twk_gamma,
	lsr_twk_gain,
	lsr_twk_bias;
  int	lsr_twk_active,
    	lsr_twk_nonlin_active;
 
  float lsr_perlin_gain;
  int   lsr_perlin_gain_active;
  float lsr_perlin_bias;
  int   lsr_perlin_bias_active;

  WTURB_LapSmConf lsr_lapsm;

  int   st_on,
	st_options;
  double st_strength,
	 st_strength_SI;
  int    st_dfss_on;
  
  int   lsr_est_typ;
  float lsr_est_thr,
	lsr_est_sig,
	lsr_est_thr_small,
	lsr_est_strength,
	lsr_est_ampl;

  int	lsr_typ,
	lsr_n_iter;
  float	lsr_pradius,
        lsr_dmax,
	lsr_lambda,
	lsr_lambda_sdf,
	lsr_taubin,
	lsr_redist,
	lsr_dt;

  int	rsurf_typ,
	rsurf_do_calc_SDF;
  float rsurf_dens_thr,
        rsurf_particle_radius,
	rsurf_drop_radius;
  int   rsurf_krn_typ,
	rsurf_do_avgpos,
	rsurf_do_upsample_orig,
	rsurf_version,
	rsurf_test;
  float rsurf_randomize_radius;
  float rsurf_nb_dist,
	rsurf_dens_soft_clip;

  int	rsurf_drop_typ;
  float rsurf_drop_sizelevel_mult;
  int   rsurf_drop_rad_rand;
  float rsurf_drop_rad_rand_sig,
	rsurf_drop_rad_rand_min,
	rsurf_drop_rad_rand_max;

  int rsurf_drop_tiny;
  float rsurf_drop_tiny_radius,
	rsurf_drop_tiny_max_radius,
	rsurf_drop_tiny_sharpness,
	rsurf_drop_tiny_gamma,
	rsurf_drop_tiny_atten;
  int   rsurf_drop_tiny_krn,
	rsurf_drop_tiny_radius_opt;

  float rsurf_twk_strength,
        rsurf_twk_offset,
	rsurf_twk_gamma,
	rsurf_twk_gain,
	rsurf_twk_bias;
  int	rsurf_twk_active,
    	rsurf_twk_nonlin_active;
  int	ww_on;
  int   ww_crv_typ;
  float ww_crv_radius,
	ww_crv_offs,
	ww_crv_post_offs,
	ww_crv_dens_thr,
	ww_crv_dens_thr2,
	ww_crv_eps,
	ww_crv_gain,
	ww_crv_gamma;
  int   ww_crv_opt;
  int   ww_crv_radius_krn_supp;
  int   ww_soot_typ;
  float ww_soot_strength;
  int   ww_mist_as_tiny_droplets;
  float ww_mist_as_tiny_droplets_strength;

  int	ww_diff_typ,
	ww_diff_num_iter;
  float ww_diff_strength,
	ww_diff_anis;

  int air_typ,
      air_opt;
  int air_extrapol_liq_vel;
  float air_extrapol_liq_vel_theta0;
  float air_flip_alpha[8];

  float air_direction[4],
	air_speed_ref,
	air_speed_profile_alpha,
	air_speed_profile_height_zero;

  int	air_compr_on;
  float air_compr_speed_of_sound;

  int air_amp_typ;
  float air_amp_strength,
	air_amp_scale;

  int air_proc_typ,
      air_proc_vers;
  float air_proc_alpha,
	air_proc_beta,
	air_proc_beta_time,
	air_proc_oct,
	air_proc_oct_leadin,
	air_proc_oct_offs,
	air_proc_seed;
  float air_proc_fscale,
	air_proc_fscale_min_cells,
	air_proc_fscale_time;
  float air_proc_ampl[4];

  int   air_proc_ampl_hf_on,
	air_proc_ampl_hf_adapt,
	air_proc_ampl_hf_adapt_with_main;
  float air_proc_ampl_hf_gain,
	air_proc_ampl_hf_maxscale,
	air_proc_ampl_hf_adapt_velref_SI[4],
	air_proc_ampl_hf_adapt_maxdist;
  int   air_proc_ampl_hf_adapt_opt;

  float air_proc_obst_nbdist,
	air_proc_obst_proj_lf_scale;
  unsigned air_proc_options;

  int	air_mist_typ;
  float air_mist_strength,
        air_mist_gravity,
	air_mist_dissolve;

  WTURB_RsurfSbConf rsurf_sb;

  WTURB_LapSmConf rsurf_lapsm;
  int	rsurf_lapsm_vd_on,
	rsurf_lapsm_vd_use_dfss;
  float rsurf_lapsm_vd_fit[4];

  int rsurf_ldn_typ;	// Particle local density normalization
  float rsurf_ldn_sig;

  int rsurf_est_on,
      rsurf_est_do_apply_after_surfsmooth;
  float rsurf_est_thr,
	rsurf_est_sig,
	rsurf_est_thr_small,
	rsurf_est_strength,
	rsurf_est_ampl;
  
  int gp_active,
      gp_type;

  int kill_isolated_nonfluid_part;
  float kill_isolated_nonfluid_part_dist;

  int drag_on;
  double drag_strength, drag_strength_SI,
	drag_strength_obst, drag_strength_obst_SI,
	drag_strength_obst_dist;

  int drag_vcfl_on;
  float drag_vcfl_val;

  int drag_dfss_on,
      drag_dfss_cubic;
  float drag_dfss_strength,
	drag_dfss_blob_radius_SI;

  int drag_hfe_on;
  float drag_hfe_enrg_min,
	drag_hfe_enrg_max;
  float drag_hfe_ampl;
  int   drag_hfe_opt;

  int	drop_on;
  unsigned drop_opt;
  float drop_thr_curv,
	drop_thr_hfe;

  float drop_gen_timescale;
  int drop_gen_dfss_on;
  float drop_gen_dfss_strength,
	drop_gen_dfss_bup_thr,
	drop_gen_dfss_bup_time_factor;

  int drop_gen_hfe_on;
  float drop_gen_hfe_min,
	drop_gen_hfe_max;

  int   drop_rr_on;
  float drop_rr_sig,
	drop_rr_min,
	drop_rr_max,
	drop_rr_mean,
	drop_rr_base_radius_SI;
  int   drop_rr_typ;
  float	drop_sizelevel_mult;

  int	drop_evap;
  float drop_evap_strength,
	drop_evap_delay,
	drop_evap_maxrate;

  int	drop_hfe;
  float drop_hfe_decay;

  int   drop_mult,
	drop_mult_old_vers;
  float drop_mult_disp;
  int   drop_mult_disp_opt;
  float	drop_mult_disp_vel,
	drop_mult_em_radius,
	drop_mult_em_length,
	drop_mult_em_vel_ref_SI,
	drop_mult_em_vel_ref_range[4];
  float drop_mult_sigmaN;
  int   drop_mult_em_vel_ref_use_range,
	drop_mult_em_vel_ref_use_energy,
	drop_mult_use_proc_disp_vel,
	drop_mult_em_time_jitter;
  float drop_mult_em_offset;

  int   drop_mult_proc_vel_typ;
  float drop_mult_proc_vel_fscale,
	drop_mult_proc_vel_fscale_min,
	drop_mult_proc_vel_leadoct,
	drop_mult_proc_vel_fscale_time,
	drop_mult_proc_vel_alpha,
	drop_mult_proc_vel_ampl;

  int   drop_mult_proc_pos_typ;
  float drop_mult_proc_pos_fscale,
	drop_mult_proc_pos_leadoct,
	drop_mult_proc_pos_fscale_min,
	drop_mult_proc_pos_fscale_time,
	drop_mult_proc_pos_alpha,
	drop_mult_proc_pos_ampl;

  int	drop_dist_typ;
  float drop_dist,
	drop_dens,
	drop_dens_secondary_only,
	drop_blend_vel;
  int   drop_kill;
  float drop_kill_dist,
	drop_kill_dens,
	drop_kill_lifetime;
  int   drop_kill_decay_on,
	drop_kill_decay_opt;
  int   drop_rejoin;
  float drop_rejoin_dist,
        drop_rejoin_blend_vel,
	drop_rejoin_blend_width,
	drop_rejoin_blend_rmin,
	drop_rejoin_blend_rmax;
  int   drop_drag_on, 
	drop_drag_with_rad;
  int   drop_drag_deform_on;
  float drop_drag_deform_strength,
  	drop_drag_deform_rmin,
  	drop_drag_deform_rmax;
  int   drop_airbck;
  float drop_airbck_strength;

  double drop_drag_strength, drop_drag_strength_SI,
	drop_drag_strength_obst, drop_drag_strength_obst_SI,
	drop_drag_strength_obst_dist;
  int   drop_rndf_on,
	drop_rndf_do_resindep,
	drop_rndf_use_sqrt_dt,
	drop_rndf_opt;
  double drop_rndf_strength, drop_rndf_strength_SI;

  int   dfss_typ,
	dfss_blur_level;
  float dfss_epsilon,
	dfss_theta,
	dfss_dens_thr;

  int	rndf_on;
  float rndf_strength,
	rndf_strength_dfss;
}
WTURB_LiquidConf;

typedef struct
{
  int pic;
  int sdf;
  int rsurf;
  int drag;
  int grid_topology,
      grid_topology_obst,
      block_flags_rgb[3];
  int obst_sdf,
      obst_obj_id,
      obst_face_area,
      obst_sdf_grad_mag,
      obst_sdf_grad;
  int vis_mdens,
      vis_pf;
  int cur_on,
      cur_level,
      cur_ix,
      cur_iy,
      cur_iz;
  double roi_slice_pos_[3],
	 roi_box_pos[3],
	 roi_box_size,
	 roi_box_zoom;
  int tm_mod;
  int psv_opt,
      psv_mg,
      psv_cg;
  int psv_ar_on,
      psv_ar_opt;
  int pc_on,
      pc_stage,
      pc_opt;
  int ww;
  int vis_drop;
  int vis_turb,
      vis_turb_amp,
      vis_turb_opt;
  float vis_turb_maxval,
	vis_turb_thr,
	vis_turb_thr2,
	vis_turb_thr_liq,
	vis_turb_thr2_liq;
  int   vis_vmag_opt;
  float vis_vmag_maxval,
	vis_vmag_thr,
	vis_vmag_thr2,
	vis_vmag_thr_liq,
	vis_vmag_thr2_liq;
  int vis_dfss;

  int vis_velstats_on,
      vis_velstats_nbin;

  float vis_vel_arrow_scale;
  int vis_vel_arrow_freq;   // arrows per cell

  int vis_part, // _visParticles,
      vis_part_mod, //	 _visParticlesMod,
      vis_part_opt, //_visParticlesOpt,
      vis_part_opt3,
      vis_part_when;  
  int vis_part_maxres_cell,
      vis_part_maxres;

  int vis_pressure,
      vis_pressure_ps;
}
WTURB_VisualizationConf;

typedef struct
{
  int	   on,
	   ppc_avg,
	   ppc_min,
	   ppc_max,
	   ppc_max_zerodens;
  double   alpha,
	   alpha_dens,
	   alpha_soot;
  float    alpha_levels[8];
  int	   alpha_coarse_levels_on;

  int		p2g_kerntyp,
		p2g_kerntyp_dns;
  float	        p2g_kern_radius;

  int		p2g_sh_typ,
		p2g_sh_accumul;

  double	p2g_sh_radius,
		p2g_sh_strength;
  int		p2g_sh_d_typ,
		p2g_sh_d_accumul;
  double	p2g_sh_d_radius,
		p2g_sh_d_strength;

  int		ivb_on;
  double	ivb_strength;
  int		ivb_krn_exp;

  int		sgc_on,
		sgc_plevel0_only,
		sgc_timestep_indep;
  double        sgc_stepsize;
  float         sgc_rnd_strength,
		sgc_rnd_freq,
		sgc_max_step;
  int           sgc_min_depth;

  int		rs_ipol_meth,
		rs_ipol_vel_meth;
  double        rs_prob;

  int		prs_seed_ipol;

  int 		adv_vel_ipol,
		adv_vel_ipol_divfree,
      		adv_trc_meth,
		adv_trc_max_substeps;
  double      	adv_trc_dxmax;
  int		adv_trc_dxmax_leveldep,
    		adv_trc_drop_meth,
		adv_trc_drop_expfcalcmod;
  float         adv_trc_drop_dxmax;

  int		adv_tjit_active;
  float         adv_tjit_frac;
  int	        g2p_ipol;

  int		ipvec_typ,
		ipvec_fadewgt;

  int 		ipdivf_on,
    		ipdivf_on_G2P_,
       		ipdivf_on_advect_,
       		ipdivf_on_seed_,
		ipdivf_on_dropgen_;

  int		sm_factor,
		sm_free_velo_grid;
  float         init_velocity[3],
  		init_velocity_SI[3],
		init_velocity_bbox[6];
  int		init_velocity_opt_;
  int		init_velocity_typ,
		use_synthvel;
  double        init_velocity_params[20];
  float         init_velocity_air[3],
  		init_velocity_air_SI[3],
		init_velocity_air_bbox[6];
  int		init_velocity_air_opt;
}
WTURB_FLIPConf;

#endif /* WTURBCONF_H */
