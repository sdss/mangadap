pro manga_dap,input_number,configure_file

;check_version=check_version,$
;       dont_remove_null_templates=dont_remove_null_templates;,datacubes_or_rss=datacubes_or_rss

t0=systime(/seconds)/60./60.
c=299792.458d
;*****SETTING THE CONFIGURATION PARAMETERS DEFINED IN THE CONFIGURATION FILE
readcol,configure_file,command_line,comment='#',delimiter='%',/silent,format='A'
for i=0, n_elements(command_line)-1 do d=execute(command_line[i])
if n_elements(save_intermediate_steps) eq 0 then save_intermediate_steps=0
if n_elements(remove_null_templates) eq 0 then remove_null_templates = 1



if n_elements(instrumental_fwhm_file) ne 0 then readcol, instrumental_fwhm_file, ww,r,fwhm_ang,fwhm_kms,/silent;'../../instrumental_fwhm/instrumental_fwhm.dat'

if keyword_set(check_version) then begin
   readcol,total_filelist,root_name_vector,velocity_initial_guess_vector,$
          velocity_dispersion_initial_guess_vector,ellipticity_vector,position_angle_vector,fibre_number_vector,Reff_vector,mode_vector,/silent, format='A,F,F,F,F,F,F,A',comment='#'
   root_name = root_name_vector[input_number]
   mode = mode_vector[input_number]
   output_dir='results/'+mode+'/'+root_name+'/'+root_name+'_'
   output_idlsession=output_dir+'mdap_session.idl'
   check_previous_session = file_test(output_idlsession)
   if check_previous_session eq 1 then restore,output_idlsession
endif

manga_dap_version = '0.7'    ; 28 Mar 2014 by L. Coccato
;stop

readcol,total_filelist,root_name_vector,velocity_initial_guess_vector,$
        velocity_dispersion_initial_guess_vector,ellipticity_vector,position_angle_vector,fibre_number_vector,Reff_vector,mode_vector,/silent, format='A,F,F,F,F,F,F,A',comment='#'


root_name = root_name_vector[input_number]
mode = mode_vector[input_number]
datacube_dir = datacube_root_dir+mode[0]+'/'
velocity_initial_guess = velocity_initial_guess_vector[input_number]
velocity_dispersion_initial_guess = velocity_dispersion_initial_guess_vector[input_number]
ell=ellipticity_vector[input_number]
pa=position_angle_vector[input_number]
number_of_fibres=fibre_number_vector[input_number]
Reff=Reff_vector[input_number]
output_dir=output_root_dir+'results_'+mode+'/'+root_name+'/'+root_name+'_'
datacube_name=root_name+'.fits'  
if mode eq 'datacubes' then begin
   datacube_name=root_name+'.fits'  
   res = file_test(output_root_dir+'results_datacubes/'+root_name+'/',/directory)
   if res eq 0 then spawn,'mkdir '+output_root_dir+'results_datacubes/'+root_name
endif
if mode eq 'rss' then begin
   datacube_name=root_name+'_rss.fits'
   res = file_test(output_root_dir+'results_rss/'+root_name+'/',/directory)
   if res eq 0 then spawn,'mkdir '+output_root_dir+'results_rss/'+root_name
endif
output_filefits=output_dir+'high_level.fits'
output_idlsession=output_dir+'mdap_session.idl'


print, ''
print, '# WORKING ON '+root_name+' ('+mode+')'
print, ''
openw,1,output_dir+'mdap.log'


;BLOCK 0 
;*** GET MODULES VERSION AND CHECK PREVIOUS ANALYSIS ******************************
check_previous_analysis = file_test(output_filefits)
check_previous_session = file_test(output_idlsession)
if check_previous_analysis eq 1 and check_previous_session eq 1 and keyword_set(check_version) then begin
   tmp_header = headfits(output_filefits)
   manga_dap_version_previous = strcompress(sxpar(tmp_header,'DAPVER'),/remove_all)

   mdap_read_datacube_version=strcompress(sxpar(tmp_header,'BLOCK1'),/remove_all)
   mdap_read_datacube_version_previous = mdap_read_datacube_version
   mdap_spatial_binning_version=strcompress(sxpar(tmp_header,'BLOCK2'),/remove_all)
   mdap_spatial_binning_version_previous = mdap_spatial_binning_version
   mdap_log_rebin_version=strcompress(sxpar(tmp_header,'BLOCK3'),/remove_all)
   mdap_log_rebin_version_previous=mdap_log_rebin_version
   mdap_spectral_fitting_version=strcompress(sxpar(tmp_header,'BLOCK4'),/remove_all)
   mdap_spectral_fitting_version_previous=mdap_spectral_fitting_version
   mdap_measure_indices_version=strcompress(sxpar(tmp_header,'BLOCK5'),/remove_all)
   mdap_measure_indices_version_previous=mdap_measure_indices_version
   mdap_spatial_radial_binning_version=strcompress(sxpar(tmp_header,'BLOCK6'),/remove_all)
   mdap_spatial_radial_binning_version_previous=mdap_spatial_radial_binning_version

endif else begin
   manga_dap_version_previous = '0'

   mdap_read_datacube_version='0'
   mdap_read_datacube_version_previous='0'
   mdap_spatial_binning_version='0'
   mdap_spatial_binning_version_previous='0'
   mdap_log_rebin_version='0'
   mdap_log_rebin_version_previous='0'
   mdap_spectral_fitting_version='0'
   mdap_spectral_fitting_version_previous='0'
   mdap_measure_indices_version='0'
   mdap_measure_indices_version_previous='0'
   mdap_spatial_radial_binning_version='0'
   mdap_spatial_radial_binning_version_previous='0'
endelse


mdap_read_datacube,version=mdap_read_datacube_version
mdap_spatial_binning,version=mdap_spatial_binning_version
mdap_log_rebin,version=mdap_log_rebin_version
mdap_spectral_fitting,version=mdap_spectral_fitting_version
mdap_measure_indices,version=mdap_measure_indices_version
mdap_spatial_radial_binning,version=mdap_spatial_radial_binning_version
;junk = temporary(chek_version)
;**********************************************************************************
 

;; KINEMETRY CHECK, remove it after the tests
; restore, output_idlsession
; goto, kinemetry_step
;;


execute_all_modules=0
if manga_dap_version gt manga_dap_version_previous then execute_all_modules=1

;goto, inizia_da_qui
; BLOCK 1
;*** READING THE INPUT DATACUBE* **************************************************
if mdap_read_datacube_version gt mdap_read_datacube_version_previous or execute_all_modules eq 1 then begin
   mdap_read_datacube,datacube_dir+datacube_name,data,error,wavelength,x2d,y2d,signal,noise,cdelt1,cdelt2,header_2d,lrange=w_range_for_sn_computation,/keep_original_step,$
     x2d_reconstructed=x2d_reconstructed,y2d_reconstructed=y2d_reconstructed,signal2d_reconstructed=signal2d_reconstructed,noise2d_reconstructed=noise2d_reconstructed,$
    number_of_fibres=number_of_fibres
      if n_elements(signal2d_reconstructed) eq 0 then signal2d_reconstructed = signal   ;same as signal if datacube is read, not RSS.
      if n_elements(x2d_reconstructed) eq 0 then x2d_reconstructed = x2d                ;same as x2d if datacube is read, not RSS.
      if n_elements(y2d_reconstructed) eq 0 then y2d_reconstructed = y2d                ;same as y2d if datacube is read, not RSS.
      if n_elements(noise2d_reconstructed) eq 0 then noise2d_reconstructed = noise      ;same as noise if datacube is read, not RSS.
   sxaddpar,header_2d,'BLOCK1',mdap_read_datacube_version,'mdap_read_datacube version'
   execute_all_modules=1

  ;-- I compute the signal and the noise over a different wavelength range (redshifted)
  ; The mdap_read_datacube is run again, but only the signal and noise output are stored
  if n_elements(w_range_for_sn_computation_for_gas) ne 0 then begin
     w_range_for_sn_computation_for_gas_redshifted = w_range_for_sn_computation_for_gas*(1.+velocity_initial_guess[0]/299792.5)
     mdap_read_datacube,datacube_dir+datacube_name,junk,junk,junk,junk,junk,signal_for_gas,noise_for_gas,junk,junk,header_junk,$
               lrange=w_range_for_sn_computation_for_gas_redshifted,/keep_original_step,/use_total
      junk = 0.
  endif else begin
     signal_for_gas = signal
     noise_for_gas = noise
     w_range_for_sn_computation_for_gas = w_range_for_sn_computation
  endelse
  ;--

endif

sxaddpar,header_2d,'DAPVER',manga_dap_version,'manga_dap version',BEFORE  ='BLOCK1'
printf,1,'[INFO] mdap_read_datacube ver '+max([mdap_read_datacube_version,mdap_read_datacube_version_previous])
sz=size(data)
if sz[0] eq 3 then printf,1,'[INFO] datacube '+root_name+' read; size: '+mdap_stc(sz[1],/integer)+'x'+mdap_stc(sz[2],/integer)
if sz[0] eq 2 then printf,1,'[INFO] RSS file '+root_name+' read; Nfibres: '+mdap_stc(sz[1],/integer)

;**** definition of some variables
MW_extinction = 0.                    ; will be replaced by appropriate extension in the input file, or total_filelist
;fwhm_instr=wavelength*0.+2.73         ; will be replaced by appropriate extension in the input file
if n_elements(instrumental_fwhm_file) ne 0 then begin
   fwhm_instr=interpol(fwhm_ang,ww,wavelength) 
   fwhm_instr_kmsec_matrix = fltarr(2,n_elements(ww))
   fwhm_instr_kmsec_matrix[0,*]=ww
   fwhm_instr_kmsec_matrix[1,*]=fwhm_kms
 endif else begin
    fwhm_instr=wavelength*0.+2.54
    ;fwhm_instr_kmsec_matrix = fltarr(2,n_elements(wavelength))
    ;fwhm_instr_kmsec_matrix[0,*]=wavelength
    ;fwhm_instr_kmsec_matrix[1,*]=mean(fwhm_instr/wavelength*300000.)
endelse

mask_range=[5570.,5590.,5800.,5850.]  ; will be replaced by mask extension in the input file

;*** END OF BLOCK 1 *************************************************************
;if save_intermediate_steps eq 1 then save,filename=root_name+mode+'_block1.idl',/variables 
;stop




;BLOCK 2
;*** SPATIAL BINNING *************************************************************
if mode eq 'rss' then begin
   sn1=sn1_rss
   sn2=sn2_rss
   sn3=sn3_rss
   sn_thr_tpl=sn_thr_tpl_rss
   sn_thr_str=sn_thr_str_rss
   sn_thr_ems=sn_thr_ems_rss
   if n_elements(sn_calibration_rss) ne 0 then sn_calibration = sn_calibration_rss
endif
if mode eq 'datacubes' then begin
   sn1=sn1_datacubes
   sn2=sn2_datacubes
   sn3=sn3_datacubes
   sn_thr_tpl=sn_thr_tpl_datacubes;4.
   sn_thr_str=sn_thr_str_datacubes
   sn_thr_ems=sn_thr_ems_datacubes
   if n_elements(sn_calibration_datacubes) ne 0 then sn_calibration = sn_calibration_datacubes
endif

if mdap_spatial_binning_version gt mdap_spatial_binning_version_previous  or execute_all_modules eq 1 then begin
   mdap_spatial_binning,data,error,signal,noise,sn1,x2d,y2d,cdelt1,cdelt2,spatial_binning_tpl,spectra_tpl,errors_tpl,$
                        xbin_tpl,ybin_tpl,area_bins_tpl,bin_sn_tpl,sn_thr=sn_thr_tpl,$
                         x2d_reconstructed=x2d_reconstructed,y2d_reconstructed=y2d_reconstructed,$
                         nelements_within_bin=nelements_within_bin_tpl,sn_calibration=sn_calibration,$
                         user_bin_map=user_bin_map_spatial_binning_1,weight_for_sn=weight_for_sn;,/plot
   sxaddpar,header_2d,'BLOCK2',mdap_spatial_binning_version,'mdap_spatial_binning version'
   execute_all_modules=1
endif
;stop
printf,1,'[INFO] mdap_spatial_binning ver '+max([mdap_spatial_binning_version,mdap_spatial_binning_version_previous])
printf,1,'[INFO] datacube '+ root_name+' spatial binning 1. SN= '+mdap_stc(sn1,/integer)+' Nbins: '+mdap_stc(n_elements(xbin_tpl),/integer)
if mdap_spatial_binning_version gt mdap_spatial_binning_version_previous  or execute_all_modules eq 1 then begin
   mdap_spatial_binning,data,error,signal,noise,sn2,x2d,y2d,cdelt1,cdelt2,spatial_binning_str,spectra_str,errors_str,$
                        xbin_str,ybin_str,area_bins_str,bin_sn_str,sn_thr=sn_thr_str,$
                        x2d_reconstructed=x2d_reconstructed,y2d_reconstructed=y2d_reconstructed,$
                        nelements_within_bin=nelements_within_bin_str,sn_calibration=sn_calibration,$
                        user_bin_map=user_bin_map_spatial_binning_2,weight_for_sn=weight_for_sn;,/plot
   execute_all_modules=1 
endif
printf,1,'[INFO] datacube '+ root_name+' spatial binning 2. SN= '+mdap_stc(sn2,/integer)+' Nbins: '+mdap_stc(n_elements(xbin_str),/integer)

;stop

if mdap_spatial_binning_version gt mdap_spatial_binning_version_previous  or execute_all_modules eq 1 then begin
   mdap_spatial_binning,data,error,signal_for_gas,noise_for_gas,sn3,x2d,y2d,cdelt1,cdelt2,spatial_binning_ems,spectra_ems,errors_ems,$
                        xbin_ems,ybin_ems,area_bins_ems,bin_sn_ems,sn_thr=sn_thr_ems,$
                        x2d_reconstructed=x2d_reconstructed,y2d_reconstructed=y2d_reconstructed,$
                        nelements_within_bin=nelements_within_bin_ems,sn_calibration=sn_calibration,$
                        user_bin_map=user_bin_map_spatial_binning_3,weight_for_sn=weight_for_sn;,/plot
   execute_all_modules=1
endif
printf,1,'[INFO] datacube '+ root_name+' spatial binning 3. SN= '+mdap_stc(sn3,/integer)+' Nbins: '+mdap_stc(n_elements(xbin_ems),/integer)
;stop

;*** END OF BLOCK 2 *************************************************************




; BLOCK 3
;*** LOG-REBIN, STELLAR TEMPLATE CONVOLUTION, SELECTION OF WAVELENGTH RANGE ******
;stop

fwhm_stars = wavelength*0.+2.54   ; To be used for MARCS templates.
fwhm_diff=sqrt(fwhm_instr^2 - fwhm_stars^2)
indici = where(finite(fwhm_diff) ne 1)
if indici[0] ne -1 then fwhm_diff[indici] = 0
if mdap_log_rebin_version gt mdap_log_rebin_version_previous   or execute_all_modules eq 1 then begin
   mdap_log_rebin,spectra_tpl,errors_tpl,wavelength,stellar_library_spatial_binning_1,fwhm_diff,$
        log_spc_tpl,log_err_tpl,log_wav_tpl,library_log_tpl,log_wav_library_tpl,$
        input_velscale=velscale,wave_range=trim_wav_range_spatial_binning_1,/gal_wavelength_log_step,/quiet
   sxaddpar,header_2d,'BLOCK3',mdap_log_rebin_version,'mdap_log_rebin version'
   execute_all_modules=1
endif


printf,1,'[INFO] mdap_log_rebin ver '+max([mdap_log_rebin_version,mdap_log_rebin_version_previous])
printf,1,'[INFO] datacube '+root_name+'log_rebin 1'
if mdap_log_rebin_version gt mdap_log_rebin_version_previous   or execute_all_modules eq 1 then begin
   mdap_log_rebin,spectra_str,errors_str,wavelength,stellar_library_spatial_binning_2,fwhm_diff,$
        log_spc_str,log_err_str,log_wav_str,library_log_str,log_wav_library_str,$
        input_velscale=velscale,wave_range=trim_wav_range_spatial_binning_2,/gal_wavelength_log_step,/quiet
   execute_all_modules=1
endif
printf,1,'[INFO] datacube '+root_name+'log_rebin 2'
if mdap_log_rebin_version gt mdap_log_rebin_version_previous  or execute_all_modules eq 1 then begin
   mdap_log_rebin,spectra_ems,errors_ems,wavelength,stellar_library_spatial_binning_3,fwhm_diff,$
        log_spc_ems,log_err_ems,log_wav_ems,library_log_ems,log_wav_library_ems,$
        input_velscale=velscale,wave_range=trim_wav_range_spatial_binning_3,/gal_wavelength_log_step,/quiet
   execute_all_modules=1
endif
printf,1,'[INFO] datacube '+root_name+'log_rebin 3'
;
;save,filename='block3.idl',/variables
;

library_log_tpl = library_log_tpl/10.^25.
library_log_str = library_log_str/10.^25.
library_log_ems = library_log_ems/10.^25.
if save_intermediate_steps eq 1 then save,filename=root_name+mode+'_block3.idl',/variables 
;*** END OF BLOCK 3 *************************************************************


;
;    ##############################################
;    #####      HIGH LEVEL DATA PRODUCTS      #####
;    ##############################################
;
;stop
;help, spatial_binning_ems
; BLOCK 4
;*** SPECTRAL-FITTING PROCEDURES *************************************************

;--full spectral fitting 

star_kin_starting_guesses = fltarr(n_elements(xbin_tpl),4)
star_kin_starting_guesses[*,0] = velocity_initial_guess[0];str guess for velocity
star_kin_starting_guesses[*,1]=velocity_dispersion_initial_guess[0];200.  ;str guess for sigma
gas_kin_starting_guesses = star_kin_starting_guesses[*,0:1];fltarr(n_elements(xbin_str),2)
gas_kin_starting_guesses[*,1]=50.  ;str guess for sigma


if mdap_spectral_fitting_version gt mdap_spectral_fitting_version_previous or execute_all_modules eq 1 then begin
   mdap_spectral_fitting,log_spc_tpl,log_err_tpl,log_wav_tpl,library_log_tpl,log_wav_library_tpl,velscale,$
        stellar_kinematics_tpl,stellar_kinematics_tpl_err,stellar_weights_tpl,emission_line_kinematics_tpl,emission_line_kinematics_tpl_err,$
         emission_line_kinematics_tpl_individual,emission_line_kinematics_tpl_individual_err,$
        emission_line_intens_tpl,emission_line_intens_tpl_err,emission_line_fluxes_tpl,emission_line_fluxes_tpl_err,emission_line_equivW_tpl,emission_line_equivW_tpl_err,wavelength_input=exp(log_wav_library_tpl),$
        wavelength_output_tpl,best_fit_model_tpl,galaxy_minus_ems_fit_model_tpl,best_template_tpl,best_template_LOSVD_conv_tpl,reddening_tpl,reddening_tpl_err,residuals_tpl,$
        star_kin_starting_guesses=star_kin_starting_guesses,gas_kin_starting_guesses=gas_kin_starting_guesses,$
        emission_line_file=emission_line_file_spatial_binnin_1,fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix,$
        extra_inputs=spectra_fittin_parameters_patial_binning_1,mask_range=mask_range,external_library=external_library,/quiet ;,$
   ;extra_inputs=['MOMENTS=4','DEGREE=-1','MDEGREE=4']
   junk = size(emission_line_fluxes_tpl);need to save all flux maps (warning: sgandalf computes intensities, not fluxes)
   for i = 0, junk[1]-1 do emission_line_intens_tpl[i,*]=emission_line_intens_tpl[i,*]/area_bins_tpl[i]
   for i = 0, junk[1]-1 do emission_line_intens_tpl_err[i,*]=emission_line_intens_tpl_err[i,*]/sqrt(area_bins_tpl[i])
   for i = 0, junk[1]-1 do emission_line_fluxes_tpl[i,*]=emission_line_fluxes_tpl[i,*]/area_bins_tpl[i]
   for i = 0, junk[1]-1 do emission_line_fluxes_tpl_err[i,*]=emission_line_fluxes_tpl_err[i,*]/sqrt(area_bins_tpl[i])

   ;--retain only templates with positive weights for the next fits
   if remove_null_templates[0] eq 1 then begin
       wp = total(stellar_weights_tpl,1)
       library_log_tpl = temporary(library_log_tpl[*,where(wp gt 0)])
       library_log_str = temporary(library_log_str[*,where(wp gt 0)])
       library_log_ems = temporary(library_log_ems[*,where(wp gt 0)])
       stellar_weights_tpl = temporary(stellar_weights_tpl[*,where(wp gt 0)])
    endif
   execute_all_modules=1
   ;--

   ;--calculate the real S/N of binned spectra
   bin_sn_tpl_real = bin_sn_tpl
   indxx = where(wavelength_output_tpl ge w_range_for_sn_computation[0] and wavelength_output_tpl le w_range_for_sn_computation[1])
   for i = 0, n_elements(bin_sn_tpl_real) -1 do begin
       mdap_calculate_spectrum_sn,best_fit_model_tpl[i,indxx],residuals_tpl[i,indxx],wavelength_output_tpl[indxx],sn_per_angstorm,/rms
       bin_sn_tpl_real[I] = sn_per_angstorm[0]
   endfor
   ;--
   sxaddpar,header_2d,'BLOCK4',mdap_spectral_fitting_version,'mdap_spectral_fitting_version'
endif
printf,1,'[INFO] mdap_spectral_fitting ver '+max([mdap_spectral_fitting_version,mdap_spectral_fitting_version_previous])
printf,1,'[INFO] datacube '+root_name+' spectral fitting 1'
;--
;goto, block5
if save_intermediate_steps eq 1 then save,filename=root_name+mode+'_block4a.idl',/variables 
;stop

;--stellar kinematics
;   interpolate stellar kinematics results over ems grid to get starting guesses from the previous results
mdap_create_starting_guesses,stellar_kinematics_tpl,xbin_tpl,ybin_tpl,x2d,y2d,xbin_str,ybin_str,star_kin_starting_guesses,$
               velocity_initial_guess[0],velocity_dispersion_initial_guess[0],0.,0.,/h3h4
mdap_create_starting_guesses,emission_line_kinematics_tpl,xbin_tpl,ybin_tpl,x2d,y2d,xbin_str,ybin_str,gas_kin_starting_guesses,$
               velocity_initial_guess[0],50.,0.,0.

if mdap_spectral_fitting_version gt mdap_spectral_fitting_version_previous  or execute_all_modules eq 1 then begin




   ; extra FIT with only two kinematic moments, if required
   if n_elements(spectra_fittin_parameters_patial_binning_2_two_moments) ne 0 then begin
   mdap_spectral_fitting,log_spc_str,log_err_str,log_wav_str,library_log_str,log_wav_library_str,velscale,$
        stellar_kinematics_str_two_moments,stellar_kinematics_str_two_moments_err,stellar_weights_str,emission_line_kinematics_str,emission_line_kinematics_str_err,$
         emission_line_kinematics_str_individual,emission_line_kinematics_str_individual_err,$
        emission_line_intens_str,emission_line_intens_str_err,emission_line_fluxes_str, emission_line_fluxes_str_err,emission_line_equivW_str,emission_line_equivW_str_err,wavelength_input=exp(log_wav_library_str),$
        wavelength_output_str,best_fit_model_str,galaxy_minus_ems_fit_model_str,best_template_str,best_template_LOSVD_conv_str,reddening_str,reddening_str_err,residuals_str,$
        star_kin_starting_guesses=star_kin_starting_guesses,gas_kin_starting_guesses=gas_kin_starting_guesses,fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix,$
        emission_line_file=emission_line_file_spatial_binnin_2,extra_inputs=spectra_fittin_parameters_patial_binning_2_two_moments,$
        mask_range=mask_range,external_library=external_library,/quiet
   endif
   ;

   mdap_spectral_fitting,log_spc_str,log_err_str,log_wav_str,library_log_str,log_wav_library_str,velscale,$
        stellar_kinematics_str,stellar_kinematics_str_err,stellar_weights_str,emission_line_kinematics_str,emission_line_kinematics_str_err,$
         emission_line_kinematics_str_individual,emission_line_kinematics_str_individual_err,$
        emission_line_intens_str,emission_line_intens_ems_str,emission_line_fluxes_str, emission_line_fluxes_str_err,emission_line_equivW_str,emission_line_equivW_str_err,wavelength_input=exp(log_wav_library_str),$
        wavelength_output_str,best_fit_model_str,galaxy_minus_ems_fit_model_str,best_template_str,best_template_LOSVD_conv_str,reddening_str,reddening_str_err,residuals_str,$
        star_kin_starting_guesses=star_kin_starting_guesses,gas_kin_starting_guesses=gas_kin_starting_guesses,fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix,$
        emission_line_file=emission_line_file_spatial_binnin_2,extra_inputs=spectra_fittin_parameters_patial_binning_2,$
        mask_range=mask_range,external_library=external_library,/quiet
;extra_inputs=['MOMENTS=4','DEGREE=-1','MDEGREE=4']
   junk = size(emission_line_fluxes_str);need to save all flux maps
   for i = 0, junk[1]-1 do emission_line_intens_str[i,*]=emission_line_intens_str[i,*]/area_bins_str[i]
   for i = 0, junk[1]-1 do emission_line_intens_str_err[i,*]=emission_line_intens_str_err[i,*]/sqrt(area_bins_str[i])
   for i = 0, junk[1]-1 do emission_line_fluxes_str[i,*]=emission_line_fluxes_str[i,*]/area_bins_str[i]
   for i = 0, junk[1]-1 do emission_line_fluxes_str_err[i,*]=emission_line_fluxes_str_err[i,*]/sqrt(area_bins_str[i])

   ;calculate the real S/N of binned spectra
   bin_sn_str_real = bin_sn_str
   indxx = where(wavelength_output_str ge w_range_for_sn_computation[0] and wavelength_output_str le w_range_for_sn_computation[1] )
   for i = 0, n_elements(bin_sn_str_real) -1 do begin
       mdap_calculate_spectrum_sn,best_fit_model_str[i,indxx],residuals_str[i,indxx],wavelength_output_str[indxx],sn_per_angstorm,/rms
       bin_sn_str_real[I] = sn_per_angstorm[0]
   endfor

   execute_all_modules = 1
endif
printf,1,'[INFO] datacube '+root_name+' spectral fitting 2'
;--
if save_intermediate_steps eq 1 then save,filename=root_name+mode+'_block4b.idl',/variables 
;stop
;--emission lines kinematics
;interpolate stellar kinematics results over ems grid to get starting guesses from the previous results.
mdap_create_starting_guesses,stellar_kinematics_str,xbin_str,ybin_str,x2d,y2d,xbin_ems,ybin_ems,star_kin_starting_guesses,$
               velocity_initial_guess[0],velocity_dispersion_initial_guess[0],0.,0.,/h3h4
mdap_create_starting_guesses,emission_line_kinematics_str,xbin_str,ybin_str,x2d,y2d,xbin_ems,ybin_ems,gas_kin_starting_guesses,$
               velocity_initial_guess[0],50.,0.,0.

wavelength_input=exp(log_wav_ems)
if mdap_spectral_fitting_version gt mdap_spectral_fitting_version_previous  or execute_all_modules eq 1 then begin
   mdap_spectral_fitting,log_spc_ems,log_err_ems,log_wav_ems,library_log_ems,log_wav_library_ems,velscale,$
        stellar_kinematics_ems,stellar_kinematics_ems_err,stellar_weights_ems,emission_line_kinematics_ems,emission_line_kinematics_ems_err,$
         emission_line_kinematics_ems_individual,emission_line_kinematics_ems_individual_err,$
        emission_line_intens_ems,emission_line_intens_ems_err,emission_line_fluxes_ems,emission_line_fluxes_ems_err,emission_line_equivW_ems,emission_line_equivW_ems_err,wavelength_input=wavelength_input,$
        wavelength_output_rest_frame_log,best_fit_model_ems,galaxy_minus_ems_fit_model_ems,best_template_ems,best_template_LOSVD_conv_ems,reddening_ems,reddening_ems_err,residuals_ems,$
        star_kin_starting_guesses=star_kin_starting_guesses,gas_kin_starting_guesses=gas_kin_starting_guesses,$
        emission_line_file=emission_line_file_spatial_binnin_3,$
        extra_inputs=spectra_fittin_parameters_patial_binning_3,/rest_frame_log,fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix,$
        mask_range=mask_range,external_library=external_library,/quiet;,$
;        extra_inputs=['MOMENTS=4','MDEGREE=4','DEGREE=-1','reddening=[0.01,0.01]','LAMBDA=exp(loglam_gal)']
   junk = size(emission_line_fluxes_ems) ;need to save all flux maps (warning: sgandalf computes intensities, not fluxes)
   for i = 0, junk[1]-1 do emission_line_intens_ems[i,*]=emission_line_intens_ems[i,*]/area_bins_ems[i]
   for i = 0, junk[1]-1 do emission_line_intens_ems_err[i,*]=emission_line_intens_ems_err[i,*]/sqrt(area_bins_ems[i])
   for i = 0, junk[1]-1 do emission_line_fluxes_ems[i,*]=emission_line_fluxes_ems[i,*]/area_bins_ems[i]
   for i = 0, junk[1]-1 do emission_line_fluxes_ems_err[i,*]=emission_line_fluxes_ems_err[i,*]/sqrt(area_bins_ems[i])

   execute_all_modules = 1


   ;calculate the real S/N of binned spectra
   bin_sn_ems_real = bin_sn_ems
   indxx = where(wavelength_output_rest_frame_log ge w_range_for_sn_computation_for_gas[0]  and wavelength_output_rest_frame_log le w_range_for_sn_computation_for_gas[1] )
   for i = 0, n_elements(bin_sn_ems_real) -1 do begin
       mdap_calculate_spectrum_sn,best_fit_model_ems[i,indxx],residuals_ems[i,indxx],wavelength_output_rest_frame_log[indxx],sn_per_angstorm,/rms
       bin_sn_ems_real[i] = sn_per_angstorm[0]
   endfor
   junk = temporary(sn_per_angstorm)
endif
printf,1,'[INFO] datacube '+root_name+' spectral fitting 3'
if save_intermediate_steps eq 1 then save,filename=root_name+mode+'_block4c.idl',/variables
;*********************************************************************************


;block5:
; BLOCK 5
;*** MEASURE OF THE EQUIVALENT WIDTH OF ABSORPTION LINE INDICES
;*******************


if n_elements(instrumental_fwhm_file) ne 0 then begin
   fwhm_instr=interpol(fwhm_ang,ww,wavelength_output_tpl) 
 endif else begin
    fwhm_instr=wavelength_output_tpl*0.+2.54
endelse

;resolution_x_lick=[4000.,4400.,4900.,5400.,6000.]   ;CONTROLLARE
;lick_fwhm_y=[11.5,9.2,8.4,8.4,9.8]                     ;CONTROLLARE
;lick_resolution_tmp=interpol(lick_fwhm_y,resolution_x_lick,wavelength_output_tpl)
;rrr=poly_fit(wavelength_output_tpl,lick_resolution_tmp,4,yfit=lick_resolution)
miles_resolution = fwhm_instr*0.+2.54
fwhm_diff_indices=sqrt(miles_resolution^2-fwhm_instr^2)   ;fwhm in angstrom
indici = where(finite(fwhm_diff_indices) NE 1)
if indici[0] ne -1 then fwhm_diff_indices[indici] = 0.


if mdap_measure_indices_version gt mdap_measure_indices_version_previous or execute_all_modules eq 1 then begin
   print, 'measuring indices on '+mdap_stc(n_elements(best_template_tpl[*,0]),/integer)+' spectra'
   mdap_measure_indices,absorption_line_indices,wavelength_output_tpl,galaxy_minus_ems_fit_model_tpl,$
             best_template_tpl,best_template_LOSVD_conv_tpl,stellar_kinematics_tpl[*,0],residuals_tpl,$
             fwhm_diff_indices,abs_line_indices,abs_line_indices_errors,abs_line_indices_template,abs_line_indices_template_losvd,$
             dir=output_dir;,/noplot
   sxaddpar,header_2d,'BLOCK5',mdap_measure_indices_version,'mdap_measure_indices version'
   execute_all_modules = 1

   

endif
printf,1,'[INFO] datacube '+root_name+' indices measured on '+mdap_stc(n_elements(best_template_tpl[*,0]),/integer)+' spectra'


if save_intermediate_steps eq 1 then save,filename=root_name+mode+'_block5.idl',/variables
;*********************************************************************************

;inizia_da_qui:
;restore,'ma001_145.125955+21.253809rss_block5.idl'
; BLOCK 6
;*** RADIAL PROFILE OF MEASURED QUANTITIES ***************************************
;  ;-- building error array from fit residuals 

if mdap_spatial_radial_binning_version gt mdap_spatial_radial_binning_version_previous or execute_all_modules eq 1 then begin
;  mdap_get_error_from_residual,residuals_ems,galaxy_minus_ems_fit_model_ems,input_errors
sz_ems=size(log_err_ems)
nl=n_elements(wavelength_output_rest_frame_log)
input_errors=fltarr(sz_ems[2],nl)
log_step_gal=LOG_WAV_EMS[1]-LOG_WAV_EMS[0]
for i = 0, sz_ems[2]-1 do begin
   rf_gal_lam = exp(log_wav_ems-stellar_kinematics_ems[i,0]/velscale*(log_step_gal))
   input_errors[i,*]=interpol(log_err_ems[*,i],rf_gal_lam,wavelength_output_rest_frame_log)
endfor

;-- radial binning 
;v0.2spatial_binning_scheme
stop
mdap_spatial_radial_binning,bin_sn_ems_real,x2d_reconstructed,y2d_reconstructed,spatial_binning_ems,xbin_ems,ybin_ems,ell,pa,$
      galaxy_minus_ems_fit_model_ems,input_errors,wavelength_output_rest_frame_log,$
      spatial_binning_rad,r_bin,r_bin_lo,r_bin_up,r2d_bin,r2d_bin_lo,r2d_bin_up,radially_binned_spectra,radially_binned_errors,$
      output_lrange=trim_wav_range_radial_binning,output_wav=output_wav,n_elements_bin=nelements_within_bin_radial,$
      low_radial_bins_user_inputs=low_radial_bins_user_inputs,upper_radial_bins_user_inputs=upper_radial_bins_user_inputs,$
      Reff=Reff,PSFsize=PSFsize,add_default_bins=add_default_bins

;--
   printf,1,'[INFO] datacube '+root_name+' radial binning: ',mdap_stc(n_elements(r_bin),/integer),' bins'

 
;-- spectral fit of radially binned spectra
   star_kin_starting_guesses_rbin = fltarr(n_elements(r_bin),4)
   gas_kin_starting_guesses_rbin  = fltarr(n_elements(r_bin),2)
   star_kin_starting_guesses_rbin[*,1]=velocity_dispersion_initial_guess[0] ;200.  ;str guess for sigma
   gas_kin_starting_guesses_rbin[*,1]=50.
   loglam_gal_rbin=alog(output_wav)
   
   mdap_spectral_fitting,radially_binned_spectra,radially_binned_errors,loglam_gal_rbin,library_log_ems,log_wav_library_ems,velscale,$
          stellar_kinematics_rbin,stellar_kinematics_rbin_err,stellar_weights_rbin,emission_line_kinematics_rbin,emission_line_kinematics_rbin_err,$
         emission_line_kinematics_rbin_individual,emission_line_kinematics_rbin_individual_err,$
          emission_line_intens_rbin,emission_line_intens_rbin_err,emission_line_fluxes_rbin,emission_line_fluxes_rbin_err,emission_line_equivW_rbin,emission_line_equivW_rbin_err,wavelength_input=output_wav,$
          wavelength_output_rbin,best_fit_model_rbin,galaxy_minus_ems_fit_model_rbin,best_template_rbin,best_template_LOSVD_conv_rbin,reddening_rbin,reddening_rbin_err,residuals_rbin,$
          star_kin_starting_guesses=star_kin_starting_guesses_rbin,gas_kin_starting_guesses=gas_kin_starting_guesses_rbin,$
          emission_line_file=emission_line_file_radial_binning,$
          extra_inputs=spectra_fittin_parameters_patial_binning_readial,fwhm_instr_kmsec_matrix=fwhm_instr_kmsec_matrix/3.,$
         ;extra_inputs=['MOMENTS=4','DEGREE=-1','mdegree=4','reddening=[0.01]','LAMBDA=exp(loglam_gal)'],$
         range_v_star=[-50.,50.],range_v_gas=[-50.,50.],mask_range=mask_range,external_library=external_library,/quiet ;,$
   printf,1,'[INFO] datacube '+root_name+' radial binning: spectral fitting' 
 
 
  ;calculate the real S/N of binned spectra
   bin_sn_rad_real = r_bin
   indxx = where(wavelength_output_rbin ge w_range_for_sn_computation[0]  and wavelength_output_rbin le w_range_for_sn_computation[1] )
   for i = 0, n_elements(bin_sn_rad_real) -1 do begin
       mdap_calculate_spectrum_sn,best_fit_model_rbin[i,indxx],residuals_rbin[i,indxx],wavelength_output_rbin[indxx],sn_per_angstorm,/rms
       bin_sn_rad_real[i] = sn_per_angstorm[0]
   endfor
   junk = temporary(sn_per_angstorm)
;--
;stop
;--measurement of abs indices
 ;  lick_resolution_tmp=interpol(lick_fwhm_y,resolution_x_lick,wavelength_output_rbin)
 ;  rrr=poly_fit(wavelength_output_rbin,lick_resolution_tmp,4,yfit=lick_resolution)
 ;  fwhm_diff_indices=sqrt(double(lick_resolution)^2.-double(fwhm_instr)^2.)*0. ;fwhm in angstrom

   if n_elements(instrumental_fwhm_file) ne 0 then begin
      fwhm_instr=interpol(fwhm_ang,ww,wavelength_output_rbin) 
   endif else begin
      fwhm_instr=wavelength_output_rbin*0.+2.54
   endelse
   miles_resolution = fwhm_instr*0.+2.54
   fwhm_diff_indices=sqrt(miles_resolution^2-fwhm_instr^2) ;fwhm in angstrom
   indici = where(finite(fwhm_diff_indices) NE 1)
   if indici[0] ne -1 then fwhm_diff_indices[indici] = 0.

   mdap_measure_indices,absorption_line_indices,wavelength_output_rbin,galaxy_minus_ems_fit_model_rbin,$
        best_template_rbin,best_template_LOSVD_conv_rbin,stellar_kinematics_rbin[*,0],residuals_rbin,$
        fwhm_diff_indices,abs_line_indices_rbin,$
        abs_line_indices_errors_rbin,abs_line_indices_template_rbin,abs_line_indices_template_losvd_rbin,dir=output_dir+'rbin_';,/noplot
;--
   sxaddpar,header_2d,'BLOCK6',mdap_spatial_radial_binning_version,'mdap_spatial_radial_binning_version'
   execute_all_modules = 1
   printf,1,'[INFO] datacube '+root_name+' radial binning: measured indices'
endif
if save_intermediate_steps eq 1 then save,filename=root_name+mode+'_block6a.idl',/variables


;- Rotation curves (kinemetry)

;stars
mdap_do_kinemetry,signal2d_reconstructed,x2d_reconstructed,y2d_reconstructed,xbin_str,ybin_str,$
                  STELLAR_KINEMATICS_STR[*,0],STELLAR_KINEMATICS_STR_ERR[*,0],$
                  PA_kin_str,PA_kin_std_str, q_kin_str,q_kin_std_str, Vsyst_str, Vsyst_std_str, Rad_kin_str, Vrot_str, Vrot_err_str,Vexp_str, Vexp_err_str;,$
                  ;Xcenter_used_for_stellar_rot_curve,Ycenter_used_for_stellar_rot_curve

;gas
mdap_do_kinemetry,signal2d_reconstructed,x2d_reconstructed,y2d_reconstructed,xbin_ems,ybin_ems,$
                  emission_line_kinematics_ems[*,0],emission_line_kinematics_EMS_err[*,0],$
                  PA_kin_ems, PA_kin_std_ems, q_kin_ems,q_kin_std_ems, Vsyst_ems,Vsyst_std_ems, Rad_kin_ems, Vrot_ems, Vrot_err_ems,Vexp_ems, Vexp_err_ems;,$
                  ;Xcenter_used_for_gas_rot_curve,Ycenter_used_for_gas_rot_curve

;- radial profiles (lambda, V/sigma, sigma?) -- still to test
if n_elements(radii_for_rprofiles) eq 0 then radii_rprofiles=Rad_kin_str else radii_rprofiles=radii_for_rprofiles
mdap_do_k_rprofiles,radii_rprofiles,STELLAR_KINEMATICS_STR[*,0]-Vsyst_str[0],STELLAR_KINEMATICS_STR_ERR[*,0],STELLAR_KINEMATICS_STR[*,1],STELLAR_KINEMATICS_STR_ERR[*,1],$
                    xbin_str,ybin_str,total(log_spc_str,1),ell,pa,$
                    lambda_profile,lambda_profile_err,vsigma_profile,vsigma_profile_err,sigma_profile,sigma_profile_err
 
;-
;*********************************************************************************
 
;-- SAVING HIGH LEVEL SCIENCE PRODUCTS RESULTS IN THE OUTPUT DATACUBE.... to be optimized
;   


NINDICES = n_elements(ABS_LINE_INDICES[0,*])
mdap_read_indices_definitions,absorption_line_indices,indices=indices
readcol,emission_line_file_spatial_binnin_3,cnt,ln_name,ln_wav,comment='#',format='I, A, A',/silent
NLINES = n_elements(cnt)

writefits,output_filefits,signal2d_reconstructed,header_2d
k=1
mdap_add_fits_layer,output_filefits,noise2d_reconstructed,k,'EXTNAME','noise'
k=k+1
mdap_add_fits_layer,output_filefits,spatial_binning_tpl,k,'EXTNAME','Binning_map_1'
k=k+1
 stringa = ["X","Y","area_bin","StoN","Nelements"]
 for i = 0, NINDICES-1 do begin 
    stringa = [stringa,indices[i].name,indices[i].name+"_err"]
 endfor 
 stringa2="{"
 for i = 0, n_elements(stringa)-2 do  stringa2 = stringa2+stringa[i]+':0.,'
 stringa2 = stringa2+stringa[i]+':0.}'
 d = execute('str='+stringa2)
 p1=replicate(str,n_elements(xbin_tpl))
 p1.x=xbin_tpl
 p1.y=ybin_tpl
 p1.area_bin = AREA_BINS_TPL
; p1.reddening = reddening_tpl
 p1.StoN = bin_sn_tpl_real
 p1.Nelements=nelements_within_bin_tpl
 for i = 0, NINDICES-1 do begin 
    d=execute('p1.'+indices[i].name+'=ABS_LINE_INDICES[*,i]')
    d=execute('p1.'+indices[i].name+'_err=ABS_LINE_INDICES_ERRORS[*,i]')
 endfor 
 mwrfits,p1,output_filefits,/silent
 h1=HEADFITS(output_filefits,EXTEN = k)                                         
 sxaddpar,h1,'EXTNAME','Binning_1_data'
 modfits,output_filefits,0,h1,exten_no=k

; storing stellar kinematics
k=k+1
mdap_add_fits_layer,output_filefits,spatial_binning_str,k,'EXTNAME','Binning_map_2'
k=k+1
stringa = ["X","Y","area_bin","StoN","Nelements","Vel","Vel_err","Disp","Disp_err","H3","H3_err","H4","H4_err","Chi2_DOF","Vel_2moms","Vel_2moms_err","Disp_2moms","Disp_2moms_err"]
stringa2="{"
 for i = 0, n_elements(stringa)-2 do  stringa2 = stringa2+stringa[i]+':0.,'
 stringa2 = stringa2+stringa[i]+':0.}'
 d = execute('str='+stringa2)
 p2=replicate(str,n_elements(xbin_str))
 p2.x=xbin_str
 p2.y=ybin_str
 p2.area_bin = AREA_BINS_STR
 p2.StoN = bin_sn_str_real
 p2.Nelements=nelements_within_bin_str
 p2.Vel=STELLAR_KINEMATICS_STR[*,0]
 p2.Vel_err=STELLAR_KINEMATICS_STR_ERR[*,0]
 p2.Disp=STELLAR_KINEMATICS_STR[*,1]
 p2.Disp_err=STELLAR_KINEMATICS_STR_ERR[*,1]
 p2.H3=STELLAR_KINEMATICS_STR[*,2]
 p2.H3_err=STELLAR_KINEMATICS_STR_ERR[*,2]
 p2.H4=STELLAR_KINEMATICS_STR[*,3]
 p2.H4_err=STELLAR_KINEMATICS_STR_ERR[*,3]
 p2.Chi2_DOF=STELLAR_KINEMATICS_STR[*,4]

 if n_elements(spectra_fittin_parameters_patial_binning_2_two_moments) eq 0 then begin ; if the fit with only 2 moments is not required, I copy the results with 4 moments
     p2.Vel_2moms=STELLAR_KINEMATICS_STR[*,0]
     p2.Vel_2moms_err=STELLAR_KINEMATICS_STR_ERR[*,0]
     p2.Disp_2moms=STELLAR_KINEMATICS_STR[*,1]
     p2.Disp_2moms_err=STELLAR_KINEMATICS_STR_ERR[*,1]
 endif else begin                                                                      ; ... otherwise, I save also the results with the fit considering only 2 moments
     p2.Vel_2moms=STELLAR_KINEMATICS_STR_TWO_MOMENTS[*,0]
     p2.Vel_2moms_err=STELLAR_KINEMATICS_STR_TWO_MOMENTS_ERR[*,0]
     p2.Disp_2moms=STELLAR_KINEMATICS_STR_TWO_MOMENTS[*,1]
     p2.Disp_2moms_err=STELLAR_KINEMATICS_STR_TWO_MOMENTS_ERR[*,1]
 endelse


 mwrfits,p2,output_filefits,/silent
 h1=HEADFITS(output_filefits,EXTEN = k)                                         
 sxaddpar,h1,'EXTNAME','Binning_2_data'
 modfits,output_filefits,0,h1,exten_no=k

; storing emission lines kinematics, fluxes, and equivalent widths
k=k+1
mdap_add_fits_layer,output_filefits,spatial_binning_ems,k,'EXTNAME','Binning_map_3'
k=k+1
 stringa = ["X","Y","area_bin","StoN","Nelements","Vel","Vel_err","Disp","Disp_err","reddening_star","reddening_star_err","reddening_gas","reddening_gas_err","Chi2_DOF"]
 for i = 0, NLINES-1 do begin 
    stringa = [stringa,ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_intens',ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_intens_err']
 endfor 
 for i = 0, NLINES-1 do begin 
    stringa = [stringa,ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_flux',ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_flux_err']
 endfor 
 for i = 0, NLINES-1 do begin 
    stringa = [stringa,ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_EW',ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_EW_err']
 endfor 
 for i = 0, NLINES-1 do begin 
    stringa = [stringa,ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Vel',ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Vel_err']
 endfor 
 for i = 0, NLINES-1 do begin 
    stringa = [stringa,ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Sig',ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Sig_err']
 endfor 

 stringa2="{"
 for i = 0, n_elements(stringa)-2 do  stringa2 = stringa2+stringa[i]+':0.,'
 stringa2 = stringa2+stringa[i]+':0.}'
 d = execute('str='+stringa2)
 p3=replicate(str,n_elements(xbin_ems))
 p3.x=xbin_ems
 p3.y=ybin_ems
 p3.area_bin = AREA_BINS_EMS
 p3.StoN = bin_sn_ems_real
 p3.Nelements=nelements_within_bin_ems
 p3.Vel=EMISSION_LINE_KINEMATICS_EMS[*,0]
 p3.Vel_err=EMISSION_LINE_KINEMATICS_EMS_ERR[*,0]
 p3.Disp=EMISSION_LINE_KINEMATICS_EMS[*,1]
 p3.Disp_err=EMISSION_LINE_KINEMATICS_EMS_ERR[*,1]
 p3.reddening_star = reddening_ems[*,0]
 p3.reddening_star_err = reddening_ems_err[*,0]
 p3.reddening_gas = reddening_ems[*,1]
 p3.reddening_gas_err = reddening_ems_err[*,1]
 p3.Chi2_DOF=STELLAR_KINEMATICS_EMS[*,4]


 for i = 0, NLINES-1 do begin 
    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_intens=EMISSION_LINE_INTENS_EMS[*,i]')
    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_intens_err=EMISSION_LINE_INTENS_EMS_ERR[*,i]')
 endfor 
 for i = 0, NLINES-1 do begin 
    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_flux=EMISSION_LINE_FLUXES_EMS[*,i]')
    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_flux_err=EMISSION_LINE_FLUXES_EMS_ERR[*,i]')
 endfor 
 for i = 0, NLINES-1 do begin 
    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_EW=EMISSION_LINE_EQUIVW_EMS[*,i]')
    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_EW_err=EMISSION_LINE_EQUIVW_EMS_ERR[*,i]')
 endfor 
 for i = 0, NLINES-1 do begin 
    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Vel=emission_line_kinematics_ems_individual[*,i,0]')
    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_VEl_err=emission_line_kinematics_ems_individual_err[*,i,0]')
 endfor 
 for i = 0, NLINES-1 do begin 
    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Sig=emission_line_kinematics_ems_individual[*,i,1]')
    d=execute('p3.'+ln_name[i]+'_'+mdap_stc(round(float(ln_wav[i])))+'_Sig_err=emission_line_kinematics_ems_individual_err[*,i,1]')
 endfor 
 mwrfits,p3,output_filefits,/silent
 h1=HEADFITS(output_filefits,EXTEN = k)                                         
 sxaddpar,h1,'EXTNAME','Binning_3_data'
 modfits,output_filefits,0,h1,exten_no=k

;storing  radial binning
k=k+1
mdap_add_fits_layer,output_filefits,spatial_binning_rad,k,'EXTNAME','Binning_map_radial'
k=k+1
stringa = ["amaj","amaj_lo","amaj_up","StoN","Nelements","Disp","Disp_err","Chi2_DOF"]
for i = 0, NINDICES-1 do begin 
   stringa = [stringa,indices[i].name+"_rad",indices[i].name+"_rad_err"]
endfor 
stringa2="{"
 for i = 0, n_elements(stringa)-2 do  stringa2 = stringa2+stringa[i]+':0.,'
 stringa2 = stringa2+stringa[i]+':0.}'
 d = execute('str='+stringa2)
 p4=replicate(str,n_elements(r_bin))
 p4.amaj=r_bin
 p4.amaj_lo=r_bin_lo
 p4.amaj_up=r_bin_up
 p4.StoN=bin_sn_rad_real
 p4.Nelements=nelements_within_bin_radial
 p4.Disp=stellar_kinematics_rbin[*,1]
 p4.Disp_err=stellar_kinematics_rbin_err[*,1]
 p4.Chi2_DOF=stellar_kinematics_rbin[*,4]

 for i = 0, NINDICES-1 do begin 
    d=execute('p4.'+indices[i].name+'_rad=ABS_LINE_INDICES_RBIN[*,i]')
    d=execute('p4.'+indices[i].name+'_rad_err=ABS_LINE_INDICES_ERRORS_RBIN[*,i]')
 endfor 
 mwrfits,p4,output_filefits,/silent
 h1=HEADFITS(output_filefits,EXTEN = k)                                         
 sxaddpar,h1,'EXTNAME','Binning_radial_data'
 modfits,output_filefits,0,h1,exten_no=k

;storing rotation curve of stars
k=k+1
;stringa = ["a_rot_curve","PA_kin_stars","PA_kin_std_stars","q_kin_stars","q_kin_std_stars","Vsyst_stars","Vsyst_std_stars", "Vrot_stars", "Vrot_err_stars","Vexp_stars", "Vexp_err_stars","x0","y0"]
stringa = ["a_rot_curve","PA_kin_stars","PA_kin_std_stars","q_kin_stars","q_kin_std_stars","Vsyst_stars","Vsyst_std_stars", "Vrot_stars", "Vrot_err_stars","Vexp_stars", "Vexp_err_stars"]
stringa2="{"
for i = 0, n_elements(stringa)-2 do  stringa2 = stringa2+stringa[i]+':0.,'
stringa2 = stringa2+stringa[i]+':0.}'
 d = execute('str='+stringa2)
p5=replicate(str,n_elements(Rad_kin_str))
p5.a_rot_curve=Rad_kin_str
p5.PA_kin_stars=replicate(PA_kin_str,n_elements(Rad_kin_str))
p5.PA_kin_std_stars=replicate(PA_kin_std_str,n_elements(Rad_kin_str))
p5.q_kin_stars=replicate(q_kin_str,n_elements(Rad_kin_str))
p5.q_kin_std_stars=replicate(q_kin_std_str,n_elements(Rad_kin_str))
p5.Vsyst_stars=Vsyst_str
p5.Vsyst_std_stars=Vsyst_std_str
p5.Vrot_stars= Vrot_str 
p5.Vrot_err_stars=Vrot_err_str
p5.Vexp_stars=Vexp_str
p5.Vexp_err_stars=Vexp_err_str
;p5.x0=replicate(Xcenter_used_for_stellar_rot_curve,n_elements(Rad_kin_str)) 
;p5.y0=replicate(Ycenter_used_for_stellar_rot_curve,n_elements(Rad_kin_str)) 
mwrfits,p5,output_filefits,/silent
h1=HEADFITS(output_filefits,EXTEN = k)                                         
sxaddpar,h1,'EXTNAME','Stars Rotation'
modfits,output_filefits,0,h1,exten_no=k


;storing rotation curve of gas
k=k+1
;stringa = ["a_rot_curve","PA_kin_gas","PA_kin_std_gas","q_kin_gas","q_kin_std_gas","Vsyst_gas","Vsyst_std_gas", "Vrot_gas", "Vrot_err_gas","Vexp_gas", "Vexp_err_gas","x0","y0"]
stringa = ["a_rot_curve","PA_kin_gas","PA_kin_std_gas","q_kin_gas","q_kin_std_gas","Vsyst_gas","Vsyst_std_gas", "Vrot_gas", "Vrot_err_gas","Vexp_gas", "Vexp_err_gas"]
stringa2="{"
for i = 0, n_elements(stringa)-2 do  stringa2 = stringa2+stringa[i]+':0.,'
stringa2 = stringa2+stringa[i]+':0.}'
 d = execute('str='+stringa2)
p6=replicate(str,n_elements(Rad_kin_ems))
p6.a_rot_curve=Rad_kin_ems
p6.PA_kin_gas=replicate(PA_kin_ems,n_elements(Rad_kin_ems))
p6.PA_kin_std_gas=replicate(PA_kin_std_ems,n_elements(Rad_kin_ems))
p6.q_kin_gas=replicate(q_kin_ems,n_elements(Rad_kin_ems))
p6.q_kin_std_gas=replicate(q_kin_std_ems,n_elements(Rad_kin_ems))
p6.Vsyst_gas=Vsyst_ems
p6.Vsyst_std_gas=Vsyst_std_ems
p6.Vrot_gas= Vrot_ems 
p6.Vrot_err_gas=Vrot_err_ems
p6.Vexp_gas=Vexp_ems
p6.Vexp_err_gas=Vexp_err_ems
;p6.x0=replicate(Xcenter_used_for_gas_rot_curve,n_elements(Rad_kin_ems)) 
;p6.y0=replicate(Ycenter_used_for_gas_rot_curve,n_elements(Rad_kin_ems)) 
mwrfits,p6,output_filefits,/silent
h1=HEADFITS(output_filefits,EXTEN = k)                                         
sxaddpar,h1,'EXTNAME','GAS rotation'
modfits,output_filefits,0,h1,exten_no=k

;storing lambda, v/sigma and sigma radial profiles
k = k+1
stringa = ["radii_rprofiles","lambda_profile","lambda_profile_err","vsigma_profile","vsigma_profile_err","sigma_profile","sigma_profile_err"]
stringa2="{"
for i =0, n_elements(stringa)-2 do stringa2=stringa2+stringa[i]+':0.,'
stringa2 = stringa2+stringa[i]+':0.}'
d = execute('str='+stringa2)
p7=replicate(str,n_elements(radii_rprofiles))
p7.radii_rprofiles=radii_rprofiles
p7.lambda_profile=lambda_profile
p7.lambda_profile_err=lambda_profile_err
p7.vsigma_profile=vsigma_profile
p7.vsigma_profile_err=vsigma_profile_err
p7.sigma_profile=sigma_profile
p7.sigma_profile_err=sigma_profile_err
mwrfits,p7,output_filefits,/silent
h1=HEADFITS(output_filefits,EXTEN = k)                                         
sxaddpar,h1,"EXTNAME","lambda, V/S, S"
modfits,output_filefits,0,h1,exten_no=k


printf,1,'[INFO] '+root_name+' output saved: '+output_filefits
save,/variables,filename=output_idlsession; output_dir+'mdap_session.idl'
printf,1,'[INFO] '+root_name+' output saved: '+output_idlsession

printf,1,root_name+'_'+mode+' execution time '+mdap_stc(systime(/seconds)/60./60.-t0)+' hours'
close, 1
openw,1,root_name+'_'+mode+'.done'
printf,1,root_name+'_'+mode+' execution time '+mdap_stc(systime(/seconds)/60./60.-t0)+' hours'
close, 1

;*********************************************************************************
;
;    
;   
;   ;
;    ##############################################
;    #####    MODEL DEPENDENT DATA PRODUCTS   #####
;    ##############################################
;
;T.B.D.
end

; ;uncomment these lines when performing parallelization on SCIAMA.
; PREF_SET,'IDL_CPU_TPOOL_NTHREADS',1,/COMMIT
; manga_dap,!!num!!
; EXIT
; END
