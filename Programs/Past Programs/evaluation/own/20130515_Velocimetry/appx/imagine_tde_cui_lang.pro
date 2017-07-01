;=======================================================================================================
;Name: image_filtering_2D
;Author: Lang Cui (modified from Zheng's code)
;Update: 28Jan2009
;---------------------------------------------------------
;
;Description: 1) filtering every pixel of the fast movie, so that the mode coupling could be visualized
;       directly.
;             2) calculate v_theta at specified radius
;=======================================================================================================

@movie_io
PRO imagine_TDE_cui_lang

  DEVICE,DEC=0   ; For indexed image

   file_name = 'G:\CSDXdata\2012_03_23\3.6mTorr_2kW_260A.cine' ;
   codec = 'CVID'    ; This codec works in WMP and Quicktime
   default, fileout, 'G:\CSDXdata\2012_03_23\3.6mTorr_2kW_260A.avi'

   fmin = 0
   fmax = 500


  image3d = cine_load_1(file_name = file_name, fmin = fmin, fmax = fmax, MOV)
  n_frames = MOV.nframes
  frame_rate = MOV.pps
  nx = MOV.nx
  ny = MOV.ny

  ; Time axis
  frame_num = findgen(n_frames)
  time = ( frame_num + fmin )/frame_rate*1000      ; unit is ms
  
  time_axis_position = ( frame_num + fmin )/frame_rate*1000*(float( indgen(n_frames) ))   


; ---------- subtract the mean from image --------------

  mean_image = fltarr(MOV.nx,MOV.ny)

  sum = TOTAL(image3d[*,*,0:n_frames-1], 3)
  mean_image = sum/ FLOAT(n_frames)

  image_sub = fltarr (MOV.nx, MOV.ny, n_frames)
 
  image_sub[*,*,*] = image3d[*,*,*] - mean_image[*,*, indgen(n_frames)] ; subtract the mean from image


;;  ---------------------------------- filtering signals ----------------------------------

  image_scale = SIZE(image3d, /dimensions)
;  image3d_filtered = fltarr(image_scale[0],image_scale[1],image_scale[2])   ;  no need to filter for TDE
;  image3d_filtered_low_band = fltarr(image_scale[0],image_scale[1],image_scale[2])
;  image3d_filtered_high_band = fltarr(image_scale[0],image_scale[1],image_scale[2])

  fNyq = frame_rate/2.0
  
  
; xcent = 63;60 ;60/61 ;72 ; movie frame size, change plasma center
; ycent = 68;31; 33
  xcent = 60;63;60/61 ;72 ; movie frame size, change plasma center
  ycent = 31;68; 33
  cm_per_pixel = 1.0/6.907; 1.2/8.0
  
  horizon_position = float( indgen(nx) - xcent ) * cm_per_pixel
  vertical_position = float( indgen(ny) - ycent ) * cm_per_pixel
  
  frame_start = 300;352         ;  the index num of the first frame
  frame_end = 500;550       ;  the index num of the final frame

  radius_start = LONG(3.0/cm_per_pixel)        ;  the first radius (in unit of pixel) of the band
  radius_end = LONG(3.6/cm_per_pixel)       ;  the final radius
  radius_step = radius_end - radius_start
  
  theta_step = !Pi/45.0
  ntheta=128        ; total num of points after interpolation
  
  bandpass_low = 3   
  bandpass_high = 30
  window_width = 4096.0
  sample_rate = 0.500e6
  bandpass_low_sub = LONG( 1.0e3*bandpass_low * window_width/sample_rate )
  bandpass_high_sub = LONG( 1.0e3*bandpass_high * window_width/sample_rate )
      
  V_theta = FLTARR(n_frames)     
        
                                               
;     +++++++++++++++++ time-averaged azimuthal velocity ++++++++++++++++
                        
                         
      select_on = 1 ;switch for on/off of indicating the maximum density gradient
        if select_on EQ 1 then begin
         for i_frame=frame_start, frame_end - 1L do begin
      
              data=image_sub[*,*, i_frame] ; frame        

               pixel_select, $
          ;                 -------- input variables --------
                            rad=radius_start, $ ; interested radius
                            theta_step=theta_step, $ ;step moving on circle
                            xcent=xcent, $ ; units: pixel   x coordinate of the center of image
                            ycent=ycent, $
                            data=data,$; 2D frame data
                            ntheta=ntheta, $ ; total points after interpolation
                            
          ;                  --------- output variables ----------
                            totarr=totarr, $ ; 2D frame data with the circle indicated
                            d_choose=d_choose, $ ; data on the circle w/o removing multiple counting
                            result_interp=result_interp, $ ; interpolated data on the circle (w/ multiple counting removed)
                            thetagrid=thetagrid, $; the evenly spaced theta grid for interpolation
                            realtheta=realtheta, $; the real theta for data with unique coordinates
                            sub = sub,$;1D index on the circle corresponding to uni_data
                            uni_data = uni_data, $;data with unique coordinates, but w/ zero at the end of vector
                            uni_real_data= uni_real_data, $; unique data vector with zeros at the end cut
                            n_pixel = n_pixel,$; total number of unique pixels on the circle
                            xx = xx, yy = yy ; the x and y coordinate for pixels on the circle (possible multiple counting)

      
      
              image_sub[*,*, i_frame] = totarr ; image with the circle on
              
              data_scale = SIZE(uni_data, /dimensions)
              uni_data_reformed = fltarr(data_scale[0],float(frame_end - frame_start)) 
            
              
              V_theta[i_frame] += CROSS_POWER_WEIGHTED_PHASE_DIFFERENCE_VELOCITY_imagine( $
                                 x = uni_data_reformed[34, *], $   ;the first vector used to calculate cross-power,and it will be conjugated
                                                        ; when caculating the cross-power
                                 y = uni_data_reformed[44, *], $   ;the second vector used to calculate cross-power
                                 window_width = window_width, $  ; window width used when calculating the cross-power
                                 /use_hanning, $ ; keyword to determine if using a hanning window when
                                                              ; when calculating the cross power
                                ;/remove_mean, $ ; keyword to determine if removing the mean before calculate
                                                               ; the cross-power
                      
                                  pixel_spacing = 0.25, $ ; spacing between the two interested tips
                                  sample_rate = sample_rate, $  ; sample rate of the data
                                  bandpass_low = 5.0, $ ; the lowest frequency for doing the average
                                                            ; over frequencies  ( unit: KHz )
                                  bandpass_high = 20.0, $ ; the highestest frequency for doing the average
                                                          ; over frequencies  ( unit: KHz )
                      
                                  count_unphysical = count_unphysical_temp, $ ; counting the number of cases which does not
                                                                    ; fit the physical criterion
                                  k_low_cutoff = 0.1,$  ; should use a unit corresponding to the tip_spacing
                                                          ; unit ( 'cm^-1' or 'm^-1' )
                                                          ; All the k with absolute value above this value will be considered
                                                          ; as a unphysical noise spike and will be ignored,
                                  k_high_cutoff = 12, $
                                  v_theta_vector_cross_power_weighted = v_theta_vector_cross_power_weighted_temp , $
                                                                                    ; vector contains velocity for each window
                                                                                          
                                  velocity_freq_resol = velocity_freq_resol_temp, $;frequency resolved velocity, corresponding to
                                                                                          ;frequency range from bandpass_low to 
                                                                                          ;bandpass_high
                                                                                          ; crosspower weighted average over windows
                                   frequency_vector_plot = frequency_vector_plot, $;for plotting frequency resolved velocity
                                                                                           ;runs from bandpass_low to bandpass_high
                                  /frequency_resolved     )
              
          endfor; i_frame
          
         endif
  ;      +++++++++++++++++++ end of indicating the maximal density gradient ++++++++++++++++++++++            

    
;


;;        ++++++++++++++++++++ get a solid spot for drawing +++++++++++++++++++
;          radius_density_grad = LONG(3.0/cm_per_pixel)
;          radius_shear_layer = LONG(4.1/cm_per_pixel)
;          
;          density_grad_solid_spot = fltarr(image_scale[0], image_scale[1], fmax )
;          shear_layer_solid_spot = fltarr(image_scale[0], image_scale[1], fmax )
;          
;          FOR i_row = 0, image_scale[0] - 1 DO BEGIN
;                 FOR i_col = 0, image_scale[1] - 1 DO BEGIN
;            
;                     distance_temp = round( ( ( float(i_row - xcent) )^2 + ( float(i_col - ycent) )^2 )^0.5 )
;                     if distance_temp LE radius_density_grad then begin
;                        density_grad_solid_spot[i_row, i_col, *] = 255
;                     endif
;                     
;                     if distance_temp LE radius_shear_layer then begin
;                        shear_layer_solid_spot[i_row, i_col, *] = 255
;                     endif                                   
;                                  
;                 ENDFOR ;i_col
;          ENDFOR; i_row
 
 WINDOW, 5, XSIZE=(960), YSIZE=(720)
        PLOT,   time_axis_position[frame_start:frame_end], $
              V_theta[frame_start:frame_end] , $
              ytitle='!12< !3V!S!D!7h!3!R !12>!3   (cm/s)',xtitle='r (cm)',$
                  title='average azimuthal velocity', $
                  THICK = 2.8, CHARTHICK = 2.4, CHARSIZE = 2.6, COLOR = 0, PSYM = -6, $
                  xrange = [0,10];, yrange=[-2e4, 1.2e5]

end