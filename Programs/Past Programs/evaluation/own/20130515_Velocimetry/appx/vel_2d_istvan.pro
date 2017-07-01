pro Vel_2d_Istvan

;file_name = 'G:\CSDXdata\camera data\2012_10_19\15163(m=1,1200G)\15163.cine' ;
    file_name = 'G:\CSDXdata\camera data\12216(1000G)\12216.cine' ;
   codec = 'CVID'    ; This codec works in WMP and Quicktime

   window_width = 4096; 4096
  
;   for mx=-2,2 do begin 
  fmin = -17000
;  fmin = -17000+ mx*1000
  fmax = fmin + window_width
   
  flow=2
  fhigh=50; KHZ 
;------- read in camera data-----------------------------------
  image3d = cine_load_1(file_name = file_name, fmin = fmin, fmax = fmax, MOV)
  n_frames = MOV.nframes
  frame_rate = MOV.pps
  nx = MOV.nx
  ny = MOV.ny 
     
  xcent = 65;60;63;60/61 ;72 ; movie frame size, change plasma center
  ycent = 63;31;68; 33
  cm_per_pixel = 0.5/4   ;  cm 1.0/6.907; 1.2/8.0
    
  sample_rate = 210526; 1/4.7e-6 ; frames/sec
  
 ; ---------- subtract the mean from image --------------

  mean_image = fltarr(MOV.nx,MOV.ny)

  sum = TOTAL(image3d[*,*,0:n_frames-1], 3)
  mean_image = sum/ FLOAT(n_frames)


; image_sub = fltarr (64, 16, n_frames)

camtime = findgen(n_frames)/sample_rate
image_sub_4 = fltarr (65, 2, n_frames) ; 
image_sub_6 = fltarr (65, 2, n_frames) ; 
image_sub_8 = fltarr (65, 2, n_frames) ;
image_sub_10 = fltarr (65, 2, n_frames) ;
  camx_4 =  fltarr (65,2)+1
  camy_4 =  fltarr (65,2)+1
  camx_6 =  fltarr (65,2)+1
  camy_6 =  fltarr (65,2)+1
  camx_8 =  fltarr (65,2)+1
  camy_8 =  fltarr (65,2)+1
   camx_10 =  fltarr (65,2)+1
  camy_10 =  fltarr (65,2)+1
;--------隔4个点----------------------------
 for index = 0,1 do begin
     col_y_4= 61 + 4 * index
;   image_sub_4[*,index,*] = image3d[65:126,col_y_4,*] - mean_image[65:126,col_y_4, indgen(n_frames)]; angle=0     
    image_sub_4[*,index,*] = image3d[0:64,col_y_4,*] - mean_image[0:64,col_y_4, indgen(n_frames)]  ; angle=180
 endfor
   camframes_4=image_sub_4[*,*,*]
 for j=0,1 do begin
   camy_4[0:64,j] = j*cm_per_pixel*4 ; 隔4个点
  endfor
 for m = 0,64 do begin
   camx_4[m,0:1] = m*cm_per_pixel*4
 endfor
 
 ;--------------------end--------------------
 ;--------隔6个点----------------------------
 for index = 0,1 do begin
     col_y_6= 60 + 6 * index
;   image_sub_6[*,index,*] = image3d[65:126,col_y_6,*] - mean_image[65:126,col_y_6, indgen(n_frames)]  ; angle=0   
    image_sub_6[*,index,*] = image3d[0:64,col_y_6,*] - mean_image[0:64,col_y_6, indgen(n_frames)]  ; angle=90  
 endfor
   camframes_6=image_sub_6[*,*,*]
 for j=0,1 do begin
   camy_6[0:64,j] = j*cm_per_pixel*6 ; 隔6个点
  endfor

 for m = 0,64 do begin
     camx_6[m,0:1] = m*cm_per_pixel*6
 endfor
 
 ;--------------------end---------------------
 ;--------隔8个点----------------------------
 for index = 0,1 do begin
     col_y_8= 59 + 8 * index
;     image_sub_8[*,index,*] = image3d[65:126,col_y_8,*] - mean_image[65:126,col_y_8, indgen(n_frames)] ; angle=0 
     image_sub_8[*,index,*] = image3d[0:64,col_y_8,*] - mean_image[0:64,col_y_8, indgen(n_frames)] ; angle=90    
 endfor
   camframes_8=image_sub_8[*,*,*]
 for j=0,1 do begin
   camy_8[0:64,j] = j*cm_per_pixel*8 ; 隔8个点
  endfor

 for m = 0,64 do begin
     camx_8[m,0:1] = m*cm_per_pixel*8
 endfor
 
 ;--------------------end---------------------
 ;
 ; ;--------隔10个点----------------------------
 for index = 0,1 do begin
     col_y_10= 59 + 10 * index
;     image_sub_8[*,index,*] = image3d[65:126,col_y_8,*] - mean_image[65:126,col_y_8, indgen(n_frames)] ; angle=0 
     image_sub_10[*,index,*] = image3d[0:64,col_y_10,*] - mean_image[0:64,col_y_10, indgen(n_frames)] ; angle=90    
 endfor
   camframes_10=image_sub_10[*,*,*]
 for j=0,1 do begin
   camy_10[0:64,j] = j*cm_per_pixel*10 ; 隔8个点
  endfor

 for m = 0,64 do begin
     camx_10[m,0:1] = m*cm_per_pixel*10
 endfor
 
 ;--------------------end---------------------
;  camframes=image_sub[*,*,*]
;  camtime = findgen(n_frames)/sample_rate
;  camx=  fltarr (64,16)+1
;  camy=  fltarr (64,16)+1
  

; for j=0,15 do begin
;   camy[0:63,j] = j*cm_per_pixel
;  endfor
;
; for m = 0,63 do begin
;     camx[m,0:15] = m*cm_per_pixel
; endfor
  
  tstep=50
       vel_4=vel_2d_t($
                   frames=camframes_4,$
                   x=camx_4,$
                   y=camy_4,$
                   ncor =100,$
                   tstep=tstep,$
                   time=camtime,$
                   correlation=CORRELATION,$
                   flow=flow,$ ;khz
                   fhigh=fhigh)
                   
      vel_6=vel_2d_t($
                   frames=camframes_6,$
                   x=camx_6,$
                   y=camy_6,$
                   ncor =100,$
                   tstep=tstep,$
                   time=camtime,$
                   correlation=CORRELATION,$
                   flow=flow,$ ;khz
                   fhigh=fhigh)
      vel_8=vel_2d_t($
                   frames=camframes_8,$
                   x=camx_8,$
                   y=camy_8,$
                   ncor =100,$
                   tstep=tstep,$
                   time=camtime,$
                   correlation=CORRELATION,$
                   flow=flow,$ ;khz
                   fhigh=fhigh)
;     vel_10=vel_2d_t($
;                   frames=camframes_8,$
;                   x=camx_10,$
;                   y=camy_10,$
;                   ncor =100,$
;                   tstep=tstep,$
;                   time=camtime,$
;                   correlation=CORRELATION,$
;                   flow=flow,$ ;khz
;                   fhigh=fhigh)


numb= long(window_width/tstep)
vel_vy4=total(vel_4.vy,1)/numb
vel_vy6=total(vel_6.vy,1)/numb
vel_vy8=total(vel_8.vy,1)/numb
;vel_vy10=total(vel_10.vy,1)/numb

;vel_vy2=findgen(63)

;for idx=0,62 do begin
;   vel_vy2[idx] = vel_vy1[15*idx+7] 
; endfor

  num=62 
x_pos = (findgen(63))*cm_per_pixel

  WINDOW, 3 , XSIZE=(960), YSIZE=(750)
        PLOT, x_pos[0:num],vel_vy8[0:num], $
              ytitle='!12< !3V!S!D!7h!3!R !12>!3   (cm/s)', xtitle='r (cm)',$
                  title='V_theta', $
                  THICK = 2.8, CHARTHICK = 2.4, CHARSIZE = 2.2, COLOR = 0, PSYM = -6, $
                  xrange = [0,10];, yrange=[-2e4, 1.2e5]
      oplot, x_pos[0:num-1],vel_vy6[0:num-1],$
      THICK = 2.8 ,COLOR = 50, PSYM = -6          
       oplot, x_pos[0:num],vel_vy4[0:num],$
      THICK = 2.8, COLOR = 100, PSYM = -6 
;        oplot, x_pos[0:num],vel_vy10[0:num],$
;      THICK = 2.8, COLOR = 100, PSYM = -6 
  
;endfor


     position_start = 0
     position_end = num

outfile = 'G:\CSDXdata\camera data\12216(1000G)\istvan’s\12216_180_2_50KHZ_1.dat'
 openw, lun, outfile, /get_lun ;, width = 500
  printf, lun 
  for i = position_start, position_end do begin 
    printf, lun, x_pos[i], vel_vy4[i],vel_vy6[i],vel_vy8[i],$
                  format='(f20.10,e20.10,e20.10,e20.10)'
  endfor
free_lun, lun  

 
  end