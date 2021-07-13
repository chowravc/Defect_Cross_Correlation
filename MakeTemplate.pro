;********************************************************************************************************************************
; MakeTemplate
;******************************************************************************************************************************
; DESCRIPTION
    ; Generate a 2D image template used to track point defects
    ; in freely suspended Smectic C liquid crystal films.
    ; These defects resemble Maltese crosses with 4 brushes (with crossed polarizers)
    ; or simple crosses with 2 brushes (with decrossed polarizer/analyzer),
    ; The +1 and -1 defects
    ; may have different orientations depending on the polarizer/analyzer settings.
    ; The template is cross-correlated with (part of) some image to identify
    ; the location of the defect.

; NOTES
    ; Original template procedure created by Dave Coleman ca. 2002.
    ; Modified by Joe Maclennan August 2006 to produce rotated cross if desired.
    ; Changed to be a Function and to display zoomed in template images for easier viewing.
    ; Modified to allow no. of brushes in the "cross" to be passed to the procedure.
    ; Added local common block to keep track of graphics window indexes (January 2013).

FUNCTION MakeTemplate, Ra=ra, Phi0=phi0, Nbrushes=nbrushes, Tvflag=tvflag
    ; ra is the "radius" of the square template (i.e. half the width)
    ; phi0 is the orientation of the first arm of the cross (in degrees)
    ; nbrushes is the number of bright brushes in the template
    ; tvflag determines whether the final template is displayed

; define a common block in order to save some local variables
    COMMON mtlocal, windex0, windex1, windex2, windex3

; Check arguments
    IF N_Elements(ra) EQ 0 THEN ra = 5S              ; default square template "radius"
    IF N_Elements(phi0) EQ 0 THEN phi0 = 0.0         ; default cross orientation
    IF N_Elements(nbrushes) EQ 0 THEN nbrushes = 4S  ; default no. brushes
    IF ~keyword_set(tvflag) THEN tvflag = 0B         ; the symbol ~ means logical negation

; define logical variables affecting the display
    tvflag2 = 0S             ; if set=1, then plot raw templates when computed
    cleanup = 0S             ; if set=1 then close plot windows when done
    zoomfactor = 10.0        ; magnify the template images displayed on the screen
    
    waittime = 0.1

;create template x+iy and modify it to  ((x+iy)/(SQRT(x^2 + y^2)))^(-4) for 4-brush patter
    txpixels = 2*ra+1                       ; define size of the (square) template image
    typixels = 2*ra+1
    real_z = fltarr(txpixels, typixels)
    imag_z = fltarr(txpixels, typixels)
    range = findgen(2*ra+1) - ra            ; indexed vector [-ra, ... -1, 0, 1, ... ra]
    FOR i=0, 2*ra DO BEGIN
       real_z[*,i] = range                  ; fill rows with this vector    ( "x" )
       imag_z[i,*] = range                  ; fill columns with this vector ( "y" )
    ENDFOR

    real_z[ra,ra] = 1                       ; force center values to unity
    imag_z[ra,ra] = 1
    z0  = complex(real_z,  imag_z)          ; gives x + iy
    z90 = complex(imag_z,  -real_z)         ; x+iy rotated  by -90 degrees (~ y -ix)   (was +Re, -Im but this rotates the patern clockwise ...)
    z   = z0 * COS(phi0 * !dtor) + z90 * SIN(phi0 * !dtor)           ; combine the orthogonal complex fields at angle phi0


    template = z/SQRT(real_z^2 + imag_z^2)           ; normalize  to unit vectors so have exp[i(phi - phi0)]
;    template = 1./template^4                        ; invert to create a sharp 4-fold cross
    template = 1./temporary(template)^nbrushes                         ; invert to create a sharp cross with nbrushes arms    jem april 2009
    template = (SQRT(real_z^2 + imag_z^2) lt ra) * temporary(template) ; zero the elements outside radius ra

; compute real and imaginary parts of the template
; (although seems none of these was passed back to the main program in the original code - jem)
    realT = (((FLOAT(template) + 1)/2)*255)
    imaginaryT = (((IMAGINARY(template) + 1)/2)*255)
    bothT=realT+imaginaryT

; plot the image templates (scaled up by zoomfactor)
    IF (tvflag AND tvflag2) THEN BEGIN     ; also display raw template images on screen
    ; first the Real part
      IF WindowAvailable(windex0) THEN BEGIN        ; this window already exists
        wset, windex0
      ENDIF ELSE BEGIN                              ; this window does not exist so create it
        window, /FREE,    xpos=400, ypos=500, title='Real Template', $
                              xsize=zoomfactor*txpixels, ysize=zoomfactor*typixels
        windex0 = !D.Window
      ENDELSE
      imsize, realT, x0, y0, xsize, ysize               ; fit existing image to this window
      newimage = congrid (realT, xsize, ysize)          ; resize image to computed xsize,ysize values
      loadct, 0                                         ; load color table
      tvscl, newimage, x0, y0                           ; plot scaled image positioned at x0, y0

    ; then the Imaginary part
      IF WindowAvailable(windex1) THEN BEGIN        ; this window already exists
        wset, windex1
      ENDIF ELSE BEGIN                              ; this window does not exist so create it
        window, /FREE,    xpos=410, ypos=510, title='Imaginary Template', $
                          xsize=zoomfactor*txpixels, ysize=zoomfactor*typixels
        windex1 = !D.Window
      ENDELSE
      imsize, imaginaryT, x0, y0, xsize, ysize          ; fit existing image to this window
      newimage = congrid (imaginaryT, xsize, ysize)     ; resize image to computed xsize,ysize values
      loadct, 0                                         ; load color table
      tvscl, newimage, x0, y0                           ; plot scaled image positioned at x0, y0

    ; Overlay the images:  Real plus Imaginary templates
      IF WindowAvailable(windex2) THEN BEGIN        ; this window already exists
        wset, windex2
      ENDIF ELSE BEGIN                              ; this window does not exist so create it
        window, /FREE,    xpos=420, ypos=520, title='R+I Template', $
                          xsize=zoomfactor*txpixels, ysize=zoomfactor*typixels
        windex2 = !D.Window
      ENDELSE
      imsize, bothT, x0, y0, xsize, ysize               ; fit existing image to this window
      newimage = congrid (bothT, xsize, ysize)          ; resize image to computed xsize,ysize values
      loadct, 0                                         ; load color table
      tvscl, newimage, x0, y0                           ; plot scaled image positioned at x0, y0
    ENDIF

    IF (tvflag) THEN BEGIN     ; display template images on screen
    ; the final Template
      IF WindowAvailable(windex3) THEN BEGIN        ; this window already exists
        wset, windex3
      ENDIF ELSE BEGIN                              ; this window does not exist so create it
        window, /FREE,    xpos=430, ypos=530, title='Final Template', $
                                  xsize=zoomfactor*txpixels, ysize=zoomfactor*typixels
        windex3 = !D.Window
      ENDELSE
      imsize, template, x0, y0, xsize, ysize         ; fit existing image to this window
      newimage = congrid (template, xsize, ysize)    ; resize image to computed xsize,ysize values
      loadct, 0                                        ; load color table
      tvscl, newimage, x0, y0                          ; plot scaled image positioned at x0, y0
    ENDIF

    ; IF (1) THEN WAIT, waittime

    IF (cleanup) THEN BEGIN       ; close the local plot windows
      warray = [windex0, windex1, windex2, windex3]     ; create an array of window indexes
      FOR i=0, N_Elements(warray)-1 DO BEGIN
        IF WindowAvailable(warray(i)) THEN wdelete,warray(i)
      ENDFOR
 
    ENDIF

    RETURN, template
END
