;  porp_move v.1.0

;  Movement model used for simulating harbour porpoise fine scale movement behaviour.
;  Please refer to the scientific publication for detailed documentation: 
;  Nabe-Nielsen, J., Tougaard, J., Teilmann, J., Lucke, K. & Forchhammer, M.C. (2013):
;  "How a simple adaptive foraging strategy can lead to emergent home ranges and increased
;  food intake." Oikos, 122, 1307–1316.


; The model was created as part of the project
; BRIDGES AS BARRIERS PORPOISE MODELLING: DOES THE GREAT BELT BRIDGE HINDER 
; MOVEMENT OF HARBOUR PORPOISES IN THE GREAT BELT
; funded by Fehmern Belt A/S
;

; Copyright (C) 2016, Jacob Nabe-Nielsen <jnn@bios.au.dk>
; 
; This program is free software; you can redistribute it and/or modify it 
; under the terms of the GNU General Public License version 2 and only 
; version 2 as published by the Free Software Foundation.
; 
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
; 
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


; The model was developed and tested using NetLogo version 4.1. Development ended 2011-12-05

 
; debug levels: 
;   0   no debugging
;   1   debug porp_avoid_land / water depth problems
;   2   write turning angles before and after avoiding land
;   3   debugging turning angles -- esp. loop control
;   4   debugging attraction vector caused by reference memory
;   5   debugging deterrence vector caused by working memory
;   6   debugging attraction + deterrence vectors
;   7   debugging porp_get_exp_food_val -- expected value of future food

; behavioural modes in model: 
;   0   Behaviour produced by simple non-autocorrelated Markov process (except when porps seek to avoid land); 
;       variables calibrated using dead-reckoning data.
;   1   Like model 1, but introducing autocorrelation in turning angles and speed, and sharper turning angles 
;       at low speed (i.e. a correlated random walk (CRW) model). 
;   2   Extends model 1 behaviour by introducing a desire to return to return to areas with food (i.e. a reference
;       memory, cf Van Moorter, Oikos 2009). Food is eaten in the entire cell once the porp has been there. 
;       Food doesn't affect step length (intrinsic behaviour).


;   7   Extends behaviour model 1. Optimal foraging model where animal behaviour is influenced by the probability 
;       of encountering food in each of the patches it passes. The probability of finding food is assumed to decrease
;       linearly with water depth (to make porps stay on shallow water as observed in the satellite track data).
;       Food quality (energy obtained from the food found) is assumed to be constant. Porpoises are able to learn
;       where they are most likely to find food (i.e. at which depth), but do not remember where they have been.
;   8   Builds on model 3. Enables variations in food quality. 
;       This will allow porpoises to develop different feeding strategies depending on their experience: They may
;       go for low-quality food which is encountered with a high probability at shallow water, or they may go for 
;       high-quality food (schoals of cod or herring) on deeper waters.
;   9   Like behaviour model 4, but allowing porpoises to remember where they have been (this is probably necessary
;       in order to get them to stay in the same region for a long time, and maybe it is sufficient to make them
;       shift region from time to time).



extensions [ gis ]
globals [
  my_tick
  days
  time                   ; hours since simulation start
  path                   ; path to directory for input/output, one dir for each area
  outfile                ; simulation output
  corr_logmov            ; correlation in movement distance in CRW +
  corr_angle             ; correlation in direction in CRW +
  bathy_data 
  vt                     ; resultant attraction vector, resulting from reference memory of food availability (model >=2)
  total_food_eaten       ; Used for testing whether optimal foraging strategy is used
  salinity_data
  sediment_data
  bot_ns_data
  bot_ew_data
  food_prob01_data
  xllcorner 
  yllcorner
  deploy_x               ; x-coord of porp 0
  deploy_y               ; y-coord of porp 0
  turn_right             ; random variable; turn right if = 1, left if = -1
  min_depth
  movlgt_list            ; list of 100 move lengths for porpoise 0
  angle_list             ; list of 100 turning angles for porpoise 0
  max_movlgt             ; monotonously increasing, move length/100m
  mean_movlgt
  use_exp_food_val       ; get more attracted to the CRW path if food was found recently
  CRW_contrib            ; length of vector pointing in direction predicted by CRW
  MR_contrib             ; length of vector pointing in direction of remembered food
 ]
 
breed [
  porps porp 
]

patches_own [ 
  bathymetry
  salinity
  sediment
  bot-ns-current
  bot-ew-current
  food-prob01             ; probability of finding food in a patch
  food_level              ; amount of food in a patch
]

porps_own [ 
  energy-level           ; Porpoises get energy by eating and loose energy by moving.
  porp_length
  prev_angle             ; Last turning angle (not the heading!)
  pres_angle             ; Present turning angle
  prev_logmov            ; Previous Log10 (move length /100 m)
  pres_logmov            ; Present Log10 (move length /100 m)
  enough-water-ahead     ; Turn to avoid land if false
  pos_list               ; Coordinates of previous positions -- latest positions in left end
  ptt                    ; Argos id-number -- the number of the simulated porpoise
  sex
  p-length
  p-weight
  ; Vars added in model 2
  ref-mem-strength-list  ; Memory decay with time (logistic decrease, function of rR)
  work-mem-strength-list ; Memory of poor quality patches, decays with time (logistic decrease, function of rW)
  work-mem-updated       ; Whether a list of working memory has been updated
  stored-util-list       ; Up_t ; Remembered feeding success (after memory decay) -- latest positions left
  VE-total               ; total value of food expected to be found in the future
]


; Landscape variables:

to landsc_setup
  ;; (for this model to work with NetLogo's new plotting features,
  ;; __clear-all-and-reset-ticks should be replaced with clear-all at
  ;; the beginning of your setup procedure and reset-ticks at the end
  ;; of the procedure.)
  __clear-all-and-reset-ticks
  reset-timer
  ; Note that setting the coordinate system here is optional, as
  ; long as all of your datasets use the same coordinate system.
  ; gis:load-coordinate-system (word "data/" projection ".prj")
  ; Load datasets:
  set path word "raster-data/" area      ; note that "/" may be a problem on Windows
  set path word path "/"
  set bathy_data gis:load-dataset word path "bathy.asc"
  set salinity_data gis:load-dataset word path "salinity.asc"
  set sediment_data gis:load-dataset word path "sediment.asc"
  set bot_ns_data gis:load-dataset word path "bot_ns.asc"
  set bot_ew_data gis:load-dataset word path "bot_ew.asc"
  set food_prob01_data gis:load-dataset word path "food-prob01.asc"

  ; set variables:
  ; adjust lower left corner:
  if (area = "greatbelt") [
    set xllcorner 550673  ; header for "greatbelt.asc", position in m
    set yllcorner  6100242
  ]
  if (area = "anholt") [
    set xllcorner 575473
    set yllcorner 6220242
  ]
  if (area = "homogeneous") [
    set xllcorner 575473
    set yllcorner 6220242
  ]
  if (area = "fehmarn") [
    set xllcorner 603473  ; header for "greatbelt.asc", position in m
    set yllcorner 5988242
  ]

  ; Set the world envelope to the union of all of our dataset's envelopes
  gis:set-world-envelope (gis:envelope-union-of 
    (gis:envelope-of bathy_data)
    (gis:envelope-of salinity_data)
    (gis:envelope-of sediment_data)
    (gis:envelope-of bot_ns_data)
    (gis:envelope-of bot_ew_data)
    (gis:envelope-of food_prob01_data)
    )
  ; This is the preferred way of copying values from a raster dataset
  ; into a patch variable: in one step, using gis:apply-raster.
  gis:apply-raster bathy_data bathymetry
  gis:apply-raster salinity_data salinity
  gis:apply-raster sediment_data sediment
  gis:apply-raster bot_ns_data bot-ns-current
  gis:apply-raster bot_ew_data bot-ew-current
  gis:apply-raster food_prob01_data food-prob01
  
  ; set amount of food -- if there is a chance that food is present
  ask patches [ ifelse food-prob01 > 0 [ set food_level maxU ] [ set food_level food-prob01 ] ]

  landsc_display
  let tmp-t word "Setup time: " timer
  print word tmp-t " sec"
end


to landsc_display  ; Updates the display variable
  no-display
  if (disp-var = "bathymetry") [ 
    let min-bathymetry gis:minimum-of bathy_data
    let max-bathymetry gis:maximum-of bathy_data
    ask patches
    [ ; note the use of the "<= 0 or >= 0" technique to filter out 
      ; "not a number" values, as discussed in the documentation.
      set pcolor 39
      if (bathymetry <= 0) or (bathymetry >= 0)
      [ set pcolor scale-color blue bathymetry max-bathymetry min-bathymetry ] 
    ]
  ]
  if (disp-var = "salinity") [ 
    let min-salinity gis:minimum-of salinity_data
    let max-salinity gis:maximum-of salinity_data
    ask patches
    [ 
      set pcolor 39
      if (salinity <= 0) or (salinity >= 0)
      [ set pcolor scale-color yellow salinity max-salinity min-salinity ] 
    ]
  ]
  if (disp-var = "sediment") [ 
    let min-sediment gis:minimum-of sediment_data
    let max-sediment gis:maximum-of sediment_data
    ask patches
    [ 
      set pcolor 39
      if (sediment <= 0) or (sediment >= 0)
      [ set pcolor scale-color green sediment max-sediment min-sediment ] 
    ]
  ]
  if (disp-var = "NS current") [ 
    let min-bot-ns gis:minimum-of bot_ns_data
    let max-bot-ns gis:maximum-of bot_ns_data
    ask patches
    [ 
      set pcolor 39
      if (bot-ns-current <= 0) or (bot-ns-current >= 0)
      [ set pcolor scale-color grey bot-ns-current max-bot-ns min-bot-ns ] 
    ]
  ]
  if (disp-var = "EW current") [ 
    let min-bot-ew gis:minimum-of bot_ew_data
    let max-bot-ew gis:maximum-of bot_ew_data
    ask patches
    [ 
      set pcolor 39
      if (bot-ew-current <= 0) or (bot-ew-current >= 0)
      [ set pcolor scale-color grey bot-ew-current max-bot-ew min-bot-ew ] 
    ]
  ]
  if (disp-var = "food_level") [ 
    let min-food_level 0 ; gis:minimum-of food_prob01_data
    let max-food_level maxU ; gis:maximum-of food_prob01_data
    ask patches
    [ 
      set pcolor 39
      if (food_level <= 0) or (food_level >= 0) [
      ; [ set pcolor scale-color green food_level (1.0 * max-food_level) (1.0 * min-food_level) ] 
        if ( food_level = 0 ) [ set pcolor white ]
        if ( food_level > 0 and food_level <= 0.1 * maxU ) [ set pcolor 48 ]
        if ( food_level > 0.1 * maxU and food_level <= 0.25 * maxU ) [ set pcolor 67 ]
        if ( food_level > 0.25 * maxU and food_level <= 0.5 * maxU ) [ set pcolor 65 ]
        if ( food_level > 0.5 * maxU and food_level < 1 * maxU ) [ set pcolor 63 ]
        if ( food_level = maxU ) [ set pcolor 61 ]
      ]
    ]
  ]
  display
end


; Porpoise variables

to porps_setup_ref
  ; reference porpoises (with who = 0) -- deployment information
  if (area = "greatbelt") [
    set deploy_x ( 606399.6 - xllcorner ) / 100  ; start-position, in pixels
    set deploy_y ( 6148048 - yllcorner ) / 100
    set ptt "2000-04542"
    set p-length 116
    set sex "F"
  ]
  if (area = "fehmarn") [
    set deploy_x ( 636899.6 - xllcorner ) / 100  ; start-position, in pixels
    set deploy_y ( 6108048 - yllcorner ) / 100
    set ptt "F01"
    set sex "NA"
  ]
  if (area = "anholt") [
    set deploy_x ( 619399.6 - xllcorner ) / 100  ; start-position, in pixels
    set deploy_y ( 6148048 - yllcorner ) / 100
    set ptt "A01"
    set sex "NA"
  ]
  if (area = "homogeneous") [
    set deploy_x ( 619399.6 - xllcorner ) / 100  ; start-position, in pixels
    set deploy_y ( 6148048 - yllcorner ) / 100
    set ptt "H01"
    set sex "NA"
  ]
  setxy deploy_x deploy_y
  set color red
  let tmp word "pttid " ptt
  print word tmp " deployed (red dot)"
end


to porp_inspect_0
  inspect porp 0
  print word "porp 0, true x-cor: " (( [ xcor ] of porp 0 ) * 100 + xllcorner ) 
  print word "porp 0, true y-cor: " (( [ ycor ] of porp 0 ) * 100 + yllcorner )
end


to porp_ref_mem_turn
  ; Move towards places visited previously if food was found there and they aren't too far away or forgotten.
  let bb ( [food_level] of patch-here )  ; benthic food species. The stored intrisic patch utility for t=0. Initially it is either 0, 1, or -9999, but grows logistically after food is eaten
  set total_food_eaten ( bb + total_food_eaten )
  if not ( abs(bb) >= 0 ) [  
    ; There are errors in food availability -- sometimes Na is calculated even though depth is > 0. Catch error here
    set bb 0
    if ( debug = 4 ) [ 
      print "Replaced NaN food value with 0"
      print patch-here
    ]
  ]
  set stored-util-list fput bb stored-util-list
  
  ; Update reference memory strength for past locations
  let max-mem 0.999
  set ref-mem-strength-list fput max-mem ref-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009
  let ii 1
  while [ ii < length ref-mem-strength-list ] [  
    let MRPt item (ii - 1) ref-mem-strength-list  ; (reference memory for patch p at time t)
    let reduced-mem MRPt - ref-mem-decay * (1 - MRPt) * MRPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter
    set ref-mem-strength-list replace-item ii ref-mem-strength-list reduced-mem
    set ii ii + 1
  ]

  ; Set patch value for each past location -- perceived patch utility (= reference memory x intrinsic patch utility (stuff eaten)), divided by DIRECT distance
  let perceived-util-list [ ]
  set perceived-util-list lput (item 0 stored-util-list * max-mem) perceived-util-list
  let tmp list (0) (0)
  let attr-vector-list [ ] ; each element in the list is a list with an x-and a y-direction. Vector for first element (this place) has length 0, the others have length 1
  set attr-vector-list lput tmp attr-vector-list
  let one-attr-vector [ ]
  let vector-lgt 0
  set ii 1
  let dist-to-foodpos 0
  while [ ii < length pos_list ] [
    set dist-to-foodpos (distancexy ( item 0 (item ii pos_list) ) ( item 1 (item ii pos_list) ))
    ifelse (dist-to-foodpos < 1E-20 )
      [ set perceived-util-list lput 9999 perceived-util-list ]      ; arbitrary large value for close dist
      [ set perceived-util-list lput ( (item ii stored-util-list) * (item ii ref-mem-strength-list) / dist-to-foodpos ) perceived-util-list ]
    ; = utility * memory / distance
    ; Create attraction vectors; unit-vectors pointing towards the patches in memory
    set one-attr-vector list ((item 0 (item ii pos_list)) - xcor)  ((item 1 (item ii pos_list)) - ycor)
    ; make sure that it works with wrapping landscapes:
    if ( item 0 one-attr-vector > world-width / 2 ) [ set one-attr-vector replace-item 0 one-attr-vector ( item 0 one-attr-vector - world-width ) ]
    if ( item 0 one-attr-vector < (- world-width / 2 ) ) [ set one-attr-vector replace-item 0 one-attr-vector ( item 0 one-attr-vector + world-width ) ]
    if ( item 1 one-attr-vector > world-height / 2 ) [ set one-attr-vector replace-item 1 one-attr-vector ( item 1 one-attr-vector - world-height ) ]
    if ( item 1 one-attr-vector < (- world-height / 2 ) ) [ set one-attr-vector replace-item 1 one-attr-vector ( item 1 one-attr-vector + world-height ) ]
    set vector-lgt sqrt ( item 0 one-attr-vector * item 0 one-attr-vector + item 1 one-attr-vector * item 1 one-attr-vector )
    if vector-lgt = 0 [
      if ( debug = 4 ) [ 
        show word "attr-vector-lgt = " vector-lgt
        print "skipping to next porp"
      ]
      stop
    ]
    set one-attr-vector replace-item 0 one-attr-vector ((item 0 one-attr-vector) / vector-lgt)
    set one-attr-vector replace-item 1 one-attr-vector ((item 1 one-attr-vector) / vector-lgt)
    set attr-vector-list lput one-attr-vector attr-vector-list
    set ii ii + 1
  ]
  ; Use following line to rescale value to sum to 1.
  ; I DON'T DO THIS !!! WHY SHOULD PORPOISES BE ATTRACTED TO OLD POINTS FAR AWAY AT ALL? THEN MEM DECAY DOESN'T MAKE SENSE
  ; Unlike Van Moorter I DON'T account for directional persistance, which is part of the default behaviour for the animal.
  ; The resultant attraction vector (v_t) is the sum of the attractions to all patches in memory.
  ; let val-sum sum perceived-util-list
  ; set ii 0
  ; if val-sum != 0 [
  ;   while [ ii < length perceived-util-list ] [
  ;     set perceived-util-list replace-item ii perceived-util-list (( item ii perceived-util-list ) / val-sum )
  ;     set ii ii + 1
  ;   ]
  ; ]
  ; Calculate resultant attraction vector vt as sum of products of individual values and attraction vectors (eqn 5). May have length != 1
  set ii 1  ; no attraction to current pos (t=0)
  let vt-x 0
  let vt-y 0
  while [ ii < length pos_list ] [
    set vt-x vt-x + item ii perceived-util-list * item 0 ( item ii attr-vector-list )
    set vt-y vt-y + item ii perceived-util-list * item 1 ( item ii attr-vector-list )
    set ii ii + 1
  ]
  if ( debug = 4 ) [ 
    type word "Food here: " bb
    type ",  Attr.v: "
    let attr-vect list vt-x (",")
    set attr-vect lput vt-y attr-vect
    print attr-vect
    if (not ( abs(vt-x) >= 0)) [  ; catch missing values
      write "Perc.util: "
      print perceived-util-list
    ]
  ]
  set vt list vt-x vt-y
  
  ; Remove items in distant past to increase execution speed
  if ( length ref-mem-strength-list > memory-max ) [ set ref-mem-strength-list remove-item memory-max ref-mem-strength-list ]
  if ( length stored-util-list > memory-max ) [ set stored-util-list remove-item memory-max stored-util-list ]
end  ; end porp_ref_mem_turn


to porp_work_mem_turn
   ; Influences direction moved in std-move through vector 'vt'
   ; This procedure MUST be called after porp_ref_mem_turn, as it uses the stored-util-list (experienced food at patch) calculated there,
   ; and because it updates the vt vector which is later used in std.move (adds to the previously calculated vt).

  ; Update working memory strength (short-term memory) for past locations
  let max-mem 0.999
  let MWPt 0
  set work-mem-strength-list fput max-mem work-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009

  let ii 1
  if ( work-mem-updated = false ) [
    while [ ii < length work-mem-strength-list ] [  
      set MWPt item (ii - 1) work-mem-strength-list  ; (working memory for patch p at time t)
      let reduced-mem MWPt - work-mem-decay * (1 - MWPt) * MWPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter 
      set work-mem-strength-list replace-item ii work-mem-strength-list reduced-mem
      set ii ii + 1
    ]
  ]
  set work-mem-updated true
  ; Presently no need to multiply with stored-util-list to get perceived utility -- the animal knows that all food is eaten in the patch.

  let tmp list (0) (0)
  let deter-vector-list [ ] ; each element in the list is a list with an x-and a y-direction. Vector for first element (0; this place) has length 0, the others have length 1
  set deter-vector-list lput tmp deter-vector-list 
  ;  show word "deter-vector-list: " deter-vector-list
  let one-deter-vector [ ]
  let vector-lgt 0
  set ii 1
  ;  let dist-to-foodpos 0
  while [ ii < length pos_list ] [
    ; Create deterrence vectors; unit-vectors pointing towards the patches in memory
    set one-deter-vector list ((item 0 (item ii pos_list)) - xcor)  ((item 1 (item ii pos_list)) - ycor)
    ; make sure that it works with wrapping landscapes:
    if ( item 0 one-deter-vector > world-width / 2 ) [ set one-deter-vector replace-item 0 one-deter-vector ( item 0 one-deter-vector - world-width ) ]
    if ( item 0 one-deter-vector < (- world-width / 2 ) ) [ set one-deter-vector replace-item 0 one-deter-vector ( item 0 one-deter-vector + world-width ) ]
    if ( item 1 one-deter-vector > world-height / 2 ) [ set one-deter-vector replace-item 1 one-deter-vector ( item 1 one-deter-vector - world-height ) ]
    if ( item 1 one-deter-vector < (- world-height / 2 ) ) [ set one-deter-vector replace-item 1 one-deter-vector ( item 1 one-deter-vector + world-height ) ]
    set vector-lgt sqrt ( item 0 one-deter-vector * item 0 one-deter-vector + item 1 one-deter-vector * item 1 one-deter-vector )
    if vector-lgt = 0 [
      if ( debug = 5 ) [ 
        show word "deter-vector-lgt = " vector-lgt
        print "Haven't moved, skipping to next porp"
      ]
      stop
    ]
    set one-deter-vector replace-item 0 one-deter-vector ((item 0 one-deter-vector) / vector-lgt)
    set one-deter-vector replace-item 1 one-deter-vector ((item 1 one-deter-vector) / vector-lgt)
    set deter-vector-list lput one-deter-vector deter-vector-list
    set ii ii + 1
  ]

  ; Calculate resultant deterrence vector vtd as sum of products of individual values and deterrence vectors
  set ii 1  ; no deterrence from current pos (t=0)
  let vtd-x 0
  let vtd-y 0
  while [ ii < length pos_list ] [
    set vtd-x vtd-x + inertia-const * item ii work-mem-strength-list * item 0 ( item ii deter-vector-list )
    set vtd-y vtd-y + inertia-const * item ii work-mem-strength-list * item 1 ( item ii deter-vector-list )
    set ii ii + 1
  ]

  if ( debug = 5 ) [ 
    print word "work-mem: " work-mem-strength-list
    type "Deter.v: "
    let deter-vect list vtd-x (",")
    set deter-vect lput vtd-y deter-vect
    print deter-vect
    if (length pos_list > 1) [ print word "pos. before: " item 1 pos_list ]
    print word "pos. now: " item 0 pos_list
    print ""
    ; Checked -- works, at least with length pos_list = 2
  ]
  if ( debug = 6 ) [ 
    type "Attr.v.1: "
    let attr-vect list ( precision (item 0 vt) 3 ) (",")
    set attr-vect lput ( precision (item 1 vt) 3 )  attr-vect
    print attr-vect
    ;print "       (before deterr)"
  ]
  
  ; vtd points towards the previous position, must be subtracted from vt
  set vt replace-item 0 vt ( item 0 vt - vtd-x )
  set vt replace-item 1 vt ( item 1 vt - vtd-y )

  if ( debug = 6 ) [ 
    type "Attr.v.2: "
    let attr-vect list ( precision (item 0 vt) 3 ) (",")
    set attr-vect lput ( precision (item 1 vt) 3 )  attr-vect
    print attr-vect
    ; print "      (incl deterr)"
    print ""
  ]

  ; Remove items in distant past to increase execution speed
  if ( length work-mem-strength-list > memory-max ) [ set work-mem-strength-list remove-item memory-max work-mem-strength-list ]
end ; end porp_work_mem_turn


to porp_get_exp_food_val
  ; Calculate the expaected value (VE-total) of the food to be found in the future based on food found in recent positions x the working memory
  ; Uses the values of the patches in "stored-util-list", calculated in porp_ref_mem_turn

  ; Update working memory strength (short-term memory) for past locations
  let max-mem 0.999
  let MWPt 0
  set work-mem-strength-list fput max-mem work-mem-strength-list   ; Nearly perfect memory of what just happened, cf Van Moorter et al 2009
  let ii 1
  if ( work-mem-updated = false ) [  ; list may have been updated in porp_work_mem_turn
    while [ ii < length work-mem-strength-list ] [ 
      set MWPt item (ii - 1) work-mem-strength-list  ; (working memory for patch p at time t)
      let reduced-mem MWPt - work-mem-decay * (1 - MWPt) * MWPt ; logistic decrease in memory with time, see eqn 2 in Van Moorter 
      set work-mem-strength-list replace-item ii work-mem-strength-list reduced-mem
      set ii ii + 1
    ]
  ]
  set work-mem-updated true

  set ii 1
  set VE-total 0
  let max-i min list ( length work-mem-strength-list ) ( memory-max )
  while [ ii < max-i ] [ 
    set VE-total VE-total + item (ii - 1) work-mem-strength-list * item (ii - 1) stored-util-list
    set ii ii + 1
  ]

  if ( debug = 7 ) [ 
    print ""
    ; print word "stored-util-list: " stored-util-list
    ; print word "work-mem-strength-list: " work-mem-strength-list
    print word "VE-total: " VE-total
  ]
end ; end porp_get_exp_food_val


to porp_std_move
  ; Movements based on dead-reckoning data:
  ; global vars: corr_logmov 0.94 and corr_angle 0.26
  ; ### turning angle
  let prev-mov 10 ^ prev_logmov
  let pres_heading heading
  set pres_angle 999
  let j 1
  let tmp-angle 0
  ifelse ( prev_angle < 0 ) [ set tmp-angle prev_angle - 24 ] [ set tmp-angle prev_angle + 24 ]  ; for increasing mean turning angle
  while [ abs (pres_angle) > 180 ]  [ 
    set pres_angle ( tmp-angle * (- corr_angle) + random-normal 0 38 )                  ; Autoreg can't be used for estimating param. as estimated turns are changed if on shallow water. 
    set j j + 1
    if (j = 200) [
      set pres_angle ( pres_angle * 90 / (abs pres_angle))
      if ( debug = 3 ) [ 
        print word "exiting loop 1, ang=" pres_angle
      ]
    ]
  ]
  let sign 0
  ifelse pres_angle < 0 [ set sign -1 ] [ set sign 1 ]
  set pres_angle abs pres_angle ; add the sign again later
  ; Make angle decrease linearly with mov-dist
  let go-on true
  set j 1
  let rnd 0
  while [ go-on ]  [ 
    set rnd random-normal 96 28      ; draws the number to be added to pres_angle
    if ( prev-mov <= 5.50 ) [ 
      set pres_angle pres_angle + rnd - ( rnd * prev-mov / 5.50 )
    ]
    if ( pres_angle < 180 ) [ set go-on false ]  ; remember that turning angle is unsigned here
    set j j + 1
    if (j = 200) [
      set pres_angle ( random 20 + 90 )
      set go-on false
      if ( debug = 3 ) [ 
        print word "exiting loop 2, ang=" pres_angle
      ]
    ]
  ]
  ; if ( abs pres_angle > 55 and abs pres_angle < 180 ) [ set pres_angle (1 - random-float 0.32) * pres_angle ]  ; make turning angle dist more leptokurtic
  set pres_angle pres_angle * sign
  let angle-before-avoid-land pres_angle ; for printing later using debug 2
  right pres_angle
  let angle-turned-right pres_angle ; for updating prev_angle at end of porp_std_move
  set pres_angle 0

  ; ### distance
  set pres_logmov 999
  let porp-max-dist 1.18                                                                                        ; log10 ( max distance a porpoise can travel per half-hour )
  set j 1
  while [ pres_logmov > porp-max-dist ] [ 
    set pres_logmov ( corr_logmov * prev_logmov + random-normal 0.42 0.48 ) 
    set j j + 1
    if (j = 200) [
      if (pres_angle = 0) [set pres_angle pres_angle + 0.00001]
      set pres_angle ( pres_angle * 90 / (abs pres_angle))
      if ( debug = 3 ) [ 
        print word "exiting loop 3, ang=" pres_angle
      ]

    ]    
  ]    ; Mean pres_mov should be x.x x100 m 
  let pres_mov ( 10 ^ pres_logmov )        ; This is what is plotted in the histogram
  ;
  ; Turn to avoid swimming on land if necessary:
  set enough-water-ahead false
  let count-i 0
  while [ not enough-water-ahead ] [
    porp_check_depth
    if (not enough-water-ahead) [ porp_avoid_land ]
    set pres_mov ( 10 ^ pres_logmov )                                                      ; because pres_logmov may have changed in the porp_avoid_land procedure
    right pres_angle                                                                       ; angle to turn -- pres_angle -- is changed in porp_avoid_land
    set angle-turned-right (angle-turned-right + pres_angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    set pres_angle 0
    set count-i count-i + 1
    if (count-i = 100) [ 
      set enough-water-ahead true 
      if ( debug = 1 ) [ 
        print "caught water-ahead loop"
      ]
    ]
  ]
  ; test depth again, avoid_beh = 5:
  porp_check_depth
  if (not enough-water-ahead) [
    let prev-heading heading
    let p max-one-of neighbors [ bathymetry ]
    face p
    set angle-turned-right (angle-turned-right + pres_angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    set pres_mov 1                                                                                            ; move 100 m towards deeper patch
    if ( debug = 1 ) [ 
      let tmp5 list "beh =  5 ; tck " my_tick
      write tmp5 
    ]
  ]
  ;
  ; slow down if turning sharply:
  if ( pres_mov > 10 and (abs angle-turned-right) > 90 ) [ set pres_mov  pres_mov / 5  ] 
  if ( pres_mov > 7 and (abs angle-turned-right) > 50 ) [ set pres_mov  pres_mov / 2  ] 
  ;
  ; Change direction if attracted / repelled by certain areas (model >= 2)
  let total-dx 0
  let total-dy 0
  if ( not use_exp_food_val ) [
    set total-dx (dx * pres_mov) + (item 0 vt)      ; vt isn't used in porp_std_move till here
    set total-dy (dy * pres_mov) + (item 1 vt)      ; note that dx is change in x if taking ONE step forward
    facexy (xcor + total-dx) (ycor + total-dy)
  ]
  if ( use_exp_food_val ) [
    set CRW_contrib inertia-const + pres_mov * VE-total       ; length of vector pointing in direction predicted by CRW
    set MR_contrib sqrt ( (item 0 vt) * (item 0 vt) + (item 1 vt) * (item 1 vt) )     ; length of vector pointing in direction of remembered food
    set total-dx (dx * CRW_contrib) + (item 0 vt)
    set total-dy (dy * CRW_contrib) + (item 1 vt)
    facexy (xcor + total-dx) (ycor + total-dy)                ; really not needed, it already points that way
  ]
  ; Store turn for calc of turning angle in next step:
  ; let total-turn heading - pres_heading   ; total change in heading, including all adjustments till here. 'pres_heading' was calc in beginning of porp_std_move
  let total-turn subtract-headings heading pres_heading   ; total change in heading, including all adjustments till here. 'pres_heading' was calc in beginning of porp_std_move
  ;
  ; Move: 
  fd pres_mov  ; movement length isn't affected by presence of food
  ;
  if ( debug = 2 ) [ 
    if ( my_tick = 0 ) [ print "dist angle-before-avoid-land angle-turned-right x y" ]
    let tmp-var2 (round ( (10 ^ prev_logmov) * 100) / 100)    ; THIS IS IMPORTANT -- the porp turns before it moves, so turning angle is affected by previous moving dist
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 angle-before-avoid-land
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 angle-turned-right
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 ( xcor * 100 + xllcorner )
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 ( ycor * 100 + yllcorner )
    print tmp-var2
  ]
  if ( debug = 7 ) [ 
    print word "CRW_contrib: " list ( (dx * (inertia-const + pres_mov * VE-total)) ) ( (dy * (inertia-const + pres_mov * VE-total)) )
    print word "MR_contrib: " vt
    print word "dx, dy (after): " list (total-dx) (total-dy)
    print word "heading (after): " heading
    print word "total-turn: " heading
  ]   
  ;
  ; Remember current moves for the next iteration
  ; if attraction to food alters the movement angle (i.e. vt != 0), this isn't remembered for next step
  ; set prev_angle angle-turned-right  ; so the additional turn due to attraction to food does not influence turning angle in next step 
  set prev_angle total-turn  ; so the additional turn due to attraction to food DOES influence turning angle in next step 
  set prev_logmov log pres_mov 10  ; total steplength, resulting from vt + pres_mov
  ;
  ; test depth one last time, avoid_beh = 6 - move back on same track:
  if (not ([ bathymetry ] of patch-here > 0) ) [
    let prev-heading heading
    facexy (item 0 item 1 pos_list) (item 1 item 1 pos_list)
    set angle-turned-right (angle-turned-right + pres_angle)
    if (angle-turned-right > 180) [ set angle-turned-right angle-turned-right - 360 ]
    if (angle-turned-right < -180) [ set angle-turned-right angle-turned-right + 360 ]
    setxy (item 0 item 1 pos_list) (item 1 item 1 pos_list)                              ; move 100 m towards deeper patch
    if ( debug = 1 ) [ 
      ; print "; -> "
      let tmp6 list "beh =  6 � ; tck " my_tick
      set tmp6 word tmp6 " ; "
      set tmp6 word tmp6 angle-turned-right
      print word tmp6 " degr."
    ]
  ]
  ; update position list:
  let pres_pos list xcor ycor
  set pos_list fput pres_pos pos_list
  if ( length pos_list > memory-max ) [ set pos_list remove-item memory-max pos_list ]   
end ;  end porp_std_move


to porp_eat_food
  ; reduce food to 0 in the patch that the porp just left. Increase food level in cells with food-prob01 > 0 AFTERWARDS in order to calc. stored-util-list correctly.
  ; Eat all food in the patch in one go and change patch colour to reflect food availability.
  ask patch (item 0 (item 1 pos_list)) (item 1 (item 1 pos_list)) [ ; item 0 pos_list was the last added element, i.e. the current position
    ifelse food_level > 0 
    [ 
      set food_level 0.01
      set pcolor yellow
    ]
    [ ]
  ]

end ; end porp_eat_food



; Statistics and plots
 
to my_update_plots ; update histograms and other plots 
  ; move length
  set movlgt_list fput ( 10 ^ ( [ pres_logmov ] of porp 0 ) ) movlgt_list
  if ( length movlgt_list > 100 ) [ set movlgt_list remove-item 100 movlgt_list ]   ; base histogram on list of fixed length
  set max_movlgt max fput max_movlgt movlgt_list
  set mean_movlgt mean movlgt_list
  set-current-plot "movlgt-porp-0"
  histogram movlgt_list
  ; turning angle
  set angle_list fput ( [ prev_angle ] of porp 0 ) angle_list
  if ( length angle_list > 100 ) [ set angle_list remove-item 100 angle_list ]   ; base histogram on list of fixed length
  set-current-plot "angle-porp-0"
  histogram angle_list
  ; abs turning angle vs distance moved, porp 0 -- SOMETHING WRONG WITH PLOT
  ;  set-current-plot "angle-vs-dist"
  ;  plotxy (item 1 movlgt_list) ( abs (item 0 angle_list) )    ; the latest angle is a function of the movelength in prev. step -- hence steplgt and angle are extracted from diff places in lists
  ; Plot memory related moves, models �2
  set-current-plot "memory-move-porp-0"
  set-current-plot-pen "CRW_contrib"
  plot CRW_contrib
  set-current-plot-pen "MR_contrib-x100"
  plot MR_contrib * 100
end


; File I/O

to file_setup
  let repl-counter 0                            ; different replicates of the same porpoise
  let go-on true
  while [go-on] [
    set repl-counter repl-counter + 1
    let repl-counter-2 word "1000" repl-counter
    set repl-counter-2 substring repl-counter-2 (length repl-counter-2 - 3) (length repl-counter-2)
    let file-tmp word path "output/"
    set file-tmp word file-tmp "sim-"
    set file-tmp word file-tmp  [ ptt ] of porp 0
    set file-tmp word file-tmp "-mod"
    set file-tmp word file-tmp model
    set file-tmp word file-tmp "-"
    set file-tmp word file-tmp repl-counter-2
    set outfile word file-tmp ".txt"
    if (not file-exists? outfile) [ 
      set go-on false
      file-open outfile
    ]
    if (repl-counter = 300) [ 
      set go-on false 
      user-message ( "The desired number of replicates has been obtained" )
    ]
  ]
  file-print ("animal ptt pop sex length weight x y bathy time rW rR rU Umax") ; header line (space-separated). Note that "bathy" is non-standard
end


to file_write_line
  ; Ask porpoise to write variables that should be imported into an R track object:
  ; "animal ptt pop sex length weight x y year month day hour minute second":
  file-open outfile ; append
  ; add entries
  file-write "NA" ; not followed by CR
  file-write ptt
  file-write "IDW"
  file-write sex 
  file-write p-length
  file-write p-weight
  file-write (( [ xcor ] of porp 0 ) * 100 + xllcorner )
  file-write (( [ ycor ] of porp 0 ) * 100 + yllcorner )
  file-write [ bathymetry ] of patch-here
  file-write time
  file-write work-mem-decay ; rW
  file-write ref-mem-decay ; rR
  file-write food_growth_rate ; rU
  file-write maxU
  file-print ""    ; followed by CR
  file-close
end


to file_loop  ; corresponds to the wrt button -- makes all output files
  set write_data true
  repeat 10 [
    porps_setup
    repeat 15000 [  ; max my_tick, i.e. nearly a year
      go
    ]
  ]
  set write_data false
end


to porps_setup
  clear-turtles
  clear-output
  clear-drawing   ; clears pendown tracks
  clear-all-plots
  reset-ticks
  set corr_logmov 0.94
  set corr_angle 0.26
  set vt list 0 0
  set my_tick 0
  set days 0
  set time 0
  set max_movlgt 0
  set mean_movlgt 0
  set use_exp_food_val false
  set CRW_contrib -9
  set MR_contrib -9 
  set min_depth 1               ; water depth where porps very rarely swim (according to Jakob)
  set turn_right 1
  set movlgt_list list (0) (0)  ; two inputs required...  list used for histogramming
  set movlgt_list remove-item 0 movlgt_list
  set movlgt_list remove-item 0 movlgt_list
  set angle_list list (0) (0)  ; two inputs required...  list used for histogramming
  set angle_list remove-item 0 angle_list
  set angle_list remove-item 0 angle_list
  create-porps n-porps
  ask porps [ 
    set prev_logmov 0.8 ; unit = patchsize, i.e. times 100 m
    set prev_angle 10
    set deploy_x random-xcor
    set deploy_y random-ycor
    set enough-water-ahead true
    setxy deploy_x deploy_y
    let pos list deploy_x deploy_y
    set pos_list list (0) (pos)
    set pos_list remove-item 0 pos_list
    set ref-mem-strength-list [ ]
    set work-mem-strength-list [ ]
    set work-mem-updated false
    set VE-total 0
    set stored-util-list  [ ]
    if ( not ( [bathymetry] of patch-here > 1 ) ) [ die ]
    set color orange
    set shape "circle"
    set size 5
    if (who = 0) [ porps_setup_ref ]
    ifelse track [ pen-down ] [ pen-up ]
    set pen-size 0.1
  ]
  carefully [ 
    set movlgt_list lput ( 10 ^ ( [ prev_logmov ] of porp 0 ) ) movlgt_list
  ]
  [
    set movlgt_list lput 1 movlgt_list
  ]
  if ( not (is-turtle? ( porp 0 )) ) [ porps_setup ]    ; strange that I have to do this... the porpoise isn't always created in the first go
  
  ; update food level -- also done in landsc_setup, but convenient to repeat it here
  ask patches [ 
    ifelse food-prob01 > 0 [ set food_level maxU ] [ set food_level food-prob01 ] 
  ]
  landsc_display ; update displayed amount of food etc

end

to porp_move
  if ( model = 0 ) [ porp_markov_move ]
  if ( model = 1 ) [ porp_std_move ]
  if ( model = 2 ) [ 
    set work-mem-updated false
    set use_exp_food_val true
    porp_ref_mem_turn  ; get attracted to places where food was found. Influences direction moved in std-move through vector 'vt'

    set vt replace-item 0 vt ( item 0 vt * B )  ; Weight of reference memory part
    set vt replace-item 1 vt ( item 1 vt * B )
  
    porp_get_exp_food_val
    porp_std_move
    porp_eat_food      ; food level increases in 'go'
  ]
  if not ( [ bathymetry ] of patch-here > 0 ) [ 
    follow-me
    beep
    user-message "Error, no water"
  ]
end


to porp_markov_move
  ; Movements based on dead-reckoning data -- first calc distance, then turning angle
  set pres_logmov ( 0.5 + random-normal 0 0.25 ) 
  let pres_mov ( 10 ^ pres_logmov )
  set pres_angle random-normal 0 40
  if ( abs pres_angle > 60 ) [ set pres_angle (1 + random-float 0.5) * pres_angle ]  ; make angle dist more leptokurtic
  right pres_angle
  ;
  ; Turn to avoid swimming on land if necessary:
  ; ( section copied in from porp_std_move)
  let dd ceiling ( pres_mov / 0.25 )  ; number of 25-m steps to check water depth at
  let goto_avoid_land false
  if (not ( [ bathymetry ] of patch-ahead pres_mov >= min_depth ) ) [ set goto_avoid_land true ]
  repeat dd [
    if ( not ( [ bathymetry ] of patch-ahead ( dd * 0.25 ) >= min_depth ) ) [ set goto_avoid_land true ]  ; must use "not >= " rather than " < " for catching NaN's
    set dd dd - 1
  ]
  if ( goto_avoid_land ) [ porp_avoid_land ]
  set pres_mov ( 10 ^ pres_logmov )  ; because pres_logmov may have changed in the porp_avoid_land procedure
  ;  let ready-to-move true
  ; test again:
  set dd ceiling ( pres_mov / 0.1 )  ; number of 10-m steps to check water depth at
  let ee 0
  let depth-list list (0) ( [bathymetry] of patch-ahead pres_mov )
  set depth-list remove-item 0 depth-list
  repeat dd [
    set ee ee + 1
    set depth-list lput ( [ bathymetry ] of patch-ahead ( ee * 0.1 ) ) depth-list
  ]
  if ( not (( length depth-list ) = length ( filter [ ? > 0 ] depth-list )) ) [ ; i.e. if some items on the list aren't < 0
    uphill bathymetry
    if ( debug = 1 ) [ 
      show word "Tick = " my_tick
      show word "Moved to deeper patch, depth = " ([bathymetry] of patch-here) 
    ]
  ]
  ;
  ; move
  ;  if (ready-to-move) [ 
    fd pres_mov
  ;     ]
  ; Remember current moves for the next iteration
  set pres_logmov log pres_mov 10
  set prev_angle pres_angle
  set prev_logmov pres_logmov
end



to porp_avoid_land
  ; If shallow water ahead, turn right or left depending on where water is deeper. Turn as little as possible.
  ; Don't do the turning here, but change angle to be turned in porp_std_move or porp-markov-mov.
  ; Note that the emergency procedure "avoid_beh 5" is found in porp_std_move
  let rand_ang random 10
  let avoid_beh 0
  let pres_mov ( 10 ^ pres_logmov )
  let bath_l [ bathymetry ] of patch-left-and-ahead (40 + rand_ang) pres_mov
  let bath_r [ bathymetry ] of patch-right-and-ahead (40 + rand_ang) pres_mov
  ; alternative kinds of evasive behaviour: 
  if ( bath_r >= min_depth or bath_l >= min_depth ) [
    set avoid_beh 1  ; evasive behaviour type 1
    ifelse ( bath_r >= min_depth and bath_l >= min_depth ) 
      [ ifelse ( bath_r >= bath_l ) ; comparison can be true only if neither bath_r or bath_l are NaN, i.e. if both are > min_depth
        [ set pres_angle pres_angle + (40 + rand_ang) ]
        [ set pres_angle pres_angle - (40 + rand_ang) ]
      ]
      [ ifelse ( bath_r >= min_depth ) 
        [ set pres_angle pres_angle + (40 + rand_ang) ]
        [ set pres_angle pres_angle - (40 + rand_ang) ]
      ]
  ]
  ; else try turning more aprubtly ( = 70 deg )
  if not ( bath_r >= min_depth or bath_l >= min_depth ) [
    set avoid_beh 2  ; evasive behaviour type 2
    set bath_l [ bathymetry ] of patch-left-and-ahead (70 + rand_ang) pres_mov
    set bath_r [ bathymetry ] of patch-right-and-ahead (70 + rand_ang) pres_mov
    if ( bath_r >= min_depth or bath_l >= min_depth ) [
      ifelse ( bath_r >= min_depth and bath_l >= min_depth ) 
        [ ifelse ( bath_r >= bath_l ) ; comparison can be true only if neither bath_r or bath_l are NaN, i.e. if both are > min_depth
          [ set pres_angle pres_angle + (70 + rand_ang) ]
          [ set pres_angle pres_angle - (70 + rand_ang) ]
        ]
        [ ifelse ( bath_r >= min_depth ) 
          [ set pres_angle pres_angle + (70 + rand_ang) ]
          [ set pres_angle pres_angle - (70 + rand_ang) ]
        ]
    ]
  ]
  if not ( bath_r >= min_depth or bath_l >= min_depth ) [
    set avoid_beh 3  ; evasive behaviour type 3
    set bath_l [ bathymetry ] of patch-left-and-ahead (120 + rand_ang) pres_mov
    set bath_r [ bathymetry ] of patch-right-and-ahead (120 + rand_ang) pres_mov
    if ( bath_r >= min_depth or bath_l >= min_depth ) [
      ifelse ( bath_r >= min_depth and bath_l >= min_depth ) 
        [ ifelse ( bath_r >= bath_l ) ; comparison can be true only if neither bath_r or bath_l are NaN, i.e. if both are > min_depth
          [ set pres_angle pres_angle + (120 + rand_ang) ]
          [ set pres_angle pres_angle - (120 + rand_ang) ]
        ]
        [ ifelse ( bath_r >= min_depth ) 
          [ set pres_angle pres_angle + (120 + rand_ang) ]
          [ set pres_angle pres_angle - (120 + rand_ang) ]
        ]
    ]
  ]  
  if not ( bath_r >= min_depth or bath_l >= min_depth ) [
    ; if everything else fails, turn around
    set avoid_beh 4  ; evasive behaviour type 4
    let j 0
    porp_check_depth
    while [ not enough-water-ahead and j < length pos_list ] [
      facexy (item 0 (item j pos_list)) (item 1 (item j pos_list))  ; each item on pos_list contains a list with a x and a y-coordinate
      setxy (item 0 (item j pos_list)) (item 1 (item j pos_list))
      set j j + 1
      porp_check_depth
      if (j = 20) [ set enough-water-ahead true ]
    ]
  ]
  if ( debug = 1 ) [ 
    let tmp-list list ("beh =") avoid_beh 
    set tmp-list lput ("; tck =") tmp-list
    set tmp-list lput my_tick tmp-list
    write tmp-list 
    let tmp7 word "; " round pres_angle
    print word tmp7 " degr." 
  ]
end


to porp_check_depth
  ; Check that there is enough water at all steplengths ahead, set enough-water-ahead to false if < min_depth
  set enough-water-ahead true
  let pres_mov ( 10 ^ pres_logmov )                                                                           ; because pres_logmov may have changed in the porp_avoid_land procedure
  let dd ceiling ( pres_mov / 0.1 )                                                                           ; number of 10-m steps to check water depth at
  let ee 0
  let depth-list list (0) ( [ bathymetry ] of patch-ahead pres_mov )
  set depth-list remove-item 0 depth-list
  repeat dd [
    set ee ee + 1
    set depth-list lput ( [ bathymetry ] of patch-ahead ( ee * 0.1 ) ) depth-list
  ]
  if ( not (( length depth-list ) = length ( filter [ ? > 0 ] depth-list )) ) [                               ; i.e. if some depths on the list aren't > 0
    set enough-water-ahead false
  ]
end


; Main

to go

  if ( my_tick = 0 ) [
    reset-timer
    if ( write_data ) [ 
      file_setup
      print word "Writing to: " outfile
    ]
    if (debug = 1) [ 
      print ""
      print "Debug avoidance behaviour near land:"
    ]

    if (debug = 2) [ 
      print ""
      print "Write turning angles before/after react. to land:"
    ]
    if (debug = 3) [ 
      print ""
      print "Debugging CRW related turns (mod >=1):"
    ]
    if ( (debug = 4) and (model >= 2) ) [ 
      print ""
      print "Debugging attraction vector (mod >=2):"
    ]
    if ( (debug = 5) and (model >= 2) ) [ 
      print ""
      print "Debugging deterrence vector (mod >=2):"
    ]
    if ( (debug = 5) and (model >= 2) ) [ 
      print ""
      print "Debugging attraction + deterrence vectors (mod >=2):"
    ]
  ]

  ask porps [
    if write_data [ ; must write deployment pos before moving
      ask ( porp 0 ) [ file_write_line ] 
    ]
    porp_move
  ]

  my_update_plots
  set my_tick my_tick + 1
  set days my_tick / 48
  set time my_tick / 2     ; updates at half-hour intervals
  ;tick                    ; slows things down a lot


  ; make amount of food grow logistically (only daily update to speed things up, but done 48 times (once every 30-min) to make food grow right). Updated 2011-12-05
  if ( remainder days food-upd-interval ) = 0 [
    ask patches with [ food-prob01 > 0 and food_level < maxU ] [
      if ( food_level < 0.01 ) [ set food_level 0.01 ]    ; The minimum food level has a strong impact on how fast food gets back
      let f_lev food_level + ( food_growth_rate * food_level * ( 1 - food_level / maxU ) )
      if (abs (f_lev - food_level) > 0.001) [
        repeat 47 [ 
          set f_lev f_lev + ( food_growth_rate * food_level * ( 1 - food_level / maxU ) )
        ]
      ]  ; If the food level is really low, let food grow 48 times -- like growing every half-hour step, only faster
      set food_level f_lev
    ]
    if ( not track ) [ landsc_display ] ; updates the colours
  ]


  if ( my_tick = 14999 ) [ ; 15000 half-hour intervals is slightly longer track than recorded for pttid 2000-04542
    let tmp word "Time (15000 ticks): " timer
    print word tmp "sec"
    let fnm word ( random 1000 ) ".png"
    export-interface fnm
    print " " ; blank line
    if (write_data) [ file-close ]
    stop
  ]

end
