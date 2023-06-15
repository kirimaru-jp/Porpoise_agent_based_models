
#region "Test"

using BenchmarkTools

# [https://julialang.org/blog/2022/08/julia-1.8-highlights/#typed_globals]
# global x::Float64 = 1
global x = 1
function f1()
	global x = 0
end

@code_warntype f1()
x

# global x = 2.8
x::Float64 = 2.8

function f()
	s = 0.0
	for val in x
		s = s + val
	end
	return s
end
@code_warntype f()

@btime f()

#############################


# Struct definition
struct Foo
	w::Float64
	y::Int
	x::Tuple{Int32,UInt32}
end
a = Array{Foo}(undef, 5);
a[2]
a .= @set a[2].w = 5.0

#############################


y = [0.0 0; 0 0]

function change!(x)
	x[2,2] = 1
end

change!(y)
y

#endregion


#region "Load packages"

using CSV, Tables, BenchmarkTools
using Distributions, Random

# [https://discourse.julialang.org/t/reading-big-ascii-files/22697/3]
# [https://discourse.julialang.org/t/reading-huge-csv-files/15924/5]
# [https://stackoverflow.com/questions/54410030/read-csv-into-array]
# [https://discourse.julialang.org/t/redefine-struct-when-working-with-repl/25942/2]


"""
    read_asc(fil, skip, n_col)
	Read asc file into a matrix data:
		- fil: file path
		- skip: number of header line to skip
		- n_col: number of column to read
"""
function read_asc(fil, skip, n_col)
	f = CSV.File(fil, skipto=skip, delim=' ') |> Tables.matrix
	return Float64.(f[:,1:n_col])
end

# Read relevant environmental data
bathy_data::Matrix{Float64} = read_asc("raster-data/greatbelt/bathy.asc", 7, 1000)
food_prob01_data::Matrix{Float64} = read_asc("raster-data/greatbelt/food-prob01.asc", 7, 1000)
salinity_data::Matrix{Float64} = read_asc("raster-data/greatbelt/salinity.asc", 7, 1000)
sediment_data::Matrix{Float64} = read_asc("raster-data/greatbelt/sediment.asc", 7, 1000)
bot_ns_data::Matrix{Float64} = read_asc("raster-data/greatbelt/bot_ns.asc", 7, 1000)
bot_ew_data::Matrix{Float64} = read_asc("raster-data/greatbelt/bot_ew.asc", 7, 1000)


# Create the list of underwater patches
##TODO: avoid using temporal Vector{Any} p; currently don't know how to push! to Vector{CartesianIndex{2}}
# [https://discourse.julialang.org/t/eachindex-for-individual-axes-of-multidimensional-arrays/81264/2]
p = []
for i in axes(bathy_data, 1)
	for j in axes(bathy_data, 2)
		if bathy_data[i, j] > 0
			push!(p, CartesianIndex(i,j))
		end
	end
end
patches::Vector{CartesianIndex{2}} = Vector{CartesianIndex{2}}(undef, length(p))
for i in eachindex(p)
	patches[i] = p[i]
end


#endregion


#region "Global variables/porps_setup()"

## Several global variables

food_upd_interval::Int64 = 10			# food update interval
food_growth_rate::Float64 = 0.2			# food growth rate
maxU::Float64 = 1.7						# maximum food level
min_depth::Float64 = 1.0 				# min depth that porpoise can swim
n_porp::Int64 = 4						# number of porpoise
inertia_const::Float64 = 0.001			# Inertia constant; the animal’s tendency to keep moving using CRW irrespective of foraging success.
Nx::Int64, Ny::Int64 = size(bathy_data)	# dimension of bathymetric data


corr_logmov::Float64 = 0.94			# correlation in movement distance in CRW +
corr_angle::Float64 = 0.26			# correlation in direction in CRW +
vt::Vector{Float64} = [0.0, 0.0]	# resultant attraction vector, resulting from reference memory of food availability (model >=2)

use_exp_food_val::Bool= false		# get more attracted to the CRW path if food was found recently


# [https://discourse.julialang.org/t/how-to-obtain-indices-of-an-array-satisfying-boolean-condition/37780/19]
# [https://stackoverflow.com/questions/74562790/what-is-the-point-of-copy]
#
#; update food level -- also done in landsc_setup, but convenient to repeat it here
# ask patches [ 
#	 ifelse food-prob01 > 0 [ set food_level maxU ] [ set food_level food-prob01 ] 
# ]
food_level::Matrix{Float64} = copy(food_prob01_data)
food_level[findall(food_prob01_data .> 0.0)] .= maxU;


# global my_tick = 0
# global days = 0
# global time = 0						# hours since simulation start

# global path							# path to directory for input/output, one dir for each area
# global outfile						# simulation output

# global total_food_eaten			 	# Used for testing whether optimal foraging strategy is used
# global xllcorner 
# global yllcorner
# global deploy_x						# x-coord of porp 0
# global deploy_y						# y-coord of porp 0

# global turn_right = 1				# random variable# turn right if = 1, left if = -1
# global movlgt_list = [0.0, 0]	 	# list of 100 move lengths for porpoise 0 ; two inputs required...	list used for histogramming

# global angle_list					# list of 100 turning angles for porpoise 0

# global max_movlgt = 0.0				# monotonously increasing, move length/100m
# global mean_movlgt = 0.0
# global use_exp_food_val = false	 	# get more attracted to the CRW path if food was found recently
# global CRW_contrib = -9				# length of vector pointing in direction predicted by CRW
# global MR_contrib = -9				# length of vector pointing in direction of remembered food


# set corr-logmov 0.94
# set corr-angle 0.26
# set vt list 0 0
# set my-tick 0
# set days 0
# set time 0
# set max-movlgt 0
# set mean-movlgt 0
# set use-exp-food-val false
# set CRW-contrib -9
# set MR-contrib -9 
# set min-depth 1							 ; water depth where porps very rarely swim (according to Jakob)
# set turn-right 1
# set movlgt-list list (0) (0)	; two inputs required...	list used for histogramming
# set movlgt-list remove-item 0 movlgt-list
# set movlgt-list remove-item 0 movlgt-list
# set angle-list list (0) (0)	# two inputs required...	list used for histogramming
# set angle-list remove-item 0 angle-list
# set angle-list remove-item 0 angle-list


#endregion


#region "Porpoise & patches" 

mutable struct Patch
	bathymetry::Float64
	salinity::Float64
	sediment::Float64
	bot_ns_current::Float64
	bot_ew_current::Float64
	food_prob01::Float64		#; probability of finding food in a patch
	food_level::Float64			#; amount of food in a patch
end


porps_own [ 
	# energy-level					 #; Porpoises get energy by eating and loose energy by moving.
	# porp_length
	# prev_angle						 #; Last turning angle (not the heading!)
	# pres_angle						 #; Present turning angle
	# prev_logmov						#; Previous Log10 (move length /100 m)
	# pres_logmov						#; Present Log10 (move length /100 m)
	# enough_water_ahead		 #; Turn to avoid land if false
	# pos_list							 #; Coordinates of previous positions -- latest positions in left end
	# ptt										#; Argos id-number -- the number of the simulated porpoise
	# sex
	# p-length
	# p-weight
	# # ; Vars added in model 2
	# ref-mem-strength-list	#; Memory decay with time (logistic decrease, function of rR)
	# work-mem-strength-list #; Memory of poor quality patches, decays with time (logistic decrease, function of rW)
	# work-mem-updated			 #; Whether a list of working memory has been updated
	# stored-util-list			 #; Up_t ; Remembered feeding success (after memory decay) -- latest positions left
	# VE_total							 #; total value of food expected to be found in the future
]


#endregion

18 ÷ 4

#region "porps_setup()"
# [https://discourse.julialang.org/t/initialize-a-variable-inside-a-struct/43967/7]

using StaticArrays, Setfield, Parameters, LinkedLists

memory_max::Int64 = 325 # maximum length of memory coodinates
# memory_max::Int64 = 5
@with_kw struct Porp
	xcor::Float64
	ycor::Float64
	heading::Float64 = 0
	prev_angle::Float64
	pres_angle::Float64
	prev_logmov::Float64
	pres_logmov::Float64
	enough_water_ahead::Bool
	VE_total::Float64
	# pos_list::Array{CartesianIndex{2}, memory_max}
	# pos_list_x::SVector{memory_max,Float64}
	# pos_list_y::SVector{memory_max,Float64}
	# pos_list_x::LinkedList{Float64}
	# pos_list_y::LinkedList{Float64}
end

## Init 
## [https://discourse.julialang.org/t/how-to-define-a-function-such-that-its-arguments-can-mutate/66320/2]
## [https://stackoverflow.com/a/51760506/16077905]

Random.seed!(814);
porps::Vector{Porp} = Array{Porp}(undef, n_porp)
# pos_list = [LinkedList{CartesianIndex{2}}() for i in 1:n_porp]
pos_list::Vector{LinkedList{Tuple{Float64, Float64}}} = [LinkedList{NTuple{2,Float64}}() for i in 1:n_porp]


"""
    get_linked_list_value(L, i): get ith element of linked list L
"""
function get_linked_list_value(L, i)
	return getindex(L, positiontoindex(i, L))
end


"""
    porps_setup!(porps, pos_list, patches, n_porp): setup initial paramters for harbour porpoise argents
		- porps: list of Porpoises
		- pos_list: coordinates of previous positions -- latest positions at the end of list
		- patches: list of underwater patches
		- n_porp: number of porpoise agents
"""
function porps_setup!(porps, pos_list, patches, n_porp)
	idx = rand(patches, n_porp)
	for i in 1:n_porp
		porps .= @set porps[i].xcor = idx[i][1]
		porps .= @set porps[i].ycor = idx[i][2]
		porps .= @set porps[i].heading = 0.0
		porps .= @set porps[i].prev_angle = 0.0
		porps .= @set porps[i].pres_angle = 0.0
		porps .= @set porps[i].prev_logmov = 0.0
		porps .= @set porps[i].pres_logmov = 0.0
		porps .= @set porps[i].enough_water_ahead = true
		porps .= @set porps[i].VE_total = 0.0

		# porps .= @set porps[i].pos_list_x[1] = porps[i].xcor
		# porps .= @set porps[i].pos_list_y[1] = porps[i].ycor
		# push!(pos_list[i], CartesianIndex(porps[i].xcor, porps[i].ycor))
		empty!(pos_list[i])
		push!(pos_list[i], (porps[i].xcor, porps[i].ycor))
	end
	return
end


porps_setup!(porps, pos_list, patches, n_porp)
# @code_warntype porps_setup!(porps, pos_list, patches, n_porp)
@btime porps_setup!($porps, $pos_list, $patches, $n_porp)  ##   3.325 μs (46 allocations: 17.88 KiB)


# # Set 1st porpoise 
# porps .= @set porps[1].xcor = patches[1][1]
# porps .= @set porps[1].ycor = patches[1][2]

# push!(pos_list[1], CartesianIndex(1,2))
# push!(pos_list[1], CartesianIndex(11,2))
# last(pos_list[1])
# popfirst!(pos_list[1])
# # keys(pos_list[1])
# get_linked_list_value(pos_list[1], 1)


#endregion


#region "function porp_avoid_land()"

"""
    bound(X, Y): ensure that coordinate (X, Y) lies inside model world
"""
function bound(X, Y)
	return 1 .+ [mod(X, Nx), mod(Y, Ny)]
end


"""
    patch_ahead!(p0, x, y, d, heading): reports the single patch that is the given distance ahead of the asking porpoise
		- p0: coordinate of patch ahead
		- x: coordinate x
		- y: coordinate y
		- heading: heading of porpoise
"""
function patch_ahead!(p0, x, y, d, heading)
	dx = d * sin(heading * pi/180)
	dy = d * cos(heading * pi/180)
	X = Int64(floor(x + dx))
	Y = Int64(floor(y + dy))
	p0 .= bound(X, Y)
	return
end

"""
    patch_left_and_ahead(p::Porp, d, ang): reports the single patch that is the given distance from this porpoise, in the direction turned left the given angle (in degrees) from the turtle's current heading.
	- p: porpoise
	- d: distance ahead
	- ang: turned given angle
"""
function patch_left_and_ahead(p::Porp, d, ang)
	x = p.xcor
	y = p.ycor
	heading = p.heading

	dx = d * sin((heading-ang) * pi/180)
	dy = d * cos((heading-ang) * pi/180)
	# return [Int64(floor(x + dx)), Int64(floor(y + dy))]	
	X = Int64(floor(x + dx))
	Y = Int64(floor(y + dy))
	return bound(X, Y)

end

"""
    patch_right_and_ahead(p::Porp, d, ang): reports the single patch that is the given distance from this porpoise, in the direction turned right the given angle (in degrees) from the turtle's current heading.
	- p: porpoise
	- d: distance ahead
	- ang: turned given angle
"""
function patch_right_and_ahead(p::Porp, d, ang)
	x = p.xcor
	y = p.ycor
	heading = p.heading

	dx = d * sin((heading+ang) * pi/180)
	dy = d * cos((heading+ang) * pi/180)
	# return [Int64(floor(x + dx)), Int64(floor(y + dy))]	
	X = Int64(floor(x + dx))
	Y = Int64(floor(y + dy))
	return bound(X, Y)

end

"""
    porp_avoid_land!(porp, k): let kth porpoise change direction if land ahead
	; If shallow water ahead, turn right or left depending on where water is deeper. Turn as little as possible.
	; Don't do the turning here, but change angle to be turned in porp_std_move or porp-markov-mov.
	; Note that the emergency procedure "avoid_beh 5" is found in porp_std_move
"""
function porp_avoid_land!(porp, k)
	rand_ang = rand(Uniform(0, 10), 1)[1]
	avoid_beh = 0
	pres_mov = 10 ^ porp[k].pres_logmov

	r_ang = 40 + rand_ang
	pat = patch_left_and_ahead(porp[k], pres_mov, r_ang)
	bath_l = bathy_data[pat[1], pat[2]]
	pat = patch_right_and_ahead(porp[k], pres_mov, r_ang)
	bath_r = bathy_data[pat[1], pat[2]]

	pres_angle = porp[k].pres_angle
	#; alternative kinds of evasive behaviour: 
	if bath_r >= min_depth || bath_l >= min_depth
		avoid_beh = 1	#; evasive behaviour type 1
		if bath_r >= min_depth || bath_l >= min_depth
			if bath_r >= bath_l #; comparison can be true only if neither bath_r or bath_l are NaN, i.e. if both are > min_depth
				pres_angle += r_ang
			else
				pres_angle -= r_ang
			end
		else
			if bath_r >= min_depth
				pres_angle += r_ang
			else
				pres_angle -= r_ang
			end
		end
	end
	#; else try turning more aprubtly ( = 70 deg )
	if !(bath_r >= min_depth || bath_l >= min_depth)
		avoid_beh = 2	#; evasive behaviour type 2

		r_ang = 70 + rand_ang
		p = patch_left_and_ahead(porp[k], pres_mov, r_ang)
		bath_l = bathy_data[p[1],p[2]]
		p = patch_right_and_ahead(porp[k], pres_mov, r_ang)
		bath_r = bathy_data[p[1],p[2]]

		if bath_r >= min_depth || bath_l >= min_depth
			if bath_r >= bath_l #; comparison can be true only if neither bath_r or bath_l are NaN, i.e. if both are > min_depth
				pres_angle += r_ang
			else
				pres_angle -= r_ang
			end
		else
			if bath_r >= min_depth
				pres_angle += r_ang
			else
				pres_angle -= r_ang
			end
		end
	end
	#; else try turning more aprubtly ( = 120 deg )
	if !(bath_r >= min_depth || bath_l >= min_depth)
		avoid_beh = 3	#; evasive behaviour type 3

		r_ang = 120 + rand_ang
		p = patch_left_and_ahead(porp[k], pres_mov, r_ang)
		bath_l = bathy_data[p[1],p[2]]
		p = patch_right_and_ahead(porp[k], pres_mov, r_ang)
		bath_r = bathy_data[p[1],p[2]]

		if bath_r >= min_depth || bath_l >= min_depth
			if bath_r >= bath_l #; comparison can be true only if neither bath_r or bath_l are NaN, i.e. if both are > min_depth
				pres_angle += r_ang
			else
				pres_angle -= r_ang
			end
		else
			if bath_r >= min_depth
				pres_angle += r_ang
			else
				pres_angle -= r_ang
			end
		end
	end
	porp .= @set porp[k].pres_angle = pres_angle

	#; if everything else fails, turn around
	if bath_r >= min_depth || bath_l >= min_depth
		avoid_beh = 4	#; evasive behaviour type 1
		porp_check_depth!(porp, k)

		j = 1
		while !porp[k].enough_water_ahead && j < length(pos_list[k])
			pos = get_linked_list_value(pos_list[k], j)
			if pos[1] != porp[k].xcor
				porp .= @set porp[k].heading = atan((pos[1] - porp[k].xcor)/(pos[2] - porp[k].ycor)) * 180/pi + 90
			end
			porp .= @set porp[k].xcor = pos[1]
			porp .= @set porp[k].ycor = pos[2]
			j += 1

			porp_check_depth!(porp, k)
			if j == 20
				porp .= @set porp[k].enough_water_ahead = true
			end

		end
	end

	return
end

"""
	porp_check_depth!(porp, k): Check that there is enough water at all steplengths ahead for kth porpoise, set enough_water_ahead to false if < min_depth
"""
function porp_check_depth!(porp, k)

	porp .= @set porp[k].enough_water_ahead = true
	pres_mov = 10 ^ porp[k].pres_logmov		 #; because pres_logmov may have changed in the porp_avoid_land procedure
	dd = Int64(ceil( pres_mov / 0.1 ))		 #; number of 10-m steps to check water depth at
	p1 = [1, 1]
	for ee = 1:dd
		patch_ahead!(p1, porp[k].xcor, porp[k].ycor, ee*0.1, porp[k].heading)
		if !(bathy_data[p1[1], p1[2]] > 0.0)
			porp .= @set porp[k].enough_water_ahead = false
			break
		end
	end

	return
end


# @code_warntype patch_left_and_ahead(porp, 10, 20)
@code_warntype porp_check_depth!(porps, 1)
@code_warntype porp_avoid_land!(porps, 3)


#################################

to porp_avoid_land
	#; If shallow water ahead, turn right or left depending on where water is deeper. Turn as little as possible.
	#; Don't do the turning here, but change angle to be turned in porp_std_move or porp-markov-mov.
	#; Note that the emergency procedure "avoid_beh 5" is found in porp_std_move
	let rand_ang random 10
	let avoid_beh 0
	let pres_mov ( 10 ^ pres_logmov )
	let bath_l [ bathymetry ] of patch-left-and-ahead (40 + rand_ang) pres_mov
	let bath_r [ bathymetry ] of patch-right-and-ahead (40 + rand_ang) pres_mov
	#; alternative kinds of evasive behaviour: 
	if ( bath_r >= min_depth or bath_l >= min_depth ) [
		set avoid_beh 1	#; evasive behaviour type 1
		ifelse ( bath_r >= min_depth and bath_l >= min_depth ) 
			[ ifelse ( bath_r >= bath_l ) #; comparison can be true only if neither bath_r or bath_l are NaN, i.e. if both are > min_depth
				[ set pres_angle pres_angle + (40 + rand_ang) ]
				[ set pres_angle pres_angle - (40 + rand_ang) ]
			]
			[ ifelse ( bath_r >= min_depth ) 
				[ set pres_angle pres_angle + (40 + rand_ang) ]
				[ set pres_angle pres_angle - (40 + rand_ang) ]
			]
	]
	#; else try turning more aprubtly ( = 70 deg )
	if not ( bath_r >= min_depth or bath_l >= min_depth ) [
		set avoid_beh 2	#; evasive behaviour type 2
		set bath_l [ bathymetry ] of patch-left-and-ahead (70 + rand_ang) pres_mov
		set bath_r [ bathymetry ] of patch-right-and-ahead (70 + rand_ang) pres_mov
		if ( bath_r >= min_depth or bath_l >= min_depth ) [
			ifelse ( bath_r >= min_depth and bath_l >= min_depth ) 
				[ ifelse ( bath_r >= bath_l ) #; comparison can be true only if neither bath_r or bath_l are NaN, i.e. if both are > min_depth
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
		set avoid_beh 3	#; evasive behaviour type 3
		set bath_l [ bathymetry ] of patch-left-and-ahead (120 + rand_ang) pres_mov
		set bath_r [ bathymetry ] of patch-right-and-ahead (120 + rand_ang) pres_mov
		if ( bath_r >= min_depth or bath_l >= min_depth ) [
			ifelse ( bath_r >= min_depth and bath_l >= min_depth ) 
				[ ifelse ( bath_r >= bath_l ) #; comparison can be true only if neither bath_r or bath_l are NaN, i.e. if both are > min_depth
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
		#; if everything else fails, turn around
		set avoid_beh 4	#; evasive behaviour type 4
		let j 0
		porp_check_depth
		while [ not enough_water_ahead and j < length pos_list ] [
			facexy (item 0 (item j pos_list)) (item 1 (item j pos_list))	#; each item on pos_list contains a list with a x and a y-coordinate
			setxy (item 0 (item j pos_list)) (item 1 (item j pos_list))
			set j j + 1
			porp_check_depth
			if (j = 20) [ set enough_water_ahead true ]
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
	#; Check that there is enough water at all steplengths ahead, set enough_water_ahead to false if < min_depth
	set enough_water_ahead true
	let pres_mov ( 10 ^ pres_logmov )		 #; because pres_logmov may have changed in the porp_avoid_land procedure
	let dd ceiling ( pres_mov / 0.1 )		 #; number of 10-m steps to check water depth at
	let ee 0
	let depth-list list (0) ( [ bathymetry ] of patch-ahead pres_mov )
	set depth-list remove-item 0 depth-list
	repeat dd [
		set ee ee + 1
		set depth-list lput ( [ bathymetry ] of patch-ahead ( ee * 0.1 ) ) depth-list
	]
	if ( not (( length depth-list ) = length ( filter [ ? > 0 ] depth-list )) ) [		#; i.e. if some depths on the list aren't > 0
		set enough_water_ahead false
	]
end


#endregion


#region "porp_move: model1 -> porp_markov_move()"

# porp = porps[2]

"""
    porp_markov_move!(porp, k): movement model of kth porpoise based on intrinsic exploratory behavior. Behaviour produced by simple non-autocorrelated Markov process (except when porps seek to avoid land) variables calibrated using dead-reckoning data. 
"""
function porp_markov_move!(porp, k)
	# Movements based on dead-reckoning data -- first calc distance, then turning angle
	pres_logmov = 0.5 + rand(Normal(0, 0.25), 1)[1]
	pres_angle = rand(Normal(0, 40), 1)[1]
	if abs(pres_angle) > 60 # make angle dist more leptokurtic
		pres_angle = (1 + rand(Uniform(0, 0.5), 1)[1]) * pres_angle
	end
	# right pres_angle
	heading = mod(porp[k].heading + pres_angle, 360)
	
	# Turn to avoid swimming on land if necessary
	goto_avoid_land = false
	p1 = [1,1]
	pres_mov = 10 ^ porp[k].pres_logmov
	patch_ahead!(p1, porp[k].xcor, porp[k].ycor, pres_mov, heading)
	if !(bathy_data[p1[1], p1[2]] >= min_depth)
		goto_avoid_land = true
	end
	# number of 25-m steps to check water depth at
	dd = Int64(ceil(pres_mov / 0.25))
	for d = dd:-1:1 
		patch_ahead!(p1, porp[k].xcor, porp[k].ycor, d * 0.25, heading)
		if !(bathy_data[p1[1], p1[2]] >= min_depth)
			goto_avoid_land = true
		end	
	end
	
	porp[k] = Porp(porp[k].xcor, porp[k].ycor, heading, porp[k].prev_angle, pres_angle, 
					porp[k].prev_logmov, pres_logmov, porp[k].enough_water_ahead, porp[k].VE_total)
	if goto_avoid_land
		porp_avoid_land!(porp, k)
	end
	
	pres_mov = 10 ^ porp[k].pres_logmov	#; because pres_logmov may have changed in the porp_avoid_land procedure
	#;	let ready-to-move true
	#; test again:
	dd = Int64(ceil( pres_mov / 0.1 ))	#; number of 10-m steps to check water depth at
	# depth_list = [ bathy_data[p₀[1], p₀[2]]	]
	
	xcor = porp[k].xcor
	ycor = porp[k].ycor
	heading = porp[k].heading
	for ee = 1:dd
		patch_ahead!(p1, xcor, ycor, ee*0.1, heading)
		if bathy_data[p1[1], p1[2]] < 0.0 # i.e. if some items on the list aren't < 0
			# uphill bathymetry
			x = Int64(floor(xcor))
			y = Int64(floor(ycor))
			x₀ = ifelse(x <= 1, 1, x-1)
			y₀ = ifelse(y <= 1, 1, y-1)
			x₁ = ifelse(x < Nx-1, x+1, Nx)
			y₁ = ifelse(y < Ny-1, y+1, Ny)

			q = argmax( bathy_data[x₀:x₁, y₀:y₁] )
			xcor = q[1]
			ycor = q[2]
			break
			# println(ee, ": ", bathy_data[p[1], p[2]])
		end
		# if ( debug = 1 ) [ 
		# 	show word "Tick = " my_tick
		# 	show word "Moved to deeper patch, depth = " ([bathymetry] of patch-here) 
		# ]
	end
	
	# move
	# fd pres_mov
	xcor += pres_mov * sin(heading * pi/180)
	ycor += pres_mov * cos(heading * pi/180)
	xcor, ycor = bound(xcor, ycor)
	# Remember current moves for the next iteration
	porp[k] = Porp(xcor, ycor, heading, porp[k].pres_angle, porp[k].pres_angle, 
					porp[k].pres_logmov, porp[k].pres_logmov, 
					porp[k].enough_water_ahead, porp[k].VE_total)
	
	return
		
end


@btime porp_markov_move!($porps, 1) ##   3.150 μs (54 allocations: 7.42 KiB)
# @code_warntype porp_avoid_land!(porps, 3)
# porp_markov_move!(porps, 1)


function all_porp_markov_move!(porps, n_porp, n_loop)
	for j = 1:n_loop
		for k = 1:n_porp
			porp_markov_move!(porps, k)
		end	
	end
	return
end

@btime all_porp_markov_move!($porps, $n_porp, 1500)
##   24.000 μs (314 allocations: 40.89 KiB): 1 loop
##  110.033 ms (59445 allocations: 6.19 MiB): 100 loops
##    2.937 ms (65117 allocations: 4.99 MiB): 150 loops
## 5.705 s (723250 allocations: 55.67 MiB): 1500 loops

@code_warntype all_porp_markov_move!(porps, n_porp, 150)


#################################

Threads.nthreads()
using .Threads: @threads

function th_all_porp_markov_move!(porps, n_porp, n_loop)
	for j = 1:n_loop
		@threads for k = 1:n_porp
			porp_markov_move!(porps, k)
		end	
	end
	return
end

@btime th_all_porp_markov_move!($porps, $n_porp, 1500)


function s_prefix_threads!(⊕, y::AbstractVector)
	l = length(y)
	k = ceil(Int, log2(l))
	# do reduce phase
	for j = 1:k
		for i = 2^j:2^j:min(l, 2^k)
			@inbounds y[i] = y[i - 2^(j - 1)] ⊕ y[i]
		end
	end
	# do expand phase
	for j = (k - 1):-1:1
		for i = 3*2^(j - 1):2^j:min(l, 2^k)
			@inbounds y[i] = y[i - 2^(j - 1)] ⊕ y[i]
		end
	end
	return y
end
function prefix_threads!(⊕, y::AbstractVector)
	l = length(y)
	k = ceil(Int, log2(l))
	# do reduce phase
	for j = 1:k
		@threads for i = 2^j:2^j:min(l, 2^k)
			@inbounds y[i] = y[i - 2^(j - 1)] ⊕ y[i]
		end
	end
	# do expand phase
	for j = (k - 1):-1:1
		@threads for i = 3*2^(j - 1):2^j:min(l, 2^k)
			@inbounds y[i] = y[i - 2^(j - 1)] ⊕ y[i]
		end
	end
	return y
end
A = fill(1, 500_000)
prefix_threads!(+, A)
@btime prefix_threads!(+, $A)
##   752.100 μs (938 allocations: 100.42 KiB)
@btime s_prefix_threads!(+, $A)
##   2.461 ms (0 allocations: 0 bytes)


#################################

to porp_markov_move
	#; Movements based on dead-reckoning data -- first calc distance, then turning angle
	set pres_logmov ( 0.5 + random-normal 0 0.25 ) 
	let pres_mov ( 10 ^ pres_logmov )
	set pres_angle random-normal 0 40
	if ( abs pres_angle > 60 ) [ set pres_angle (1 + random-float 0.5) * pres_angle ]	#; make angle dist more leptokurtic
	right pres_angle
	#;
	#; Turn to avoid swimming on land if necessary:
	#; ( section copied in from porp_std_move)
	let dd ceiling ( pres_mov / 0.25 )	#; number of 25-m steps to check water depth at
	let goto_avoid_land false

	if ( not( [ bathymetry ] of patch-ahead pres_mov >= min_depth ) ) [ 
	set goto_avoid_land true ]
	repeat dd [
		if ( not( [ bathymetry ] of patch-ahead ( dd * 0.25 ) >= min_depth ) ) [ 
		set goto_avoid_land true ]	# must use "not >= " rather than " < " for catching NaN's
		set dd dd - 1
	]

	if ( goto_avoid_land ) [ porp_avoid_land ]
	set pres_mov ( 10 ^ pres_logmov )	#; because pres_logmov may have changed in the porp_avoid_land procedure
	#;	let ready-to-move true
	#; test again:
	set dd ceiling ( pres_mov / 0.1 )	#; number of 10-m steps to check water depth at
	let ee 0
	let depth-list list (0) ( [bathymetry] of patch-ahead pres_mov )
	set depth-list remove-item 0 depth-list
	repeat dd [
		set ee ee + 1
		set depth-list lput ( [ bathymetry ] of patch-ahead ( ee * 0.1 ) ) depth-list
	]
	if ( not (( length depth-list ) = length ( filter [ ? > 0 ] depth-list )) ) [ #; i.e. if some items on the list aren't < 0
		uphill bathymetry
		if ( debug = 1 ) [ 
			show word "Tick = " my_tick
			show word "Moved to deeper patch, depth = " ([bathymetry] of patch-here) 
		]
	]

	# ;
	# ; move
	# ;	if (ready-to-move) [ 
	fd pres_mov
	# ;	]
	# ; Remember current moves for the next iteration
	set pres_logmov log pres_mov 10
	set prev_angle pres_angle
	set prev_logmov pres_logmov
end

#endregion


#region "function file_loop() & go()"

"""
    file_loop(n_loop): repeate the Markov movement model for all porpoises.
"""
function file_loop(n_loop)
	for i = 1:10
		porps_setup!(porps, pos_list, patches, n_porp)
		# for j = 1:15000
		for j = 1:n_loop
			for k = 1:n_porp
				porp_markov_move!(porps, k)
			end
			# my_update_plots
			my_tick = (i-1)*n_loop + j
			days = my_tick / 48
			time = my_tick / 2		 # updates at half-hour intervals
			#tick					 # slows things down a lot
	
			if mod(days, food_upd_interval) == 0
				
				for patch in patches
	
					patch_prob01 = food_level[patch]
					patch_food_level = food_prob01_data[patch]
					
					if patch_prob01 > 0 && patch_food_level < maxU
						if patch_food_level < 0.01
							patch_food_level = 0.01
						end
						f_lev = patch_food_level + 
								food_growth_rate * patch_food_level * ( 1 - patch_food_level / maxU )
						if abs(f_lev - patch_food_level) > 0.001
							for i = 1:47 
								f_lev += food_growth_rate * patch_food_level * ( 1 - patch_food_level / maxU )
							end
						end
						food_level[patch] = f_lev	
					end
	
				end
	
			end
		end
	end
	return
end


@btime file_loop(150)
# 33.000 μs (460 allocations: 178.75 KiB): only porps_setup!()
#  3.742 ms (460 allocations: 178.75 KiB): porps_setup!() + "mod(days, food_upd_interval) == 0"
# 38.030 ms (689852 allocations: 52.87 MiB): porps_setup!() + porp_markov_move!() + "mod(days, food_upd_interval) == 0"

@btime file_loop(1500)
# 518.903 ms (7103179 allocations: 545.70 MiB)

@btime file_loop(15000)
# 3.590 s (52943492 allocations: 4.08 GiB)

# @code_warntype file_loop(150)


#######################

	#to file_loop  #; corresponds to the wrt button -- makes all output files
	set write_data true
	repeat 10 [
		porps_setup
		repeat 15000 [  #; max my_tick, i.e. nearly a year
			go
		]
	]
	set write_data false
	end	
	# to go
	ask porps [
		if write_data [ # must write deployment pos before moving
			ask ( porp 0 ) [ file_write_line ] 
		]
		porp_move
	]
	my_update_plots
	set my_tick my_tick + 1
	set days my_tick / 48
	set time my_tick / 2     #; updates at half-hour intervals
	;tick                    #; slows things down a lot
	# ; make amount of food grow logistically (only daily update to speed things up, but done 48 times (once every 30-min) to make food grow right). Updated 2011-12-05
	if ( remainder days food-upd-interval ) = 0 [
		ask patches with [ food-prob01 > 0 and food_level < maxU ] [
			if ( food_level < 0.01 ) [ set food_level 0.01 ]		# The minimum food level has a strong impact on how fast food gets back
			let f_lev food_level + ( food-growth-rate * food_level * ( 1 - food_level / maxU ) )
			if (abs (f_lev - food_level) > 0.001) [
				repeat 47 [ 
					set f_lev f_lev + ( food-growth-rate * food_level * ( 1 - food_level / maxU ) ) 
				]
			]	# If the food level is really low, let food grow 48 times -- like growing every half-hour step, only faster
			set food_level f_lev
		]
		if ( not track ) [ landsc_display ] # updates the colours
	]
 
	if ( my_tick = 14999 ) [ # 15000 half-hour intervals is slightly longer track than recorded for pttid 2000-04542
		let tmp word "Time (15000 ticks): " timer
		print word tmp "sec"
		let fnm word ( random 1000 ) ".png"
		export-interface fnm
		print " " ; blank line
		if (write_data) [ file-close ]
		stop
	]

#endregion


#region "function porp_std_move()"


"""
    max_one_of_neighbors(xcor, ycor, B): report maximum bathymetric data among 8 neighbors around current porpoise agent
		- xcor, ycor: current coordinate of porpoise
		- B: bathymetric data matrix
"""
function max_one_of_neighbors(xcor, ycor, B)
	x = Int64(floor(xcor))
	y = Int64(floor(ycor))
	x₀ = ifelse(x <= 1, 1, x-1)
	y₀ = ifelse(y <= 1, 1, y-1)
	x₁ = ifelse(x < Nx-1, x+1, Nx)
	y₁ = ifelse(y < Ny-1, y+1, Ny)
	
	q = [x₀, y₀]
	if q[1] == x && q[2] == y
		if x₀ != x₁
			q[1] += 1
		end
		if y₀ != y₁
			q[2] += 1
		end
	end
	m = B[q[1], q[2]]
	for i = x₀:x₁
		for j = y₀:y₁
			if m < B[i,j] && !(i == x && j == y)
				q = [i,j]
				m = B[q[1], q[2]]
			end
		end
	end
	return q
end

# [https://github.com/NetLogo/NetLogo/blob/e35b512d129abee9bc3b81b4433d870826aa9c43/netlogo-core/src/main/agent/Turtle.java#L960]
"""
    subtractHeadings(h1, h2): Computes the difference between the given headings, that is, the number of degrees in the smallest angle by which heading h2 could be rotated to produce heading h1. A positive answer means a clockwise rotation, a negative answer counterclockwise.
"""
function subtractHeadings(h1, h2)
	if h1 < 0 || h1 >= 360
		h1 = mod(h1, 360)
	end
	if h2 < 0 || h2 >= 360
		h2 = mod(h2, 360)
	end
	diff = h1 - h2
	if diff > -180 && diff <= 180
		return diff
	elseif diff > 0
		return diff - 360
	else
		return diff + 360
	end
end

# k = 2


"""
    porp_std_move!(porps, k): movement model of kth porpoise like Markov model above, but introducing autocorrelation in turning angles and speed, and sharper turning angles at low speed (i.e. a correlated random walk (CRW) model). 
"""
function porp_std_move!(porps, k)
	
	porps .= @set porps[k].pres_angle = 999.0

	prev_logmov = porps[k].prev_logmov
	prev_mov = 10.0 ^ porps[k].prev_logmov
	pres_heading = porps[k].heading
	heading = porps[k].heading
	pres_angle = porps[k].pres_angle
	prev_angle = porps[k].prev_angle

	j = 1
	tmp_angle = 0.0
	#; for increasing mean turning angle
	if prev_angle < 0
		tmp_angle = prev_angle - 24.0
	else
		tmp_angle = prev_angle + 24.0
	end
	while abs(pres_angle) > 180.0
		# Autoreg can't be used for estimating param. as estimated turns are changed if on shallow water
		pres_angle = tmp_angle * (- corr_angle) + rand(Normal(0, 38), 1)[1]
		j += 1
		if j == 200
			pres_angle = 90 * sign(pres_angle)
		end
	end
	s = sign(pres_angle)
	pres_angle = abs(pres_angle) # add the sign again later
	#; Make angle decrease linearly with mov-dist
	go_on = true
	j = 1
	rnd = 0.0
	while go_on 
		rnd = rand(Normal(96, 28), 1)[1]      #; draws the number to be added to pres_angle
		if prev_mov <= 5.5
			pres_angle += rnd - rnd * prev_mov / 5.5
		end
		if pres_angle < 180
			go_on = false  #; remember that turning angle is unsigned here
		end
		j += 1
		if j == 200
			pres_angle = rand(0:19) + 90.0
			go_on = false
		end
	end
	#restore the sign of pres_angle
	pres_angle *= s
	angle_before_avoid_land = pres_angle #; for printing later using debug 2

	# right pres_angle
	heading += pres_angle
	angle_turned_right = pres_angle #; for updating prev_angle at end of porp_std_move
	pres_angle = 0.0
	#; ### distance
	pres_logmov = 999.0
	porp_max_dist = 1.18                 #; log10 ( max distance a porpoise can travel per half-hour )
	j = 1
	while pres_logmov > porp_max_dist 
		pres_logmov = corr_logmov * prev_logmov + rand(Normal(0.42, 0.48), 1)[1]
		j += 1
		if j == 200
			pres_angle = 90 * sign(pres_angle)
		end  
	end    #; Mean pres_mov should be x.x x100 m 
	pres_mov = 10.0^pres_logmov        #; This is what is plotted in the histogram
	#;
	#; Turn to avoid swimming on land if necessary:
	enough_water_ahead = false
	count_i = 0


	while !porps[k].enough_water_ahead
		porp_check_depth!(porps, k)
		if !porps[k].enough_water_ahead
			porp_avoid_land!(porps, k)
		end
		pres_mov = 10.0^porps[k].pres_logmov   	#; because pres_logmov may have changed in the porp_avoid_land procedure
		# right pres_angle     	    	#; angle to turn -- pres_angle -- is changed in porp_avoid_land
		porps .= @set porps[k].heading += porps[k].pres_angle
		angle_turned_right = mod(angle_turned_right + porps[k].pres_angle, 360) - 180
		porps .= @set porps[k].pres_angle = 0.0
		count_i += 1
		if count_i == 100
			porps .= @set porps[k].enough_water_ahead = true
		end
	end
	#; test depth again, avoid_beh = 5:
	porp_check_depth!(porps, k)

	xcor = porps[k].xcor
	ycor = porps[k].ycor
	heading = porps[k].heading
	pres_angle = porps[k].pres_angle
	enough_water_ahead = porps[k].enough_water_ahead


	if !enough_water_ahead
		prev_heading = heading
		# let p max-one-of neighbors [ bathymetry ]
		q = max_one_of_neighbors(xcor, ycor, bathy_data)
		# face p
		if q[2] != ycor
			heading = atan((q[1]-xcor)/(q[2]-ycor)) * 180/pi
		end
		angle_turned_right = mod(angle_turned_right + pres_angle, 360) - 180
		pres_mov = 1.0                                 #; move 100 m towards deeper patch
	end
	#; slow down if turning sharply:
	if pres_mov > 10 && abs(angle_turned_right) > 90
		pres_mov /= 5 
	end
	if pres_mov > 7 && abs(angle_turned_right) > 50
		pres_mov /= 2
	end 
	#; Change direction if attracted / repelled by certain areas (model >= 2)
	total_dx = 0.0
	total_dy = 0.0
	if !use_exp_food_val
		dx = sin(heading); dy = cos(heading)
		total_dx = dx * pres_mov + vt[1]      #; vt isn't used in porp_std_move till here
		total_dy = dy * pres_mov + vt[2]      #; note that dx is change in x if taking ONE step forward
		# facexy (xcor + total_dx) (ycor + total_dy)
		if total_dy != 0
			heading = atan(total_dx/total_dy) * 180/pi
		end
	else
		CRW_contrib = inertia_const + pres_mov * porps[k].VE_total  #; length of vector pointing in direction predicted by CRW
		MR_contrib = sqrt( vt[1]*vt[1] + vt[2]*vt[2] )     #; length of vector pointing in direction of remembered food
		dx = sin(heading); dy = cos(heading)
		total_dx = dx * CRW_contrib + vt[1]
		total_dy = dy * CRW_contrib + vt[2]
		# facexy (xcor + total_dx) (ycor + total_dy)        #; really not needed, it already points that way
		if total_dy != 0
			heading = atan(total_dx/total_dy) * 180/pi
		end
	end
	#; Store turn for calc of turning angle in next step:
	total_turn = subtractHeadings(heading, pres_heading)   #; total change in heading, including all adjustments till here. 'pres_heading' was calc in beginning of porp_std_move


	#;
	#; Move: 
	# fd pres_mov  #; movement length isn't affected by presence of food
	xcor += pres_mov * sin(heading * pi/180)
	ycor += pres_mov * cos(heading * pi/180)
	p1 = bound(xcor, ycor)
	xcor, ycor = p1[1], p1[2]
	#; Remember current moves for the next iteration
	#; if attraction to food alters the movement angle (i.e. vt != 0), this isn't remembered for next step
	#; set prev_angle angle_turned_right  ; so the additional turn due to attraction to food does not influence turning angle in next step 
	prev_angle = total_turn  #; so the additional turn due to attraction to food DOES influence turning angle in next step 
	prev_logmov = log10(pres_mov)  #; total steplength, resulting from vt + pres_mov
	#;
	#; test depth one last time, avoid_beh = 6 - move back on same track:

	if !( bathy_data[Int64(floor(xcor)), Int64(floor(ycor))] > 0 )
		prev_heading = heading
		# facexy (item 0 item 1 pos_list) (item 1 item 1 pos_list)
		pos = last(pos_list[k])
		if pos[2] != ycor
			heading = atan((pos[1] - xcor)/(pos[2] - ycor)) * 180/pi
		end
		angle_turned_right = mod(angle_turned_right + pres_angle, 360) - 180
		# setxy (item 0 item 1 pos_list) (item 1 item 1 pos_list)          #; move 100 m towards deeper patch
		xcor = pos[1]
		ycor = pos[2]
	end
	#; update position list:
	# pres_pos = (xcor, ycor)
	# pos_list = fput pres_pos pos_list
	push!(pos_list[k], (xcor, ycor))
	if length(pos_list[k]) > memory_max
		popfirst!(pos_list[k])
	end


	# porp[k] = Porp(xcor, ycor, heading, prev_angle, porp[k].pres_angle, 
	# 				prev_logmov, porp[k].pres_logmov, 
	# 				porp[k].enough_water_ahead, porp[k].VE_total)
	porps .= @set porps[k].xcor = xcor
	porps .= @set porps[k].ycor = ycor
	porps .= @set porps[k].heading = heading
	porps .= @set porps[k].prev_angle = prev_angle
	porps .= @set porps[k].prev_logmov = prev_logmov

end

# @code_warntype porp_std_move!(porps, 4)

@btime porp_std_move!($porps, 4)



# porps[4]


"""
    file_loop1(n_loop): repeate the standard movement model for all porpoises.
"""
function file_loop1(n_loop)
	for i = 1:10
		porps_setup!(porps, pos_list, patches, n_porp)
		# for j = 1:15000
		for j = 1:n_loop
			for k = 1:n_porp
				porp_std_move!(porps, k)
			end
			# my_update_plots
			my_tick = (i-1)*n_loop + j
			days = my_tick / 48
			time = my_tick / 2		 # updates at half-hour intervals
			#tick					 # slows things down a lot
	
			if mod(days, food_upd_interval) == 0
				
				for patch in patches
	
					patch_prob01 = food_level[patch]
					patch_food_level = food_prob01_data[patch]
					
					if patch_prob01 > 0 && patch_food_level < maxU
						if patch_food_level < 0.01
							patch_food_level = 0.01
						end
						f_lev = patch_food_level + 
								food_growth_rate * patch_food_level * ( 1 - patch_food_level / maxU )
						if abs(f_lev - patch_food_level) > 0.001
							for i = 1:47 
								f_lev += food_growth_rate * patch_food_level * ( 1 - patch_food_level / maxU )
							end
						end
						food_level[patch] = f_lev	
					end
	
				end
	
			end
		end
	end
	return
end


@btime file_loop1(150)
# 205.076 ms (528526 allocations: 74.52 MiB): porps_setup!() + porp_std_move!() + "mod(days, food_upd_interval) == 0"

@btime file_loop1(1500)
# 8.931 s (100520147 allocations: 16.36 GiB)

@btime file_loop1(1000)
# 5.995 s (74381535 allocations: 11.81 GiB)



####################

to porp_std_move
  #; Movements based on dead-reckoning data:
  #; global vars: corr_logmov 0.94 and corr_angle 0.26
  #; ### turning angle
  let prev_mov 10 ^ prev_logmov
  let pres_heading heading
  set pres_angle 999
  let j 1
  let tmp_angle 0
  ifelse ( prev_angle < 0 ) [ set tmp_angle prev_angle - 24 ] [ set tmp_angle prev_angle + 24 ]  #; for increasing mean turning angle
  while [ abs (pres_angle) > 180 ]  [ 
    set pres_angle ( tmp_angle * (- corr_angle) + random-normal 0 38 )                  #; Autoreg can't be used for estimating param. as estimated turns are changed if on shallow water. 
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
  set pres_angle abs pres_angle #; add the sign again later
  #; Make angle decrease linearly with mov-dist
  let go-on true
  set j 1
  let rnd 0
  while [ go-on ]  [ 
    set rnd random-normal 96 28      #; draws the number to be added to pres_angle
    if ( prev_mov <= 5.50 ) [ 
      set pres_angle pres_angle + rnd - ( rnd * prev_mov / 5.50 )
    ]
    if ( pres_angle < 180 ) [ set go-on false ]  #; remember that turning angle is unsigned here
    set j j + 1
    if (j = 200) [
      set pres_angle ( random 20 + 90 )
      set go-on false
      if ( debug = 3 ) [ 
        print word "exiting loop 2, ang=" pres_angle
      ]
    ]
  ]
  #; if ( abs pres_angle > 55 and abs pres_angle < 180 ) [ set pres_angle (1 - random-float 0.32) * pres_angle ]  ; make turning angle dist more leptokurtic
  set pres_angle pres_angle * sign
  let angle_before_avoid_land pres_angle #; for printing later using debug 2
  right pres_angle
  let angle_turned_right pres_angle #; for updating prev_angle at end of porp_std_move
  set pres_angle 0

  #; ### distance
  set pres_logmov 999
  let porp_max_dist 1.18                 #; log10 ( max distance a porpoise can travel per half-hour )
  set j 1
  while [ pres_logmov > porp_max_dist ] [ 
    set pres_logmov ( corr_logmov * prev_logmov + random-normal 0.42 0.48 ) 
    set j j + 1
    if (j = 200) [
      if (pres_angle = 0) [set pres_angle pres_angle + 0.00001]
      set pres_angle ( pres_angle * 90 / (abs pres_angle))
      if ( debug = 3 ) [ 
        print word "exiting loop 3, ang=" pres_angle
      ]

    ]    
  ]    #; Mean pres_mov should be x.x x100 m 
  let pres_mov ( 10 ^ pres_logmov )        #; This is what is plotted in the histogram
  #;
  #; Turn to avoid swimming on land if necessary:
  set enough_water_ahead false
  let count_i 0
  while [ not enough_water_ahead ] [
    porp_check_depth
    if (not enough_water_ahead) [ porp_avoid_land ]
    set pres_mov ( 10 ^ pres_logmov )   #; because pres_logmov may have changed in the porp_avoid_land procedure
    right pres_angle                 	#; angle to turn -- pres_angle -- is changed in porp_avoid_land
    set angle_turned_right (angle_turned_right + pres_angle)
    if (angle_turned_right > 180) [ set angle_turned_right angle_turned_right - 360 ]
    if (angle_turned_right < -180) [ set angle_turned_right angle_turned_right + 360 ]
    set pres_angle 0
    set count_i count_i + 1
    if (count_i = 100) [ 
      set enough_water_ahead true 
      if ( debug = 1 ) [ 
        print "caught water-ahead loop"
      ]
    ]
  ]
  #; test depth again, avoid_beh = 5:
  porp_check_depth
  if (not enough_water_ahead) [
    let prev_heading heading
    let p max-one-of neighbors [ bathymetry ]
    face p
    set angle_turned_right (angle_turned_right + pres_angle)
    if (angle_turned_right > 180) [ set angle_turned_right angle_turned_right - 360 ]
    if (angle_turned_right < -180) [ set angle_turned_right angle_turned_right + 360 ]
    set pres_mov 1                                        #; move 100 m towards deeper patch
    if ( debug = 1 ) [ 
      let tmp5 list "beh =  5 ; tck " my_tick
      write tmp5 
    ]
  ]
  #;
  #; slow down if turning sharply:
  if ( pres_mov > 10 and (abs angle_turned_right) > 90 ) [ set pres_mov  pres_mov / 5  ] 
  if ( pres_mov > 7 and (abs angle_turned_right) > 50 ) [ set pres_mov  pres_mov / 2  ] 
  #;
  #; Change direction if attracted / repelled by certain areas (model >= 2)
  let total_dx 0
  let total_dy 0
  if ( not use_exp_food_val ) [
    set total_dx (dx * pres_mov) + (item 0 vt)      #; vt isn't used in porp_std_move till here
    set total_dy (dy * pres_mov) + (item 1 vt)      #; note that dx is change in x if taking ONE step forward
    facexy (xcor + total_dx) (ycor + total_dy)
  ]
  if ( use_exp_food_val ) [
    set CRW_contrib inertia_const + pres_mov * VE_total       #; length of vector pointing in direction predicted by CRW
    set MR_contrib sqrt ( (item 0 vt) * (item 0 vt) + (item 1 vt) * (item 1 vt) )     #; length of vector pointing in direction of remembered food
    set total_dx (dx * CRW_contrib) + (item 0 vt)
    set total_dy (dy * CRW_contrib) + (item 1 vt)
    facexy (xcor + total_dx) (ycor + total_dy)            #; really not needed, it already points that way
  ]
  #; Store turn for calc of turning angle in next step:
  #; let total_turn heading - pres_heading   ; total change in heading, including all adjustments till here. 'pres_heading' was calc in beginning of porp_std_move
  let total_turn subtract-headings heading pres_heading   #; total change in heading, including all adjustments till here. 'pres_heading' was calc in beginning of porp_std_move
  #;
  #; Move: 
  fd pres_mov  #; movement length isn't affected by presence of food
  #;
  if ( debug = 2 ) [ 
    if ( my_tick = 0 ) [ print "dist angle_before_avoid_land angle_turned_right x y" ]
    let tmp-var2 (round ( (10 ^ prev_logmov) * 100) / 100)    ; THIS IS IMPORTANT -- the porp turns before it moves, so turning angle is affected by previous moving dist
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 angle_before_avoid_land
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 angle_turned_right
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 ( xcor * 100 + xllcorner )
    set tmp-var2 word tmp-var2 " "
    set tmp-var2 word tmp-var2 ( ycor * 100 + yllcorner )
    print tmp-var2
  ]
  if ( debug = 7 ) [ 
    print word "CRW_contrib: " list ( (dx * (inertia_const + pres_mov * VE_total)) ) ( (dy * (inertia_const + pres_mov * VE_total)) )
    print word "MR_contrib: " vt
    print word "dx, dy (after): " list (total_dx) (total_dy)
    print word "heading (after): " heading
    print word "total_turn: " heading
  ]   
  #;
  #; Remember current moves for the next iteration
  #; if attraction to food alters the movement angle (i.e. vt != 0), this isn't remembered for next step
  #; set prev_angle angle_turned_right  ; so the additional turn due to attraction to food does not influence turning angle in next step 
  set prev_angle total_turn  #; so the additional turn due to attraction to food DOES influence turning angle in next step 
  set prev_logmov log pres_mov 10  #; total steplength, resulting from vt + pres_mov
  #;
  #; test depth one last time, avoid_beh = 6 - move back on same track:
  if (not ([ bathymetry ] of patch-here > 0) ) [
    let prev_heading heading
    facexy (item 0 item 1 pos_list) (item 1 item 1 pos_list)
    set angle_turned_right (angle_turned_right + pres_angle)
    if (angle_turned_right > 180) [ set angle_turned_right angle_turned_right - 360 ]
    if (angle_turned_right < -180) [ set angle_turned_right angle_turned_right + 360 ]
    setxy (item 0 item 1 pos_list) (item 1 item 1 pos_list)          #; move 100 m towards deeper patch
    if ( debug = 1 ) [ 
      #; print "; -> "
      let tmp6 list "beh =  6 � ; tck " my_tick
      set tmp6 word tmp6 " ; "
      set tmp6 word tmp6 angle_turned_right
      print word tmp6 " degr."
    ]
  ]
  #; update position list:
  let pres_pos list xcor ycor
  set pos_list fput pres_pos pos_list
  if ( length pos_list > memory_max ) [ set pos_list remove-item memory_max pos_list ]   
end 
#;  end porp_std_move


#endregion