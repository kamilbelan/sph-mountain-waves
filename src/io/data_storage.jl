using HDF5
using Glob 
using Dates
using DrWatson
using SmoothedParticles

# ==============
# Directory and Metadata Setup
# ==============

"""
    initialize_run_directory(sim_params::Dict)::String

Creates (and outputs) a separate directory based on the simulation parameters, together with the DrWatson metadata.
"""

function initialize_run_directory(sim_params::Dict)::String
	# choose the parameters for the folder name
	name_keys = [:dr, :dt_rel, :t_end, :gamma_r_rel]
	name_params = filter(p -> p.first in name_keys, sim_params)
	short_name = savename(name_params)

	# extract model name safely 
	model_name_str = string(sim_params[:model])

	# generate a hex hash for uniqueness
	full_hash = string(hash(sim_params), base=16)[1:6]
	run_name = "$(short_name)_$(full_hash)"

	# build the directory tree: data/sims/modelname/YYYY-MM-DD/run_name
	date_str = Dates.format(now(), "yyyy-mm-dd")
	run_dir = datadir("sims", model_name_str, date_str, run_name)
	mkpath(run_dir)

	println("Output directory: $run_dir")

	# save metadata in the run directory
	metadata_dict = Dict(String(k) => v for (k, v) in sim_params)
	metadata_dict["module"] = model_name_str 
	@tagsave(joinpath(run_dir, "metadata.jld2"), metadata_dict)
	return run_dir
end

# ==============
# Writing data to frames
# ==============

"""
    write_file(run_dir::String, sys::ParticleSystem, frame_counter::Int64, t::Float64)

Creates a H5 frame file and writes the simulation's frame's data to it.
"""

function write_frame!(run_dir::String, sys::ParticleSystem, frame_counter::Int, t::Float64)
	frame_file = joinpath(run_dir, "frame_$(lpad(frame_counter, 4, '0')).h5")
	# pre-allocate flat arrays
	N_parts = length(sys.particles)
	pos_matrix = zeros(Float64, 3, N_parts)
	vel_matrix = zeros(Float64, 3, N_parts)
	densities = zeros(Float64, N_parts)
	densities_pert = zeros(Float64, N_parts)
	pressures = zeros(Float64, N_parts)
	pressures_pert = zeros(Float64, N_parts)
	temperatures = zeros(Float64, N_parts)
	pot_temperatures = zeros(Float64, N_parts)
	types = zeros(Float64, N_parts)

	@inbounds for (i, p) in enumerate(sys.particles)
		pos_matrix[1, i] = p.x[1]
		pos_matrix[2, i] = p.x[2]
		pos_matrix[3, i] = 0.0  

		vel_matrix[1, i] = p.v[1]
		vel_matrix[2, i] = p.v[2]
		vel_matrix[3, i] = 0.0  

		densities[i] = p.ρ
		densities_pert[i] = p.ρ′
		pressures[i] = p.P
		pressures_pert[i] = p.P′
		temperatures[i] = p.T
		pot_temperatures[i] = p.θ
		types[i] = p.type
	end

	#  write securely to HDF5
	h5open(frame_file, "w") do file
		# save metadata as HDF5 attributes 
		attributes(file)["time"] = t
		attributes(file)["frame_counter"] = frame_counter
		attributes(file)["n_particles"] = N_parts

		# save the arrays as primary datasets
		file["positions"] = pos_matrix
		file["velocities"] = vel_matrix
		file["densities"] = densities
		file["densities_pert"] = densities_pert
		file["pressures"] = pressures
		file["pressures_pert"] = pressures_pert
		file["temperatures"] = temperatures
		file["pot_temperatures"] = pot_temperatures
		file["types"] = types
	end
end

# ==============
# XDMF generation
# ==============

""" 
    generate_sph_xdmf(run_dir::String)

Finds all HDF5 frames in `run_dir and creates a XDMF file for Paraview visualisation.
"""

function generate_sph_xdmf(run_dir::String)
	# find all the HDF5 frames in the given directory
	h5_files = sort(glob("frame_*.h5", run_dir))

	if isempty(h5_files)
		println("No .h5 files found in $run_dir!")
		return
	end

	xdmf_path = joinpath(run_dir, "output.xdmf")

	open(xdmf_path, "w") do io
		# write the XDMF header 
		write(io, """
	<?xml version="1.0" ?>
	<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
	<Xdmf Version="3.0">
	<Domain>
	<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
	""")

		# loop through every frame 
		for file in h5_files
			filename = basename(file)

			# dead the metadata saved as attributes
			h5open(file, "r") do h5
				t = read(attributes(h5)["time"])
				N = read(attributes(h5)["n_particles"])

				write(io, """
	  <Grid Name="Particles" GridType="Uniform">
	  <Time Value="$t" />
	  <Topology TopologyType="Polyvertex" NumberOfElements="$N"/>

	  <Geometry GeometryType="XYZ">
	  <DataItem Dimensions="$N 3" NumberType="Float" Precision="8" Format="HDF">
	  $filename:/positions
	  </DataItem>
	  </Geometry>

	  <Attribute Name="Velocity" AttributeType="Vector" Center="Node">
	  <DataItem Dimensions="$N 3" NumberType="Float" Precision="8" Format="HDF">
	  $filename:/velocities
	  </DataItem>
	  </Attribute>

	  <Attribute Name="Density" AttributeType="Scalar" Center="Node">
	  <DataItem Dimensions="$N" NumberType="Float" Precision="8" Format="HDF">
	  $filename:/densities
	  </DataItem>
	  </Attribute>

	  <Attribute Name="Pressure" AttributeType="Scalar" Center="Node">
	  <DataItem Dimensions="$N" NumberType="Float" Precision="8" Format="HDF">
	  $filename:/pressures
	  </DataItem>
	  </Attribute>

	  <Attribute Name="Particle_Type" AttributeType="Scalar" Center="Node">
	  <DataItem Dimensions="$N" NumberType="Float" Precision="8" Format="HDF">
	  $filename:/types
	  </DataItem>
	  </Attribute>
	  </Grid>
	  """)
			end
		end

		# Close the XML tags
		write(io, """
	</Grid>
	</Domain>
	</Xdmf>
	""")
	end

	println("successfully generated $(xdmf_path)")
end

