using HDF5
using Glob 

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
