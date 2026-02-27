function biot_savart(coilVertices, targetPoints, current)
    """
    Calculate the biot-savart law at `targetPoints` due to current carying segments defined by `coilVertices`

    The `coilVertices` do NOT need to be connected. While realistic currents have to be connected (i.e. form a loop)
    this function allows for an unconnected path for current. 

    For example, allows you to calculate the magnetic field due to an infinite, stright wire 

    Args:
    coilVertices: 
        The vertices of the coil segments, size (nVertices, 3)
        Should NOT be the displacement between the vertices i.e. dl
    targetPoints:
        The points where the magnetic field is calculated, size (nTargetPoints, 3)
    current:
        The magnitude of the current in Amps

    -------------------
    Example calculation
    -------------------

    coilVertices (4 total):
    1      2      3      4

    dl           (always 1 less than number of coil vertices):
    1      2      3
    directions:
    1->2   2->3   3->4

    targetPoints (does not need to be the same size as coilVertices)
    Example 2d grid:
    1  2  3  
    4  5  6
    7  8  9

    Returned magnetic field array: (nTargetPoints, 3)
    Bx1 By1 Bz1
    Bx2 By2 Bz2
    ...
    Bx9 By9 Bz9
    """
    prefactor = (4π*1e-7/4π)*current
    vectorField = zeros(size(targetPoints))
    nCoilVertices = size(coilVertices, 1)
    nTargetPoints = size(targetPoints, 1)
    # calculate the field for each target point iteratively
    # biot-savart law says every current dl vector contributes to the b-field
    for i=1:nTargetPoints
        thisPoint = targetPoints[i,:]
        # rtilded is defined as: evaluation point - source point
        # since we have access to all source points in `coilVertices` which is a (nPoints,3) matrix,
        # we can create a (nPoints,3) matrix of the evaluation point which is just the point duplicated at every row
        # faster matrix math than looping through each source point and calculating rtilded
        # BUT the size of the evaluationPoint matrix needs to be (nPoints-1,3) because dl vectors have 1 less entry than the number of coilVertices
        thisPointVectorized = ones(nCoilVertices-1, 1) * transpose(thisPoint)
    
        rtilded_vector = thisPointVectorized .- coilVertices[1:end-1,:]
        
        rtilded_magnitude = sum(abs2, rtilded_vector, dims=2).^(1/2)
        dl = diff(coilVertices, dims=1)
        dlx = dl[:,1]
        dly = dl[:,2]
        dlz = dl[:,3]
        rx = rtilded_vector[:,1]
        ry = rtilded_vector[:,2]
        rz = rtilded_vector[:,3]
        bx = @.    (dly*rz-dlz*ry) / rtilded_magnitude^3
        by = @. -1*(dlx*rz-dlz*rx) / rtilded_magnitude^3
        bz = @.    (dlx*ry-dly*rx) / rtilded_magnitude^3
        vectorField[i,1] = sum(bx)
        vectorField[i,2] = sum(by)
        vectorField[i,3] = sum(bz)
    end

    return prefactor .* vectorField
end


function circular_loop_xy(radius, nPoints, zValue)
    """
    Create a circular loop in the xy plane at a specified z height
    """
    thetas = LinRange(0,2π,nPoints)
    xVals = []
    yVals = []
    zVals = []
    for theta in thetas
        x = radius * sin(theta)
        y = radius * cos(theta)
        append!(xVals, x)
        append!(yVals, y)      
    end
    zVals = fill(zValue, nPoints)

    return hcat(xVals, yVals, zVals)

end


# Geometry 1: circular loop of current in XY plane
coilVertices = circular_loop_xy(1.0, 100, 0.0)

# Geometry 2: straight wire in z direction
N = 100
coilVertices = hcat(range(-1,1,N), zeros(N), zeros(N))

# Geometry 3: straight wire in y direction
N = 100
coilVertices = hcat(zeros(N), range(-1,1,N), zeros(N))

# Geometry 4: straight wire in z direction
N = 100
coilVertices = hcat(zeros(N), zeros(N), range(-1,1,N))

#############################################################################################
# 2d vector plot of B-magnitude in the XZ plane
#############################################################################################
L=4
N=100
x2d = [x for x in LinRange(-L/2,L/2,N), z in range(-L/2,L/2,N)]
z2d = [z for x in LinRange(-L/2,L/2,N), z in range(-L/2,L/2,N)]
x1d = reshape(x2d, (N^2,1))
z1d = reshape(z2d, (N^2,1))
y1d = zeros(size(z1d))

targetPoints = hcat(x1d, y1d, z1d)

current = 1e-3
@time bfield = biot_savart(coilVertices, targetPoints, current);
bx1d = bfield[:,1]
by1d = bfield[:,2]
bz1d = bfield[:,3]
bx2d = reshape(bfield[:,1],(N,N))
by2d = reshape(bfield[:,2],(N,N))
bz2d = reshape(bfield[:,3],(N,N))


bMag2d = sqrt.(bx2d.^2 + by2d.^2 + bz2d.^2)
# For some reason, you need to transpose the array after the reshape. Looks like something to do with row-major vs column-major ordering
bMag2d = transpose(bMag2d)

using Plots
pythonplot()
Plots.plot()
Plots.quiver(x1d, z1d, quiver=(bx1d, bz1d))
Plots.plot()
Plots.contourf(
    LinRange(-L/2,L/2,N),
    LinRange(-L/2,L/2,N),
    log10.(bMag2d),
    c=:Reds, 
    colorbar_entry=true,
    colorbar_title="Log_10(B_mag)")
Plots.xlabel!("X")
Plots.ylabel!("Z")

#############################################################################################
# Plot a 3d Volume Slices plot using GLMakie, interactive
#############################################################################################

# target points
L = 4
N = 200

x2d = [x for x in LinRange(-L/2,L/2,N), y in range(-L/2,L/2,N), z in range(-L/2,L/2,N)]
y2d = [y for x in LinRange(-L/2,L/2,N), y in range(-L/2,L/2,N), z in range(-L/2,L/2,N)]
z2d = [z for x in LinRange(-L/2,L/2,N), z in range(-L/2,L/2,N), z in range(-L/2,L/2,N)]

x1d = reshape(x2d, (N^3,1))
y1d = reshape(y2d, (N^3,1))
z1d = reshape(z2d, (N^3,1))

targetPoints = hcat(x1d, y1d, z1d)

current = 1e-3

@time bfield = biot_savart(coilVertices, targetPoints, current);
bfield_mag = sqrt.(sum(abs2, bfield, dims=2))

bfield_mag_3d = reshape(bfield_mag, (N,N,N))
bfield_mag_3d ./= maximum(bfield_mag_3d)
b3d_log = log10.(bfield_mag_3d)

using GLMakie
fig = Figure()
ax = LScene(fig[1, 1], show_axis=false)

x = range(-L/2,L/2,N)
y = range(-L/2,L/2,N)
z = range(-L/2,L/2,N)

sgrid = SliderGrid(
    fig[2, 1],
    (label = "yz plane - x axis", range = 1:length(x)),
    (label = "xz plane - y axis", range = 1:length(y)),
    (label = "xy plane - z axis", range = 1:length(z)),
)

lo = sgrid.layout
nc = ncols(lo)


plt = volumeslices!(ax, x, y, z, b3d_log)

# connect sliders to `volumeslices` update methods
sl_yz, sl_xz, sl_xy = sgrid.sliders

on(sl_yz.value) do v; plt[:update_yz][](v) end
on(sl_xz.value) do v; plt[:update_xz][](v) end
on(sl_xy.value) do v; plt[:update_xy][](v) end

set_close_to!(sl_yz, .5length(x))
set_close_to!(sl_xz, .5length(y))
set_close_to!(sl_xy, .5length(z))

# add toggles to show/hide heatmaps
hmaps = [plt[Symbol(:heatmap_, s)][] for s ∈ (:yz, :xz, :xy)]
toggles = [Toggle(lo[i, nc + 1], active = true) for i ∈ 1:length(hmaps)]

map(zip(hmaps, toggles)) do (h, t)
    on(t.active) do active
        h.visible = active
    end
end

# cam3d!(ax.scene, projectiontype=Makie.Orthographic)