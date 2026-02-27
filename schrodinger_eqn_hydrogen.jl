using LinearAlgebra, SparseArrays, Kronecker
using Arpack
using ArnoldiMethod, LinearAlgebra, SparseArrays # for partialschur
using Plots

# Hydrogen Atom in 3D
# Set up discrete grid of sample points
# units of meters, so we are calculating on a mesh that spans 50 Bohr radii in X,Y,Z
N = 100
a0 = 5.29e-11
L = 50*a0
# L = 50
X = [x for x in range(-L/2,L/2,N), y in range(-L/2,L/2,N), z in range(-L/2,L/2,N)]
Y = [y for x in range(-L/2,L/2,N), y in range(-L/2,L/2,N), z in range(-L/2,L/2,N)]
Z = [z for x in range(-L/2,L/2,N), y in range(-L/2,L/2,N), z in range(-L/2,L/2,N)]
dx = step(range(-L/2,L/2,N))


function get_hydrogen_potential(x::Array{Float64, 3}, y::Array{Float64, 3}, z::Array{Float64, 3})
    return -dx^2 ./ (sqrt.(x.^2 + y.^2 + z.^2 .+ 1e-10))
end

# construct the Hamiltonian
V = get_hydrogen_potential(X, Y, Z)
maximum(V)

reshape(V, (N^3))
# the finite difference approximation of the laplacian is -2 on the diagonal, and +1 on the next upper and lower diagonals
# this is called a tridiagonal matrix 
D = spdiagm(
    -1 => fill(1, N-1),
    0 => fill(-2, N),
    1 => fill(1, N-1)
)

U = spdiagm(
    0 => reshape(V, (N^3))
)

id = Matrix(1I, N, N)
ħ = 6.626e-34 / (2π)
m = 9.11e-31
k = -1/2 * (ħ^2/(2*m*dx^2))
# k = -1/2

# This is called the kronecker sum which is defined as A ⊕ B = A ⊗ I + I ⊗ B
T = k .* kron(kron(D, id) .+ kron(D, id), D)

H = T + U

# use the eigenvalue solver: partial schur decomposition from ArnoldiMethod.jl
nEigen = 10
@time decomp, history = partialschur(H, nev=nEigen)

# convert to a Matrix type instead of partialschur function which returns a View into an array
evecs = decomp.Q |> Matrix


groundstate = evecs[:,1]

# # try a different solver
# @time evals, evevs = Arpack.eigs(H, nev=6)
# evecs = evevs


# the max value of the probability values tells you where to draw contour lines
maxvalue = maximum(abs2.(groundstate))

function plot_prob(energy_level, xIndex, yIndex, zIndex)
    axesValues = range(-L/2,L/2,N)

    l = @layout [a b c]

    state = evecs[:,energy_level]
    state = reshape(state, (N,N,N))
    # NORMALIZE
    # println(maximum(state))
    # state /= maximum(state)
    stateXSlice = abs2.(state[xIndex,:,:])
    println(maximum(stateXSlice))
    stateYSlice = abs2.(state[:,yIndex,:])
    println(maximum(stateYSlice))
    stateZSlice = abs2.(state[:,:,zIndex])
    println(maximum(stateZSlice))


    p1 = Plots.heatmap(axesValues, axesValues, stateXSlice, c=:viridis, xlabel="Y", ylabel="Z")
    p2 = Plots.heatmap(axesValues, axesValues, stateYSlice, c=:viridis, xlabel="X", ylabel="Z")
    p3 = Plots.heatmap(axesValues, axesValues, stateZSlice, c=:viridis, xlabel="X", ylabel="Y")


    sizePixels = 600
    figure = Plots.plot(p1, p2, p3, layout=l, size=(3*sizePixels,sizePixels), colorbar=false, left_margin=10Plots.mm, right_margin=10Plots.mm, bottom_margin=10Plots.mm)

    display(figure)
end

nState = 5
plot_prob(nState, floor(Int, N/2)-1, floor(Int, N/2)+1, floor(Int, N/2)+1)

Plots.plot()
state = evecs[:,1]
state = reshape(state, (N,N,N))
state = abs2.(state[2,:,:])
Plots.heatmap!(state)
# Plots.heatmap!(state ./ maximum(state))

stateXSlice = abs2.(state[x,:,:])



# plot a 3D isosurface plot
using GLMakie
using CairoMakie
using Plots

GLMakie.activate!()

begin
    nState = 1
    f = Figure(size = (700, 400))
    a1 = Axis3(f[1, 1], title = "Hydrogen Atom n=$(nState)")
    xx = range(-L/2,L/2,N)
    state = reshape(evecs[:,nState], N,N,N)
    prob3d = abs2.(state)

    # GLMakie.contour!(-L/2..L/2, -L/2..L/2, -L/2..L/2, prob3d, levels=3)
    GLMakie.volume!(prob3d, algorithm=:iso, alpha=1, isovalue=maxvalue/4, isorange=maxvalue/2)
end
