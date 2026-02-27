using Plots, Statistics, Printf, LinearAlgebra, Test, LaTeXStrings
using Base
######################################################################################
# single line function
######################################################################################
getSum(x, y) = x + y
x, y = 1, 2
@printf("%d + %d = %d\n", x, y, getSum(x, y))

######################################################################################
# multiple lines function
######################################################################################
function canIVote(age=16)
    if age > 18
        println("Can Vote")
    else
        println("You can't vote")
    end
end

canIVote()
######################################################################################
# in Julia, arguments are always pass by value, i.e. they cannot be changed by local scope
######################################################################################
v1 = 5
function changeV1(v1)
    v1 = 10
end
changeV1(v1)
v1

# use global to change the value of the argument
# a exclamation point at the end of a function name is a convention to tell the programmer that
# what's in the function is going to modify one of the arguments. It doesn't actually do anything special
function changeV1!(v1)
    global v1 = 10
end
changeV1!(v1)
v1

######################################################################################
# anonymous functions. Like arrow functions in javascript, inline and has no name
######################################################################################
# the function part is "x -> x^2"
v2 = map(x -> x^2, [1,2,3])
# add each corresponding element to each other
v3 = map((x, y) -> x + y, [1,2], [3,4])
# add all the numbers from 1 to 100 together
v4 = reduce(+, 1:100)
# find the longest word in a sentence
sentence = "This is a string with a short word and a long word but who knows cuz that's crazy superduperduper crazy I mean what's going on"
sArray = split(sentence)
longest = reduce((x, y) -> length(x) > length(y) ? x : y, sArray)

######################################################################################
# math formulas without unnecessary multiplication symbols
######################################################################################
x = 2
y = -1
2x^2+5y^2+2x*y+cos(5x)y
(1e-7x^2)y

2x^2

######################################################################################
# Broadcasting
######################################################################################
"""
https://julialang.org/blog/2017/01/moredots/
"""

# the `.` is for broadcasting the operation to every value in the array
# unlike python where x * y automatically performs multiplication element-wise
# multiply every value in an array by 3
# multiplication by a scalar works element-wise automatically
A = [1,2,3]
A * 3 == A .* 3  
# but for more complicated operations, it doesn't work 
A = A^2
# but this works
A = A.^2

# if you want to apply cos to each element in the array, need the `cos.()`
a = ones(3,3)
cos(a)
cos.(a)

B = fill(1.0, (2,2))
sin(B)
sin.(B)

# this is the proper exponential of a matrix
exp(B)
# this does e^x for each element x
exp.(B)

######################################################################################
# Symbols and dictionaries
######################################################################################
# symbols have `:` before them
:Matthew
println(:Matthew)

# a dictionary can use symbols
d2 = Dict(:pi=>3.14, :e=>2.718)



######################################################################################
# Subtypes and AbstractTypes
######################################################################################
"""
from https://www.youtube.com/watch?v=4giNd6HLUQg
<: means "is a subtype of"
abstract types cannot be instantiated, and it looks like they are used just as 
orginization for a collection of different types
Ex. The Ising struct is a child of the Hamiltonian type. We never instantiate a Hamiltonian type,
we just instantiate the Ising type
"""
abstract type Hamiltonian{T, S<:AbstractMatrix{<:T}} end

LinearAlgebra.ishermitian(H::Hamiltonian) = true
LinearAlgebra.eigvals(H::Hamiltonian) = eigvals(parent(H))

struct Ising{T, S<:AbstractMatrix{<:T}} <: Hamiltonian{T, S}
    H::Hermitian{T, S}
    J::Float64
    h::Float64
    L::Int
end

Base.parent(IH::Ising) = IH.H

my_h = Ising(Hermitian([1. 0.; 0. -1.]), 1., 0., 1)
println(ishermitian(my_h))
println(eigvals(my_h))

@testset for uplo in (:U, :L)
    T = Bidiagonal(dv, ev, uplo)
    @testset "Constructor and basic properties" begin
        @test size(T, 1) == size(T, 2) == n
        @test size(T) == (n, n)
        @test Array(T) == diagm(0 => dv, (uplo == :U ? 1 : -1) => ev)
        @test Bidiagonal(Array(T), uplo) == T
        @test big.(T) == T
    end
end

######################################################################################
# Python f string equivalent: 
######################################################################################
name = 4
string = "Hello, $(name)!"

######################################################################################
# Array, Matrix, and Vector explanation
######################################################################################
"""
T: data type, N: dimension,
   `Array{T,N}` type is the most-general representation in Julia
   `Vector{T}` is a special case of Array{T,N}, where N=1, i.e. one-dimensinoal array
   `Matrix{T}` is a special case of Array{T,N}, where N=2, i.e. two-dimensional array 

    Anything higher than N=2 will be the Array{T,N} type automatically

Now it's 2026, and I originally wrote this guide in 2022-2023. It looks like I got on
the wrong path trying to use the Array type to do matrix math. I don't need to do that, I can 
just use Matrix and Vector, and I definitely should not be creating arrays that have 0 dimensions.
This was probably when I was just learning numpy and I was trying to add (3,1) vectors one at a time
to create a (N,3) matrix. Now in 2026, I know the better way to do this is
    A = zeros(N,3)
    for i=1:N
        A[i,:] = [a, b, c]
    end

Also, `nothing` might be bugged after version updates since 2022 
"""
zeroArray = zeros(40)
# create an empty array of type Int64
newArray = Int64[]
# appending to the empty array just adds one value
append!(newArray, 3)
# can't assign a Float64 to this array
append!(newArray, 1.2)

# assigning the zeroArray to emptyArray makes it the same size
# Unlike in MATLAB, don't have to define a zero array with a explicit size
# and then have to make sure the array you want to add is the same size
emptyArray = zeroArray
emptyArray


# Declaring arrays with pre-determined size

# declare a (3,2) empty array of Float64
a = Array{Float64}(undef, 3, 2)
# OR
a = Array{Float64}(nothing, 3, 2)

# declare an empty row vector with 3 components
a = Array{Float64}(undef, 0, 3)

# or just use Vector type
# EQUIVALENT to Array{T,1}
a = Vector{Float64}(undef, 3)
b = Vector{Float64}(undef, 7)

# declare an empty column vector with 3 components
a = Array{Float64}(undef, 3, 0)

# vertically stack a (10,3) array to an empty array
# need to make sure the 2nd dimension is equal for vcat (1st dimension equal for hcat)
# setting the first dimension to 0 means that when you 
# stack the array there will not be an extra row
a = Array{Int32}(undef, 0, 3)
b = fill(5, (10, 3))
b
c = vcat(a, b)
c

######################################################################################
# list comprehensions
######################################################################################
# the number of index variables in the comprehension corresponds to the number of dimensions in the array
simpleComp = [i for i=1:10]
complexComp = [4*x^2 - 3.5*x + 1.3 for x in range(1,10)]
complexComp2d = [4*y+3*x for x=1:10, y=1:10]
complexComp3d = [x+y+z for x=1:10, y=1:10, z=1:10]
conditionalComp = [i+j for i=1:10, j=1:10 if i - j == 1]

######################################################################################
# Linear Algebra Stuff
######################################################################################
"""
`using LinearAlgebra` is required for most functions, need to install
https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/
"""

A = [1, 2, 3]     # this is a vector
B = [1  2  3]     # this is a matrix (no commas)
C = [1 2 3; 4 5 6] # this is a matrix

using LinearAlgebra
nDim = 6
upperDiagonal = fill(-1, nDim - 1)
diagonal = rand(Float64, nDim)
# special bidiagonal type. A matrix with a diagonal and then upper diagonal and/or lower diagonal
A = Bidiagonal(diagonal, upperDiagonal, :U)
typeof(A)
A = Symmetric(A, :U)
eigvals(A)
eigvecs(A)

firstEigenvector = eigvecs(A)[:,1]

# get the valid indices of the array
axes(A)
axes(A,1)

# TODO BUGGED 2/17/26 sort along the second dimension from smallest to largest
sort!(A, dims=2)

# no need for a special matrix multiplication symbol or function like `A@B` in Python 
BA = kron(ones(nDim), rand(Float64, nDim)')
CA = rand(Float64, nDim)
BA*CA

# Create a sparse diagonal matrix using vectors as the diagonals
using SparseArrays
A = spdiagm(
    -1 => fill(-1.0, 99),
    0  => fill(2.0, 100),
    1  => fill(-1.0, 99)
)

# get the eigenvalues of the sparse matrix using the `partialschur` decomposition from ArnoldiMethod.jl
# ArnoldiMethod.jl is Julia's version of the ARPACK library (written in Fortran, stands for ARnoldiPACKage) 
# which calculates large-scale eigenvalue problems
using ArnoldiMethod
decomp, history = partialschur(A, nev=10, tol=1e-6);
decomp
history

# convert something into a vector, i.e. collapse something to a 1D array
E = [1 2 3; 4 5 6]
vec(E)
# the same as this, but here you need to know the size of the reshaped array
reshape(E, (6))


# tensor products
⊗(x,y) = kron(x,y)
A = [1 1; 2 2]
B = [3 3; 4 4]
A ⊗ B # (4,4)

A = [1, 2, 3]
B = [4, 5, 6]
A ⊗ B # (9)

# trace
tr(A)




######################################################################################
# Julia's version of np.linalg.norm(A, axis=1)
######################################################################################
"""
Get the magnitude of a matrix with dimensions like (nPoints,3)
(bx1,by1,bz1)
(bx2,by2,bz2)
(bx3,by3,bz3)
...
(nPoints,3) -> (nPoints,1)
"""
bfield = rand(Float16, (10,3))
bfield
mag = sqrt.(sum(abs2,bfield,dims=2)) # stays as Matrix type with one dimension=1
mag = vec(sqrt.(sum(abs2,bfield,dims=2))) # Vector type
bfield_norm = bfield ./ mag
round.(sqrt.(sum(abs2,bfield_norm,dims=2)), digits=3) # this should be equal to 1

# another way with LinearAlgebra.jl
using LinearAlgebra
[norm(r) for r in eachrow(bfield)]



######################################################################################
# Create meshgrids like np.meshgrid(np.linspace(-10,10,10), np.linspace(-5,5,10))
######################################################################################

# 2D
Nx = 10
Ny = 10
xs = range(-10,10,Nx)
ys = range(-10,10, Ny)
X = [x for x in xs, y in ys]   # x varies along columns (1st dimension)
Y = [y for x in xs, y in ys]   # y varies along rows (2nd dimension)

# 3D
N = 10
L = 1
X = [x for x in range(-L/2,L/2,N), y in range(-L/2,L/2,N), z in range(-L/2,L/2,N)]
Y = [y for x in range(-L/2,L/2,N), y in range(-L/2,L/2,N), z in range(-L/2,L/2,N)]
Z = [z for x in range(-L/2,L/2,N), y in range(-L/2,L/2,N), z in range(-L/2,L/2,N)]

# here's some function I found a while ago. Looks like it creates a vector of all the points, instead of Matrix or 3D Array
meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
x,y = meshgrid(range(-5,5,11), range(-5,5,11))
x
y





######################################################################################
# PLOTTING LIBRARIES
######################################################################################
"""
It seems that when you are running a session in the REPL with multiple
plotting packages included, like Plots.jl and Makie.jl, there are naming clashes
when doing `plot`, `scatter`, `contour`, etc. So you need to put the package's name in front 
of the function like 
`Plots.scatter!()` and `Makie.scatter!()`


###############################################
Plots.jl
###############################################
Used for:
    - quiver(x,y,quiver=(u,v))        2D
    - quiver(x,y,z,quiver=(u,v,w))    3D
    - plot()
    - scatter()
    - contour()
    - contourf()

Plot options (ones that are passed in the function name)
    https://docs.juliaplots.org/latest/generated/attributes_subplot/

    xscale=:log10
    yscale=:log10
    grid=true
    colobar_title="Title"
    colorbar_entry=true (enable or disable colorbar)


  Abbreviated options:
    linecolor: lc
    linewidth: lw
    colormap: c

Colormaps:
    https://docs.juliaplots.org/dev/generated/colorschemes/
    :viridis
    :thermal
    :turbo

Plotting Backends:
    Installed with plots.jl:
        - gr(), default
        - plotly()
    Need to install before using them:
        - pythonplot() aka matplotlib
        - gaston()


Plots using gr() backend:
    - quiver(x,y,u,v)



###############################################
Makie.jl
###############################################

This is a stand-alone plotting library. It looks like it is built on 
proper graphics pipelines so it is good for 3D plotting

Backends:
    GLMakie: uses OpenGL, runs on GPU, best for 3D plots
    CairoMakie: uses Cairo, only on CPU, best for 2D plots and vector graphics
    WGLMakie: uses WebGL

Use this for:
    - volume()
    - contour3d()
    - density()

Not supported:
    - streamplot() -> bugged on 2/13/26
"""

# a simple graph
using CairoMakie
seconds = 0:0.1:2
measurements = [8.2, 8.4, 6.3, 9.5, 9.1, 10.5, 8.6, 8.2, 10.5, 8.5, 7.2,
        8.8, 9.7, 10.8, 12.5, 11.6, 12.1, 12.1, 15.1, 14.7, 13.1]

f = Figure()
ax = Axis(f[1, 1],
    title = "Experimental data and exponential fit",
    xlabel = "Time (seconds)",
    ylabel = "Value",
)
CairoMakie.scatter!(
    ax,
    seconds,
    measurements,
    color = :tomato,
    label = "Measurements"
)
CairoMakie.lines!(
    ax,
    seconds,
    exp.(seconds) .+ 7,
    color = :tomato,
    linestyle = :dash,
    label = "f(x) = exp(x) + 7",
)
CairoMakie.axislegend(position = :rb)
f
save("first_figure.png", f)


######################################################################################
# Create a 3D vector plot with Plots.jl
######################################################################################
# looks like you don't need to create 3D axes, it just works when you pass in (x,y,z, (u,v,w))
using Plots
ϕs = range(-π, π, length = 50)
θs = range(0, π, length = 25)
θqs = range(1, π - 1, length = 25)
x = vec([sin(θ) * cos(ϕ) for (ϕ, θ) = Iterators.product(ϕs, θs)])
y = vec([sin(θ) * sin(ϕ) for (ϕ, θ) = Iterators.product(ϕs, θs)])
z = vec([cos(θ) for (ϕ, θ) = Iterators.product(ϕs, θs)])
u = 0.1 * vec([sin(θ) * cos(ϕ) for (ϕ, θ) = Iterators.product(ϕs, θqs)])
v = 0.1 * vec([sin(θ) * sin(ϕ) for (ϕ, θ) = Iterators.product(ϕs, θqs)])
w = 0.1 * vec([cos(θ) for (ϕ, θ) = Iterators.product(ϕs, θqs)])
quiver(x, y, z, quiver = (u, v, w))


######################################################################################
# How to create 3D axes using Makie
######################################################################################
using Makie
using GLMakie
fig = Figure()
a1 = Axis3(fig[1, 1], aspect = (1, 1, 1), title = "aspect = (1, 1, 1)")
a2 = Axis3(fig[1, 2], aspect = (2, 1, 1), title = "aspect = (2, 1, 1)")
a3 = Axis3(fig[2, 1], aspect = (1, 2, 1), title = "aspect = (1, 2, 1)")
a4 = Axis3(fig[2, 2], aspect = (1, 1, 2), title = "aspect = (1, 1, 2)")
fig
N = 100
GLMakie.scatter!(a1, rand(N), rand(N), rand(N), color=:red)
s = range(-2π,2π,N)
z = range(-1, 1, N)
GLMakie.lines!(a1, sin.(s), cos.(s), color=:red)
GLMakie.lines!(a2, sin.(s), cos.(s), color=:red)
GLMakie.lines!(a1, sin.(s), cos.(s), z, color=:green)
GLMakie.lines!(a2, sin.(s), cos.(s), z, color=:green)
fig


######################################################################################
# GLMakie contour() vs contour3d() plot
######################################################################################
"""
Boths of these plots are plotted in 3D axes
However, 
    contour() in 3D axes plots the isosurface, i.e. constant values of a 3D array
    contour(1D, 1D, 3D)
    https://docs.makie.org/stable/reference/plots/contour

    contour3d() plots contour levels as the z height in 3D, so same as a normal contour plot in 2D
    contour3d(1D, 1D, 2D)
    https://docs.makie.org/stable/reference/plots/contour3d
"""

using GLMakie
r = range(-pi, pi, length = 21)
data2d = [cos(x) + cos(y) for x in r, y in r]
data3d = [cos(x) + cos(y) + cos(z) for x in r, y in r, z in r]

f = Figure(size = (700, 400))
a1 = Axis3(f[1, 1], title = "3D contour()")
GLMakie.contour!(a1, -pi .. pi, -pi .. pi, -pi .. pi, data3d)

a2 = Axis3(f[1, 2], title = "contour3d()")
GLMakie.contour3d!(a2, r, r, data2d, linewidth = 3, levels = 10)
f


######################################################################################
# Plot 2D contour slices in 3D using `volumeslices!()`
# this one is interactive and allows you to move the planes along their axes
######################################################################################
using GLMakie
fig = Figure()
ax = LScene(fig[1, 1], show_axis=false)

x = LinRange(0, π, 50)
y = LinRange(0, 2π, 100)
z = LinRange(0, 3π, 150)

sgrid = SliderGrid(
    fig[2, 1],
    (label = "yz plane - x axis", range = 1:length(x)),
    (label = "xz plane - y axis", range = 1:length(y)),
    (label = "xy plane - z axis", range = 1:length(z)),
)

lo = sgrid.layout
nc = ncols(lo)

vol = [cos(X)*sin(Y)*sin(Z) for X ∈ x, Y ∈ y, Z ∈ z]
plt = volumeslices!(ax, x, y, z, vol)

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

fig