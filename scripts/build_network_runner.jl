
# make sure the ChromNet package is installed
try
    @assert Pkg.installed("ChromNet") != nothing
catch
    Pkg.clone("https://github.com/slundberg/ChromNet.jl.git")
end
try
    @assert Pkg.installed("SamIO") != nothing
catch
    Pkg.clone("https://github.com/slundberg/SamIO.jl.git")
end

# This is just a simple link to the real script that is installed in the julia package.
# By keeping this file as just a link, we can maintain the script using the Julia
# package manager (i.e. Pkg.update())
include("$(Pkg.dir())/ChromNet/scripts/build_network.jl")
