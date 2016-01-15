
# This is just a simple link to the real script that is installed in the julia package.
# By keeping this file as just a link, we can maintain the script using the Julia
# package manager (i.e. Pkg.update())
include("$(Pkg.dir())/ChromNet/scripts/build_bundle.jl")
