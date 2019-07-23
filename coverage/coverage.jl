using Pkg
Pkg.instantiate()

# Only run coverage from linux 1.0 build on travis.
get(ENV, "TRAVIS_OS_NAME", "")       == "linux"   || exit()
get(ENV, "TRAVIS_JULIA_VERSION", "") == "1.0" || exit()

using Coverage

cd(joinpath(dirname(@__FILE__), "..")) do
    Coveralls.submit(process_folder())
    Codecov.submit(process_folder())
end
