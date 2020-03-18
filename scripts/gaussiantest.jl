using PyCall, Conda, ImageFiltering
using Plots

# Pycall setup:
# Conda.add("scipy")

py"""
import scipy.ndimage.filters as filt
def pyfilter(x):
    return filt.gaussian_filter1d(x,40)
"""

pyfilter(x) = py"pyfilter"(x)

data = rand(10)

plot(data; label = "Original data")
plot!(pyfilter(data); label = "Python filtering")
# "
plot!(imfilter(data, Kernel.gaussian((40,))); label = "Julia filtering")
