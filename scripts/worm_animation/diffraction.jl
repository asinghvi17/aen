image = AbstractPlotting.colorbuffer(display(sc))

using FFTW, Colors

gimg = Gray.(image)

transformed_image = fft(Float64.(gimg))[100:end-100, 100:end-100]

graytrans = Gray.(abs2.(transformed_image))

using DSP

pg = DSP.Periodograms.Periodogram2(Float64.(gimg), axes(gimg)...)
