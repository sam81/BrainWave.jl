using BrainWave
#using PlotlyJS

sf = 512
dur = 1.2
nSamp = round(Int, sf*dur)
tArr = collect(0:nSamp-1)./sf
dat = zeros(2, nSamp)
dat[1,:] = sin.(2*pi*75 .* tArr)
dat[2,:] = sin.(2*pi*120 .* tArr)

p1, f1 = getSpectrum(dat[1,:], sf, powerOfTwo=false)
p1b, f1b = getSpectrum(dat[1,:], sf, powerOfTwo=true)

p2, f2 = getSpectrum(dat[2,:], sf, powerOfTwo=false)
p2b, f2b = getSpectrum(dat[2,:], sf, powerOfTwo=true)

# s1 = scatter(;x=f1, y=10*log10.(p1))
# s1b = scatter(;x=f1b, y=10*log10.(p1b))

# s2 = scatter(;x=f2, y=10*log10.(p2))
# s2b = scatter(;x=f2b, y=10*log10.(p2b))
