module Extract
export PositiveThreshold, NegativeThreshold, DoubleThreshold, extract_wf, extract_wf_teo, noise_autocovariance

immutable NegativeThreshold{T<:Real}
    thr::T
    NegativeThreshold(thr) = new(-abs(thr))
end
NegativeThreshold(thr::Real) = NegativeThreshold{typeof(thr)}(thr)

function extract_wf(x, thr_type::NegativeThreshold, pre, post, align; upsample::Int=1, refractory::Int=post)
    nwf = 0
    wf = zeros(eltype(x), pre+post+1, nwf)
    times = zeros(Int, nwf)
    wf, times
end

end # module
