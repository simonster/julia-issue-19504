unshift!(LOAD_PATH, dirname(@__FILE__))
reload("ComputeMUA")
ComputeMUA.compute_mua()