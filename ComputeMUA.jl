module ComputeMUA
using Extract, Sessions, JLD

function compute_mua()
    sess = sessions()
    map(sess) do s
        println(s.date)
        jf = joinpath(Sessions.DATA_DIR, nevstem(s, :change_detection)*"_mua.jld")
        j = jldopen(jf, "w"; mmaparrays=true)
        try
            for ch in channels(s, :IT_VP)
                print("    $ch")
                @load joinpath(dirname(@__FILE__), "metadata.jld") time0
                data = filtereddata(s, :change_detection, ch)
                @time begin
                    threshold = 0.0
                    (_, times) = Extract.extract_wf(data, Extract.NegativeThreshold(threshold), 2, 2, false; refractory=15)
                end
                times = (times-1)/30000 + time0
                write(j, "ch$(ch.ch)", times)
            end
        finally
            close(j)
        end
    end
end
end