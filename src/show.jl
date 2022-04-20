function Base.show(io::IO, ::MIME"text/plain", v::Polar; usepi=true)
    println(io, typeof(v), ":")
    println(io, " r: ", r(v))
    # TODO: do I want to have the factor of π included?
    if usepi
        print(io, " θ: ", theta(v)/π, "π")
    else
        print(io, " θ: ", theta(v))
    end
end

function Base.show(io::IO, ::MIME"text/plain", v::Spherical; usepi=true)
    println(io, typeof(v), ":")
    println(io, " r: ", r(v))
    # TODO: do I want to have the factor of π included?
    if usepi
        println(io, " θ: ", theta(v)/π, "π")
        print(io, " ϕ: ", phi(v)/π, "π")
    else
        println(io, " θ: ", theta(v))
        print(io, " ϕ: ", phi(v))
    end
end
