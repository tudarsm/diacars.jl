# globally used constants

struct cons
    h::Float32
    c::Float32
    k::Float32
    NA::Float32
    R::Float32
    molvol::Float32
    hc_k::Float32
end

c = cons(6.62607004e-34,29979245800,1.38064852e-23,6.02214086e23,8.314,2.2413e+04,1.4388)

J = 0:100
V = 0:6