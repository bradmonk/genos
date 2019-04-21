function g = HNsigrad(z)

g = HNsig(z) .* (1-HNsig(z));