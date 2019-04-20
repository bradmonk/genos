function g = SIGRAD(z)

g = SIG(z) .* (1-SIG(z));
