# summary of data on minimal phononic band connectivities below fundamental gap
# original data is hosted at: 
#   with time-reversal:    https://gist.github.com/thchr/4c9ac711052738cf07f3857bb02c1947
#   without time-reversal: https://gist.github.com/thchr/76f6ed9861615b726fae2ce303334fe2
# and was calculated with /examples/phononic_fillings.jl

# --- minimal connectivities below fundamental gap ---
# with time-reversal
ms   = [3, 3, 3, 4, 3, 3, 4, 3, 4, 3, 4, 3, 4, 4, 4, 3, 4, 4, 4, 4, 3, 3, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 3, 4, 4, 4, 3, 4, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 3, 4, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 3, 3, 3, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 6, 6, 4, 8, 3, 4, 4, 8, 3, 4, 4, 6, 3, 4, 4, 4, 3, 4, 3, 4, 3, 6, 4, 6, 4, 6, 4, 8, 4, 4, 8, 4, 8, 4, 6, 4, 3, 4, 4, 8, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 6, 6, 3, 6, 3, 6, 3, 6, 3, 6, 3, 6, 6, 6, 6, 6, 3, 3, 6, 3, 6, 6, 6, 6, 6, 3, 6, 6, 6, 3, 6, 3, 6, 3, 6, 6, 6, 3, 3, 3, 4, 4, 3, 4, 3, 4, 3, 4, 4, 3, 4, 3, 4, 3, 4, 4, 4, 3, 3, 3, 6, 6, 6, 3, 6, 6, 4, 3, 6, 4, 8, 3, 8]
# without time-reversal
mstb = [3, 3, 3, 4, 3, 3, 4, 3, 4, 3, 4, 3, 4, 4, 4, 3, 4, 4, 4, 4, 3, 3, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 3, 4, 4, 4, 3, 4, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 3, 4, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 3, 3, 3, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 6, 3, 4, 4, 6, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 8, 4, 6, 4, 3, 4, 4, 6, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 6, 6, 6, 6, 6, 3, 3, 6, 3, 6, 6, 6, 6, 6, 3, 4, 6, 6, 3, 4, 3, 4, 3, 4, 6, 6, 3, 3, 3, 4, 4, 3, 4, 3, 4, 3, 4, 4, 3, 4, 3, 4, 3, 4, 4, 4, 3, 3, 3, 4, 4, 4, 3, 4, 6, 4, 3, 4, 4, 6, 3, 6]


# --- filling-enforced topology below minimal fundamental gap ---
# with time-reversal
fetopo_sgnums = Int[]
# without time-reversal 
fetopo_sgnums_tb = [106, 110, 126, 130, 133, 142, 220, 222, 228, 230]

# In the above cases, we also checked for cases of fragile or mixed fragile+stable filling-
# enforced topology; but didn't find any. All examples are time-broken and are of the stable
# kind (i.e. non-fragile).