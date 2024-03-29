
module WatanabeLuResults

export Msᵂᴸ, M′sᵂᴸ

# ---------------------------------------------------------------------------------------- #
# BOUNDS ON ω=0-CONNECTED (IRREGULAR) TRANSVERSE BANDS: M

# M = 2 cases
M2 = [1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 21, 22, 24, 25, 27, 28, 30, 32, 34, 35,
      37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 65, 66, 68, 70, 74, 75, 77,
      79, 80, 81, 82, 83, 84, 86, 88, 89, 93, 98, 99, 100, 101, 102, 105, 107, 108, 109,
      111, 112, 115, 119, 122, 123, 131, 141, 143, 146, 147, 148, 149, 150, 155, 156, 157,
      160, 162, 164, 166, 168, 174, 175, 177, 183, 187, 189, 191, 195, 196, 199, 200, 203,
      207, 210, 214, 215, 216, 221, 227]

# M ≥ 3 cases
M3 = unique!(
    [23, 71, 72, 97, 121, 139, 140, 197, 204, 209, 211, 217, 225, 226, 229, 69, 139,
     140, 202, 225, 226, 229, 87,  139, 140, 225, 226, 229, 144, 151, 152, 169, 172,
     178, 181, 145, 153, 154, 170, 171, 179, 180]
    )

# M ≥ 4 cases (Note: we manually corrected what I believe were an issue - in the PRL's SI,
#                    key group 4 contains 135:170, but this is not correct; checking against
#                    Bilbao, I found that it should be 135:138, 169, & 170)
M4 = unique!(
    #key group, subgroups
    [4,         11,14,17:20...,26,29,31,33,36,51:64...,76,78,90:92...,94:96...,113,114,
                127:130...,135:138...,169,170,173,176,178,179,182,185,186,193,194,198,
                205,212,213,
     67,        125,129,134,138,224,
     73,        142,206,228,230,
     85,        125,126,129,130,222,
     103,       124,130,
     104,       126,128,222,
     106,       133,135,
     110,       142,228,230,
     116,       124,130,132,138,
     117,       125,127,133,135,
     118,       126,128,134,136,222,224,
     120,       140,142,219,226,228,230,
     158,       165,184,185,188,192,193,
     159,       163,184,186,190,192,194,
     161,       167,218:220...,222,223,226,228,230,
     201,       222,224,
     208,       223,224]
    )
# There can be overlapping groups between the M≥3 and M≥4 cases;
# logically, M≥4 takes precedence over M≥3, so we remove overlaps from M3:
filter!(x->x∉(M4), M3)

# Check that we covered all the space groups:
sort(vcat(M2, M3, M4)) == 1:230 || throw("Unexpected error in loading data")

"""
   Msᵂᴸ

Vector (indexed by sgnum, 1:230) of Watanabe & Lu's bounds on M (irregular photonic band
connectivities).
"""
const Msᵂᴸ = getindex.(sort(vcat(M2 .=> 2, M3 .=> 3, M4 .=> 4), by=first), 2)

# ---------------------------------------------------------------------------------------- #
# BOUNDS ON BANDS NOT CONNECTED TO ω=0 (REGULAR) BANDS: M′

# From Table S1 of Watanabe and Lu; minimum band fillings for ω≠0 bands

"""
   M′sᵂᴸ

Vector (indexed by sgnum, 1:230) of Watanabe & Lu's bounds on M′ (regular photonic band
connectivities).
"""
const M′sᵂᴸ = getindex.(sort(vcat(
    # Symmorphic cases (always M′ = 1)
    [1,   2,   3,   5,   6,   8,   10,  12,  16,  21,  22,  23,  25,  35,  38,  42,  44, 
     47,  65,  69,  71,  75,  79,  81,  82,  83,  87,  89,  97,  99,  107, 111, 115, 119,
     121, 123, 139, 143, 146, 147, 148, 149, 150, 155, 156, 157, 160, 162, 164, 166, 168,
     174, 175, 177, 183, 187, 189, 191, 195, 196, 197, 200, 202, 204, 207, 209, 211, 215, 
     216, 217, 221, 225, 229] .=> 1,
    # Nonsymmorphic cases
    [4  => 2, 39 => 2, 66 => 2, 100 => 2, 129 => 2, 165 => 2, 201 => 2,
     7  => 2, 40 => 2, 67 => 2, 101 => 2, 130 => 4, 167 => 2, 203 => 2,
     9  => 2, 41 => 2, 68 => 2, 102 => 2, 131 => 2, 169 => 6, 205 => 4,
     11 => 2, 43 => 2, 70 => 2, 103 => 2, 132 => 2, 170 => 6, 206 => 4,
     13 => 2, 45 => 2, 72 => 2, 104 => 2, 133 => 4, 171 => 3, 208 => 2,
     14 => 2, 46 => 2, 73 => 4, 105 => 2, 134 => 2, 172 => 3, 210 => 2,
     15 => 2, 48 => 2, 74 => 2, 106 => 4, 135 => 4, 173 => 2, 212 => 4,
     17 => 2, 49 => 2, 76 => 4, 108 => 2, 136 => 2, 176 => 2, 213 => 4,
     18 => 2, 50 => 2, 77 => 2, 109 => 2, 137 => 2, 178 => 6, 214 => 4,
     19 => 4, 51 => 2, 78 => 4, 110 => 4, 138 => 4, 179 => 6, 218 => 2,
     20 => 2, 52 => 4, 80 => 2, 112 => 2, 140 => 2, 180 => 3, 219 => 2,
     24 => 2, 53 => 2, 84 => 2, 113 => 2, 141 => 2, 181 => 3, 220 => 6,
     26 => 2, 54 => 4, 85 => 2, 114 => 2, 142 => 4, 182 => 2, 222 => 2,
     27 => 2, 55 => 2, 86 => 2, 116 => 2, 144 => 3, 184 => 2, 223 => 2,
     28 => 2, 56 => 4, 88 => 2, 117 => 2, 145 => 3, 185 => 2, 224 => 2,
     29 => 4, 57 => 4, 90 => 2, 118 => 2, 151 => 3, 186 => 2, 226 => 2,
     30 => 2, 58 => 2, 91 => 4, 120 => 2, 152 => 3, 188 => 2, 227 => 2,
     31 => 2, 59 => 2, 92 => 4, 122 => 2, 153 => 3, 190 => 2, 228 => 4,
     32 => 2, 60 => 4, 93 => 2, 124 => 2, 154 => 3, 192 => 2, 230 => 8,
     33 => 4, 61 => 4, 94 => 2, 125 => 2, 158 => 2, 193 => 2,
     34 => 2, 62 => 4, 95 => 4, 126 => 2, 159 => 2, 194 => 2,
     36 => 2, 63 => 2, 96 => 4, 127 => 2, 161 => 2, 198 => 4,
     37 => 2, 64 => 2, 98 => 2, 128 => 2, 163 => 2, 199 => 4]
    ), by=first), 2)

end
