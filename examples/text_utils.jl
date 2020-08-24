module TextUtils
using Crystalline: issubdigit, issupdigit, normalizesubsup

export add_content_to_symvec_str_at_kidx, ordinal_indicator, convert_irreplabel2latex

# ----------------------------------------------------------------------------------------

function add_content_to_symvec_str_at_kidx(str::String, kidx::Integer, insert::String)

    parts = replace.(split.(strip.(str, Ref(('[', ']'))), ", "), Ref(" "=>""))
    Γslot = parts[kidx]
    parts[kidx] = insert*(length(Γslot) > 0 ? '+'*Γslot : "")

    return '['*join(parts, ", ")*']'
end

# ----------------------------------------------------------------------------------------

# print the English ordinal indicator of a decimal number (following the conventions in
# https://en.wikipedia.org/wiki/Ordinal_indicator#English)
function ordinal_indicator(n::Integer)
    n < 0 && throw("ordinal numbers must be ≥0")

    lastdigit = n % 10
    if 11 ≤ n ≤ 19          # "teens"     ⇒ 11th, 12th, 13th, ..., 19th
        return "th"
    elseif lastdigit == 1   # ending in 1 ⇒ (..)1st
        return "st"
    elseif lastdigit == 2   # ending in 2 ⇒ (..)2nd
        return "nd"
    elseif lastdigit == 3   # ending in 3 ⇒ (..)3rd
        return "rd"
    else                    # ending in 0, 4:9 ⇒ (..)0th, (..)4th, ..., (..)9th
        return "th"
    end
    
end

# ----------------------------------------------------------------------------------------

function convert_irreplabel2latex(str::AbstractString)
    buf = IOBuffer()
    previous_was_digit = false
    for c in str
        if issubdigit(c)
            previous_was_digit || write(buf, "_{")
            write(buf, normalizesubsup(c)) 
            previous_was_digit = true
        else
            previous_was_digit && (write(buf, '}'); previous_was_digit = false)
            if c ∈ ('⁺', '⁻')
                write(buf, "^{", c=='⁺' ? '+' : '-', "}")
            else
                write(buf, latexifygreek(c))
            end
        end
    end
    previous_was_digit && write(buf, '}')
    return String(take!(buf))
end

const GREEK_SYMBOL_TO_LATEX_MAP = Dict{Char, String}(
    'Γ'=>"\\Gamma", 'Ω' => "\\Omega", 'Σ' => "\\Sigma", 'Λ' => "\\Lambda",
    'Δ'=>"\\Delta"
)

function latexifygreek(c::Char)
    if c ∈ keys(GREEK_SYMBOL_TO_LATEX_MAP)
        return GREEK_SYMBOL_TO_LATEX_MAP[c]
    else
        return string(c)
    end
end

end # module