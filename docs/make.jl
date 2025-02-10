using DynamicQuantities
using DynamicQuantities.Units
using DynamicQuantities: constructorof, with_type_parameters, dimension_names
using Documenter

DocMeta.setdocmeta!(DynamicQuantities, :DocTestSetup, :(using DynamicQuantities); recursive=true)
doctest(DynamicQuantities)

readme = open(dirname(@__FILE__) * "/../README.md") do io
    read(io, String)
end

# We replace every instance of <img src="IMAGE" ...> with ![](IMAGE).
readme = replace(readme, r"<img src=\"([^\"]+)\"[^>]+>.*" => s"![](\1)")

# Then, we remove any line with "<div" on it:
readme = replace(readme, r"<[/]?div.*" => s"")

# Finally, we read in file docs/src/index_base.md:
index_base = open(dirname(@__FILE__) * "/src/index_base.md") do io
    read(io, String)
end

# And then we create "/src/index.md":
open(dirname(@__FILE__) * "/src/index.md", "w") do io
    write(io, readme)
    write(io, "\n")
    write(io, index_base)
end

makedocs(;
    modules=[DynamicQuantities, DynamicQuantities.Units],
    authors="MilesCranmer <miles.cranmer@gmail.com> and contributors",
    repo="https://github.com/SymbolicML/DynamicQuantities.jl/blob/{commit}{path}#{line}",
    sitename="DynamicQuantities.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ai.damtp.cam.ac.uk/dynamicquantities/stable",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Utilities" => "api.md",
        "Units" => "units.md",
        "Constants" => "constants.md",
        "Symbolic Units" => "symbolic_units.md",
        "Affine Units" => "affine_units.md",
        "Types" => "types.md",
    ],
    warnonly = [:missing_docs]
)

deploydocs(;
    repo="github.com/SymbolicML/DynamicQuantities.jl",
    devbranch="main"
)

# Mirror to DAMTP:
ENV["DOCUMENTER_KEY"] = ENV["DOCUMENTER_KEY_CAM"]
ENV["GITHUB_REPOSITORY"] = "ai-damtp-cam-ac-uk/dynamicquantities.git"
deploydocs(;
    repo="github.com/ai-damtp-cam-ac-uk/dynamicquantities.git",
    devbranch="main"
)
