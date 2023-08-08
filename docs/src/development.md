# Development notes
This module was created using PkgTemplates.jl and the following code:

```
using PkgTemplates


# DIA_CARS

directory = "C:\\Users\\greifenstein\\Documents\\armstrong\\Documents\\Code\\juliacars\\"

t = Template(;
    user = "tuda_rsm/cross-sections/cars/jcars",
    host = "https://git.rwth-aachen.de",
    julia=v"1.7",
    authors="Max Greifenstein, TU Darmstadt RSM",
    dir = directory,
    plugins=[
        GitLabCI(),
        Documenter{GitLabCI}(),
        License(; name="MIT"),
        !TagBot,
        !CompatHelper,
    ]
)

t("DiaCARS.jl")
```

At the moment, the ```codecov``` plugin is deactivated.

## Workflow for changing the code
If not already done, install ```Revise.jl``` and active (either automatically as suggested in the docs or by ```using Revise``` after opening a new REPL (that is, before loading any code!).
Create a new .jl file in the folder that contains ```DiaCARS``` and start working:

```
.
├── DiaCARS/
│   ├── docs
│   ├── src
│   ├── test
│   └── ...
└── testingscript.jl
```

Example content of ```testingscript.jl```:
```
using DiaCARS
T=2000
Species = "N2"
ROI = [2200 2400]
FrequencyOffset = 0
Model = "I"
P = 1

χR,χI,ω = simulateTheoreticalSusceptibility(T,Species,P,Model,ROI,FrequencyOffset)
```
Now start working on the code.
On every run, Revise will use the updated code even without closing the REPL. Check out the documentation of ```Revise.jl``` for the limitations. When satisfied, don't forget to update the documentation, commit and push.