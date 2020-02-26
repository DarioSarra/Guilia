module Guilia

using Reexport
using TableWidgets
@reexport using Interact
@reexport using Blink
@reexport using ProcessPhotometry
@reexport using JuliaDBMeta
@reexport using DataFrames
@reexport using OnlineStats
@reexport using OffsetArrays
@reexport using ShiftedArrays
@reexport using Recombinase
@reexport using StatsPlots
@reexport using StructArrays
@reexport using IndexedTables
@reexport using WeakRefStrings
@reexport using FillArrays
using Distributions
using Distances
using CategoricalArrays
using KernelDensity: pdf, kde, InterpKDE
using Images
using Tables
using BSON
using FileIO
using OrderedCollections


using Observables
import Observables: AbstractObservable

using Recombinase: offsetrange



include("utilities.jl")
include("body.jl")
include("loading.jl")
include("custom_widgets.jl")
include("plot_attributes.jl")
include("constants.jl")
include(joinpath("Editors","categorize.jl"))
include(joinpath("Editors","separate.jl"))
include(joinpath("New_Analysis","Counting.jl"))
include(joinpath("New_Analysis","Gamma.jl"))
include(joinpath("New_Analysis","Scatter_MeanSD.jl"))
include(joinpath("New_Analysis","Summing.jl"))
include(joinpath("New_Analysis","Delta_Mean.jl"))
include(joinpath("New_Analysis","Normalised_Dist.jl"))
include(joinpath("New_Analysis","Survival.jl"))
include(joinpath("Photometry","sliding_mean.jl"))
include(joinpath("Photometry","regression.jl"))
include(joinpath("Photometry","traces.jl"))
include(joinpath("Photometry","offsets.jl"))
include(joinpath("Photometry","signal_widget.jl"))
include(joinpath("Photometry","construct.jl"))
include("recombinase_gui.jl")

 export launch, launch_signal
 export carica
 export categorify_w, categorize
 export separate_w, separate
 export plot_attributes_w
 export time_offsets, events_offsets
 export collect_traces
 export sliding_f0
 export regress_trace
 export offset_window
 export categorize
 export gui4

end
