var documenterSearchIndex = {"docs":
[{"location":"visual/#Visualization","page":"Visualize","title":"Visualization","text":"","category":"section"},{"location":"visual/#Coordinate-representation","page":"Visualize","title":"Coordinate representation","text":"","category":"section"},{"location":"visual/","page":"Visualize","title":"Visualize","text":"TODO","category":"page"},{"location":"visual/#Data-processing","page":"Visualize","title":"Data processing","text":"","category":"section"},{"location":"visual/","page":"Visualize","title":"Visualize","text":"TODO","category":"page"},{"location":"visual/#Automatic-plotting","page":"Visualize","title":"Automatic plotting","text":"","category":"section"},{"location":"visual/","page":"Visualize","title":"Visualize","text":"TODO","category":"page"},{"location":"operators/#Linear-operators","page":"Linear operators","title":"Linear operators","text":"","category":"section"},{"location":"operators/#Hamiltonian-generation","page":"Linear operators","title":"Hamiltonian generation","text":"","category":"section"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"The hamiltonian function can be used to generate a Chern insulator hamiltonian.  The only essential parameter is m_lattice - it is a matrix with the size of the Chern insulator lattice which contains the value of the m parameter for each site. Its type must be CoordinateRepr. There is an alternative way to set the m_lattice parameter - if its type is a Matrix, an additional argument of type Symbol is required - a representation specifier (see Coordinate Representation for more detail).","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"Here is an example of both usages:","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"using TopologicalMarkers","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"m_lattice = ones(15, 15)\nm_lattice[6:10, 7:11] .= -1\nm_repr = CoordinateRepr(m_lattice, :n)\nH1 = hamiltonian(m_repr)\nH2 = hamiltonian(m_lattice, :n)\n\nH1 == H2 # These are equal ways to specify the m_lattice parameter","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"There is an ability to isolate some sites from other ones. To do that, a zone_mapping matrix is required.  It is a CoordinateRepr{Symbol} object, which maps lattice sites to symbols, where different symbols mean different zones,  and hoppings between sites in different zones are set to zero. You can define the mapping in the hamiltonian function call, or apply it to an existing hamiltonian matrix using the zones! function:","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"m_lattice = ones(15, 15)\n\nmapping = CoordinateRepr(fill(:zone1, 15, 15))\nmapping[6:10, 6:10] .= :zone2\n\nH1 = hamiltonian(m_lattice, :n, zone_mapping = mapping)\nH2 = hamiltonian(m_lattice, :n)\nzones!(H2, mapping)\n\nH1 == H2 # These are equal ways to specify the zone mapping","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"You also can apply some magnetic field to the Chern insulator.  This can be done by passing a function that takes the position vector and returns the gauge vector potential. You can write your own function  use one of the magnetic field macros defined: ","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"@landau for Landau gauge field\n@symm for symmetric gauge field\n@flux for flux quantum","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"You can define the field in the hamiltonian function or apply it to an existing hamiltonian matrix using the field! function:","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"m_lattice = ones(15, 15)\nB = 0.1\n\nH0 = hamiltonian(m_lattice, :n, field = (x -> [0, x[1] * B, 0]))\nH1 = hamiltonian(m_lattice, :n, field = @landau(B))\nH2 = hamiltonian(m_lattice, :n)\nfield!(H2, @landau(B))\n\nH0 == H1 == H2 # These are equal ways to specify the field","category":"page"},{"location":"operators/#Other-linear-operators","page":"Linear operators","title":"Other linear operators","text":"","category":"section"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"To calculate a density matrix for zero temperature, you can use the filled_projector function.  It accepts the hamiltonian matrix and returns a density matrix. You can also specify the Fermi level (which is zero by default).","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"You can also generate the coordinate operators using the coord_operators function. It takes an optional parameter symm which defines if the operators are defined symmetrically (i. e. the point with (0, 0) coordinates is in the center of the lattice) or not (i. e. the point with (0, 0) coordinates is in the bottom-left angle of the lattice).","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"Here is an example of usage of both functions to calculate the local Chern marker:","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"using Plots, LinearAlgebra\n\nm_lattice = ones(15, 15)\nm_lattice[6:10, 7:11] .= -1\nH = hamiltonian(m_lattice, :n)\n\nP = filled_projector(H)\nX, Y = coord_operators()\nc = 4pi * im * P * X * (I - P) * Y * P\n\nheatmap(heatmap_data(c), color=:viridis)","category":"page"},{"location":"operators/#Electric-current","page":"Linear operators","title":"Electric current","text":"","category":"section"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"You can calculate electric currents using the currents function.  It accepts the hamiltonian and the density matrix to produce a matrix with currents. Here is an example of electric currents in a magnetic field:","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"using Plots\n\nB = 0.1\nH = hamiltonian(ones(15, 15), :n, field=@symm(B))\nP = filled_projector(H)\nps, qs = quiver_data(currents(H, P) * 10 / B)\nquiver(ps, quiver = qs)","category":"page"},{"location":"operators/#LCM-current","page":"Linear operators","title":"LCM current","text":"","category":"section"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"There are several formulas which allow us to calculate currents of the local Chern marker. These formulas comply to the continuity equation, but, unlike electric currents, they are not localized. In this program, there are several macros that can interest you:","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"Modules = [TopologicalMarkers]\nPages = [joinpath(\"formulas\", \"lcm_currents.jl\")]\nFilter = f -> nameof(f) != Symbol(\"@currents\")","category":"page"},{"location":"operators/#TopologicalMarkers.@J-NTuple{4, Any}","page":"Linear operators","title":"TopologicalMarkers.@J","text":"@J(H, P, X, Y)\n\n\n\n\n\n","category":"macro"},{"location":"operators/#TopologicalMarkers.@J_best-NTuple{4, Any}","page":"Linear operators","title":"TopologicalMarkers.@J_best","text":"@J_best(H, P, X, Y)\n\n\n\n\n\n","category":"macro"},{"location":"operators/#TopologicalMarkers.@J_c-NTuple{4, Any}","page":"Linear operators","title":"TopologicalMarkers.@J_c","text":"@J_c(H, P, X, Y)\n\n\n\n\n\n","category":"macro"},{"location":"operators/#TopologicalMarkers.@J_eq-NTuple{4, Any}","page":"Linear operators","title":"TopologicalMarkers.@J_eq","text":"@J_eq(H, P, X, Y)\n\n\n\n\n\n","category":"macro"},{"location":"operators/#TopologicalMarkers.@J_inv-NTuple{4, Any}","page":"Linear operators","title":"TopologicalMarkers.@J_inv","text":"@J_inv(H, P, X, Y)\n\n\n\n\n\n","category":"macro"},{"location":"operators/#TopologicalMarkers.@J_m-NTuple{4, Any}","page":"Linear operators","title":"TopologicalMarkers.@J_m","text":"@J_m(H, P, X, Y)\n\n\n\n\n\n","category":"macro"},{"location":"operators/#TopologicalMarkers.@J_m_inv-NTuple{4, Any}","page":"Linear operators","title":"TopologicalMarkers.@J_m_inv","text":"@J_m_inv(H, P, X, Y)\n\n\n\n\n\n","category":"macro"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"Each one takes matrices of the hamiltonian, the density and coordinate operators and returns a function that calculates a current given indices of 2 sites.","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"To generate a matrix with currents, use the @currents macro:","category":"page"},{"location":"operators/","page":"Linear operators","title":"Linear operators","text":"@currents","category":"page"},{"location":"operators/#TopologicalMarkers.@currents","page":"Linear operators","title":"TopologicalMarkers.@currents","text":"@currents(currents_lambda)\n\nGenerates a matrix with currents, given a lambda/macrocall that takes lattice site indices and returns the current value\n\nExample usage:\n\ncurrents_mat = @currents @J H P X Y\n\n\n\n\n\n","category":"macro"},{"location":"scope/#TopologicalMarkers.jl","page":"Library","title":"TopologicalMarkers.jl","text":"","category":"section"},{"location":"scope/","page":"Library","title":"Library","text":"Modules = [TopologicalMarkers]","category":"page"},{"location":"scope/#Hamiltonian-matrix","page":"Library","title":"Hamiltonian matrix","text":"","category":"section"},{"location":"scope/","page":"Library","title":"Library","text":"hamiltonian\nfield!\nzones!\n@landau\n@symm\n@flux","category":"page"},{"location":"scope/#TopologicalMarkers.hamiltonian","page":"Library","title":"TopologicalMarkers.hamiltonian","text":"hamiltonian{T}((m_repr | m_lattice, repr_spec); <keyword arguments>)\n\nGenerates a Hamiltonian operator for a Chern insulator using the following formula\n\nhatH = \nsum_i m_i c^dagger_i sigma_z c_i + \nsum_x-links c^dagger_i fracsigma_z - i sigma_x2 c_j + \nsum_y-links c^dagger_i fracsigma_z - i sigma_y2 c_j + \nh c\n\nArguments\n\nm_repr: The value of m on different sites, in CoordinateRepr format\n\nAlternatively, pass a matrix and a representation specifier (see CoordinateRepr for more information)\n\npbc: Periodic boundary conditions. A Tuple{Bool, Bool}, each element sets boundary conditions for the horizontal and vertical edge respectively. Default is (false, false)\nzones: A matrix with elements of arbitrary type, which maps sites to isolated zones. The hopping members between different zones are erased. There are no isolated zones by default\nfield: A function/lambda that takes two coordinates and returns the vector potential of the magnetic field. Used to calculate phase factors on hoppings. There is no magnetic field by default\n\n\n\n\n\n","category":"function"},{"location":"scope/#TopologicalMarkers.field!","page":"Library","title":"TopologicalMarkers.field!","text":"field!(H, A[, lattice_size]; intervals = 10)\n\nApplies magnetic field to specified hamiltonian. In other words, all hoppings get multiplied on a specific phase factor that can be calculated using Peierls substitution:\n\nvarphi_ij = 2pi int_i^j A(r) cdot dr\n\nThis integral is calculated explicitly for every hopping, using the A function.\n\nArguments\n\nH: the hamiltonian matrix\nA: a function that takes a Vector representing a point and returns a Vector representing the vector potential in that point\nlattice_size: the size of the lattice the hamiltonian is defined for. If not provided, this function will use the value for the hamiltonian matrix that was created last\nintervals: the number of intervals to use when calculating the Peierls substitution phase factor\n\n\n\n\n\n","category":"function"},{"location":"scope/#TopologicalMarkers.zones!","page":"Library","title":"TopologicalMarkers.zones!","text":"zones!(H, zone_mapping[, repr_spec])\n\nDivides the Chern insulator hamiltonian into several unconnected zones. The hoppings between these zones are erased.\n\nArguments\n\nH: the hamiltonian matrix\nzone_mapping: an AbstractMatrix{Symbol} or CoordinateRepr{Symbol}. Each site is mapped to a symbol, different symbols mean different zones\nrepr_spec: if zone_mapping is an AbstractMatrix{Symbol}, this argument is a representation specifier (see CoordinateRepr docs for more information)\n\n\n\n\n\n","category":"function"},{"location":"scope/#TopologicalMarkers.@landau","page":"Library","title":"TopologicalMarkers.@landau","text":"@landau(B)\n\nGenerates a function that returns the Landau gauge vector potential. This can be used as an argument for the field! function.\n\nArguments\n\nB: the magnetic field value\n\n\n\n\n\n","category":"macro"},{"location":"scope/#TopologicalMarkers.@symm","page":"Library","title":"TopologicalMarkers.@symm","text":"@symm(B[, center])\n\nGenerates a function that returns the symmetric gauge vector potential. This can be used as an argument for the field! function.\n\nArguments\n\nB: the value of magnetic field\ncenter: the center of symmetry of the vector potential\n\n\n\n\n\n","category":"macro"},{"location":"scope/#TopologicalMarkers.@flux","page":"Library","title":"TopologicalMarkers.@flux","text":"flux(Φ[, point])\n\nGenerates a function that returns the vector potential for a flux quantum. This can be used as an argument for the field! function.\n\nArguments\n\nΦ: the value of magnetic field\npoint: the point where the flux is located\n\n\n\n\n\n","category":"macro"},{"location":"scope/#Other-linear-operators","page":"Library","title":"Other linear operators","text":"","category":"section"},{"location":"scope/","page":"Library","title":"Library","text":"coord_operators\nfilled_projector\ncurrents","category":"page"},{"location":"scope/#TopologicalMarkers.coord_operators","page":"Library","title":"TopologicalMarkers.coord_operators","text":"coord_operators([lattice_size]; <keyword arguments>)\n\nReturns a tuple of coordinate operators (i. e. hatX and hatY).\n\nArguments\n\nlattice_size: the size of the lattice\nsymmetric: defines if the operators are symmetricall (in other words, if the center of the lattice corresponds to (0, 0))\n\n\n\n\n\n","category":"function"},{"location":"scope/#TopologicalMarkers.filled_projector","page":"Library","title":"TopologicalMarkers.filled_projector","text":"filled_projector(H[, energy_thr=0])\n\nReturns a projector onto the filled states (in other words, a ground state density matrix).\n\nArguments\n\nH: the hamiltonian matrix\nenergy_thr: the Fermi energy level\n\n\n\n\n\n","category":"function"},{"location":"scope/#TopologicalMarkers.currents","page":"Library","title":"TopologicalMarkers.currents","text":"currents(H, P[, lattice_size])\n\nReturns a skew-symmetric matrix of electric currents between sites.\n\nArguments\n\nH: the hamiltonian matrix.\nP: the density matrix.\nlattice_size: the size of the lattice the hamiltonian is defined for.\n\n\n\n\n\n","category":"function"},{"location":"scope/#Unitary-evolution","page":"Library","title":"Unitary evolution","text":"","category":"section"},{"location":"scope/","page":"Library","title":"Library","text":"evolution_operator\n@evolution","category":"page"},{"location":"scope/#TopologicalMarkers.evolution_operator","page":"Library","title":"TopologicalMarkers.evolution_operator","text":"evolution_operator(H, t)\n\nCalculates the unitary evolution operator using the formula\n\n$ \\mathcal{U}(t) = e^{-\\frac{1}{i\\hbar} \\hat{H} t} $\n\nArguments\n\nH: the hamiltonian matrix\nt: the evolution time\n\n\n\n\n\n","category":"function"},{"location":"scope/#TopologicalMarkers.@evolution","page":"Library","title":"TopologicalMarkers.@evolution","text":"@evolution [rules...] for_loop\n\nGenerates an environment with defined hamiltonian and density matrices that evolve by certain laws. See Unitary evolution for more details.\n\n\n\n\n\n","category":"macro"},{"location":"scope/#Data-visualization","page":"Library","title":"Data visualization","text":"","category":"section"},{"location":"scope/","page":"Library","title":"Library","text":"CoordinateRepr\nheatmap_data\nquiver_data","category":"page"},{"location":"scope/#TopologicalMarkers.CoordinateRepr","page":"Library","title":"TopologicalMarkers.CoordinateRepr","text":"CoordinateRepr{T}\n\nA wrapper class for matrices to be conveniently plotted. \n\nArguments\n\nlattice: A matrix representing some quantity defined on the lattice (e. g. the LCM)\nrepr_spec: A symbol defining the way how the lattice sites match the matrix values:\n:c or :coord: The value for (x, y) site is A[x, y]. This is the default value\n:n of :natural: If you print out the matrix as you usually do, and then imagine \na coordinate system with its center in the bottom-left corner, this will be the mapping between   sites and matrix values\n\n\n\n\n\n","category":"type"},{"location":"scope/#TopologicalMarkers.heatmap_data","page":"Library","title":"TopologicalMarkers.heatmap_data","text":"heatmap_data(op[, lattice_size])\n\nGenerates a CoordinateRepr for values like langle r  hatmathcalO  r rangle.\n\nArguments\n\nop: the operator to find values for\nlattice_size: the size of the lattice\n\n\n\n\n\n","category":"function"},{"location":"scope/#TopologicalMarkers.quiver_data","page":"Library","title":"TopologicalMarkers.quiver_data","text":"quiver_data(currents_mat[, lattice_size]; <keyword arguments>)\n\nGenerates data for a quiver plot using a matrix with currents. The output is a tuple of two vectors with equal length: one contains arrow origins, the other one - arrow vectors.\n\nArguments\n\ncurrents_mat: a matrix with currents\nlattice_size: the size of the lattice\nthreshold: minimum value of the current to be put to output. Default is 0.1\ndist_threshold: maximum distance between sites for which the current will be evaluated. Infinite by default\nxlims and ylims: limit the area in which the currents will be evaluated. Infinite by default\n\n\n\n\n\n","category":"function"},{"location":"scope/#Plotting-shorthands","page":"Library","title":"Plotting shorthands","text":"","category":"section"},{"location":"scope/","page":"Library","title":"Library","text":"These functions may come in handy when you want to plot many marker heatmaps and electric current diagrams in one line.","category":"page"},{"location":"scope/","page":"Library","title":"Library","text":"plot_boundaries!\nplot_marker!\nplot_auto","category":"page"},{"location":"scope/#TopologicalMarkers.plot_boundaries!","page":"Library","title":"TopologicalMarkers.plot_boundaries!","text":"plot_boundaries!([pl, ]zone_mapping; <keyword arguments>)\n\nDraws boundaries between different zones described by the zone_mapping matrix.\n\nArguments\n\npl: a Plots.Plot object to visualize data on\nzone_mapping: a CoordinateRepr object that represents zone mapping. The boundaries between different zones will be drawn\n\nAll keyword arguments will be passed to the plot! function used for drawing - this can be used to change the line thickness or style, for example.\n\n\n\n\n\n","category":"function"},{"location":"scope/#TopologicalMarkers.plot_marker!","page":"Library","title":"TopologicalMarkers.plot_marker!","text":"plot_marker!(pl; <keyword arguments>)\n\nPlots complicated marker data series (heatmap, boundaries, quiver) on a single figure.\n\nArguments\n\npl: a Plots.Plot object to visualize data on\nhmap: data to be visualized on a heatmap. It can be a CoordinateRepr object (then it will be plotted directly) \n\nor a linear operator matrix (then the CoordinateRepr will be generated automatically)\n\nzone_mapping: a CoordinateRepr object that represents zone mapping. The boundaries between different zones will be drawn.\ncurrents: a matrix containing currents between sites\nxlims and ylims: objects of type Tuple{Int, Int} that define the limits of the x- and y- axes respectively\n\nAll keyword arguments with different prefixes are passed th the plot! function:\n\nhmap for the heatmap\nbounds for the boundaries\ncurrents for the quiver\n\nThis can be used to style the plot:\n\nplot_marker!(..., hmapclims=(-3, 3), boundsstyle=:dot, :currentscolor=:green)\n\n\n\n\n\n","category":"function"},{"location":"scope/#TopologicalMarkers.plot_auto","page":"Library","title":"TopologicalMarkers.plot_auto","text":"plot_auto(<arguments>; <keyword arguments>)\n\nPlots multiple heatmaps, slices or currents simultaneously.\n\nThe subplots are automatically arranged into an optimal layout.\n\nArguments\n\nEach argument can be either a CoordinateRepr object or a chain of pairs.\n\n\n\n\n\n","category":"function"},{"location":"evolution/#Unitary-evolution","page":"Unitary evolution","title":"Unitary evolution","text":"","category":"section"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"Suppose we want to study the behavior of some quantum system in time-dependent conditions. We can use the unitary evolution operator to describe how the density matrix depends on time:","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"mathcalU(t) = Tleft e^frac1ihbar int_t_0^t hatH(tau) dtau righthspace05cm\nmathcalP(t) = mathcalU(t) mathcalP_0 mathcalU^dagger (t)","category":"page"},{"location":"evolution/#The-static-evolution-function","page":"Unitary evolution","title":"The static evolution function","text":"","category":"section"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"The evolution_operator function calculates the evolution operator for a time-independent hamiltonian.  It also stores the output data, so if you call it multiple times with the same input, the output will be calculated only once.","category":"page"},{"location":"evolution/#The-evolution-macro","page":"Unitary evolution","title":"The evolution macro","text":"","category":"section"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"This macro can be quite useful if your hamiltonian depends on time or if there are multiple hamiltonians in your experiment. What you need is a function that takes the time and returns the matrix of the hamiltonian at that moment of time:","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"h(t) = hamiltonian(ms, field=@symm(Bf * t / τ))","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"Here h(t) describes the magnetic field being adiabatically turned on.","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"Let's take a look at the example 2 from the Home section. The @evolution macro creates a block where the hamiltonian and density matrices are evaluated for the given time interval.  It takes two arguments - a list/vector with evolution specifiers and a for-loop that iterates over the time interval:","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"a = Animation()\n\n@evolution [\n    :ham => h => H,\n    P0 => h => P\n] for time in time_domain\n    cur = currents(H, P)\n    plot_auto(\"Local density\" => P => cur * 40, hmapclims=(0.98, 1.02))\n    frame(a)\nend","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"Let us make it clear what an evolution specifier is: it is a pair chain with three arguments:","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"arg1 => arg2 => arg3","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"There are two options what these arguments can be:","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"The first argument is a density matrix at the moment t = 0. The second is the function describing the hamiltonian time dependence, and the third one is the name for the density matrix at the t = texttime moment - a variable with such a name will be defined inside the loop body.\nThe first argument is the :ham symbol. Then the second one is still the hamiltonian function, and the third one is the name for the hamiltonian matrix at the t = texttime moment.","category":"page"},{"location":"evolution/","page":"Unitary evolution","title":"Unitary evolution","text":"So, in the example before in the for-loop body H stands for h(time), and P is the evolved P0 density matrix.","category":"page"},{"location":"#TopologicalMarkers.jl","page":"Home","title":"TopologicalMarkers.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A package that simplifies calculation of different topological markers.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install it, simply copy this line to julia's REPL and execute it:","category":"page"},{"location":"","page":"Home","title":"Home","text":"]add https://github.com/aryavorskiy/TopologicalMarkers/","category":"page"},{"location":"#Examples","page":"Home","title":"Examples","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Let us take a Chern insulator, set the m parameter to -1 in the middle of the lattice and 1 everywhere else, and then evaluate the local Chern marker using both traditional and Streda formulas (i. e. the linear response of the local density to the magnetic field):","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TopologicalMarkers\nusing Plots\n\nm_lattice = ones(25, 25)\nm_lattice[11:15, 11:15] .= -1\nH = hamiltonian(m_lattice, :c)\nP = filled_projector(H)\nX, Y = coord_operators()\nch = -4π * im * P * X * P * Y * P\n\nB = 0.01\nHb = hamiltonian(m_lattice, :c, field=@landau(B))\nPb = filled_projector(Hb)\nstr = (Pb - P) / B\n\nplot_auto(\"LCM\" => ch, \"Streda\" => str, \n    hmapclims=(-1.5, 1.5), currentscolor=:yellow, control_site=(13, 13), markercolor=:brown)","category":"page"},{"location":"","page":"Home","title":"Home","text":"See (Hamiltonian generation)[@ref] and (Visualization)[@ref] for detailed explanation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Another good example is a problem where unitary evolution is used.  Here we create an animation of the local density changing during adiabatic magnetic field-on:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TopologicalMarkers\nusing Plots\n\nms = CoordinateRepr(ones(15, 15))\nBf = 0.01\nτ = 30\ntime_domain = 0:0.5:τ\n\nH0 = hamiltonian(ms)\nP0 = filled_projector(H0)\nh(t) = hamiltonian(ms, field=@symm(Bf * t / τ))\na = Animation()\n@evolution [\n    :ham => h => H,\n    P0 => h => P\n] for t in time_domain\n    cur = currents(H, P)\n    plot_auto(\"Local density\" => P => cur * 100, hmapclims=(0.9, 1.1))\n    frame(a)\nend\n\ngif(a, \"example_animation.gif\", fps=10)","category":"page"},{"location":"","page":"Home","title":"Home","text":"See (Unitary evolution)[@ref] for detailed explanation.","category":"page"}]
}
