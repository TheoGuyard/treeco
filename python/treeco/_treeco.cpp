#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "treeco/Geometry.hpp"
#include "treeco/LDTree.hpp"
#include "treeco/Problem.hpp"
#include "treeco/Problem/Explicit.hpp"
#include "treeco/Problem/Knapsack.hpp"
#include "treeco/Problem/Maxcut.hpp"
#include "treeco/Problem/Tsp.hpp"
#include "treeco/IO.hpp"
#include "treeco/Dynprog.hpp"
#include "treeco/Types.hpp"

namespace py = pybind11;

using namespace treeco;


// Forward declarations
void bind_module_geometry(py::module_ &m);
void bind_module_problem(py::module_ &m);
void bind_module_voronoi(py::module_ &m);
void bind_module_dynprog(py::module_ &m);
void bind_module_tree(py::module_ &m);
void bind_module_ldtree(py::module_ &m);
void bind_module_io(py::module_ &m);


PYBIND11_MODULE(_treeco, treeco) {
    treeco.doc() = "Treeco - Linear Decision Tree Policies for Combinatorial Optimization";

    py::module_ geometry = treeco.def_submodule("geometry");
    py::module_ problem  = treeco.def_submodule("problem");
    py::module_ voronoi  = treeco.def_submodule("voronoi");
    py::module_ dynprog  = treeco.def_submodule("dynprog");
    py::module_ tree     = treeco.def_submodule("tree");
    py::module_ ldtree   = treeco.def_submodule("ldtree");
    py::module_ io       = treeco.def_submodule("io");

    bind_module_geometry(geometry);
    bind_module_problem(problem);
    bind_module_voronoi(voronoi);
    bind_module_dynprog(dynprog);
    bind_module_tree(tree);
    bind_module_ldtree(ldtree);
    bind_module_io(io);
}

void bind_module_geometry(py::module_ &m) {
    
    py::enum_<Relation>(m, "Relation")
        .value("lt", Relation::LT)
        .value("le", Relation::LE)
        .value("eq", Relation::EQ)
        .value("ge", Relation::GE)
        .value("gt", Relation::GT)
        .value("rt", Relation::RT)
        .value("rf", Relation::RF)
        ;

    py::class_<ConstrData>(m, "ConstrData")
        .def_property_readonly("a", [](const ConstrData &self) { return std::get<0>(self); })
        .def_property_readonly("b", [](const ConstrData &self) { return std::get<1>(self); })
        .def_property_readonly("rel", [](const ConstrData &self) { return std::get<2>(self); })
        ;

    py::bind_vector<Domain>(m, "Domain");

    m.def("positive_orthant", &positiveOrthant, py::arg("dimension"));
    m.def("non_negative_orthant", &nonNegativeOrthant, py::arg("dimension"));
    m.def("negative_orthant", &negativeOrthant, py::arg("dimension"));
    m.def("non_positive_orthant", &nonPositiveOrthant, py::arg("dimension"));
}

void bind_module_problem(py::module_ &m) {
    
    py::class_<Problem, std::shared_ptr<Problem>>(m, "Problem")
        .def_property_readonly("dimension", &Problem::dimension)
        .def("get_feasible_set", &Problem::getFeasibleSet)
        .def("get_cost_domain", &Problem::getCostDomain)
        .def("sample_cost", static_cast<RealVector (Problem::*)() const>(&Problem::sampleCost))
        .def("sample_cost", static_cast<RealVector (Problem::*)(Index) const>(&Problem::sampleCost), py::arg("seed"))
        ;

    py::class_<Explicit, Problem, std::shared_ptr<Explicit>>(m, "Explicit")
        .def(py::init<const std::vector<BinaryVector>&>(), py::arg("feasible_set"))
        .def("get_feasible_set", &Explicit::getFeasibleSet, py::call_guard<py::gil_scoped_release>());

    py::class_<Knapsack, Problem, std::shared_ptr<Knapsack>>(m, "Knapsack")
        .def(py::init<Index>(), py::arg("num_items"))
        .def(py::init<const RealVector&, double>(), py::arg("weights"), py::arg("capacity"))
        .def_property_readonly("num_items", &Knapsack::numItems)
        .def_property_readonly("weights", &Knapsack::weights, py::return_value_policy::copy)
        .def_property_readonly("capacity", &Knapsack::capacity, py::return_value_policy::copy)
        .def("get_feasible_set", &Knapsack::getFeasibleSet, py::call_guard<py::gil_scoped_release>())
        .def("get_cost_domain", &Knapsack::getCostDomain)
        .def("sample_cost", static_cast<RealVector (Problem::*)() const>(&Problem::sampleCost))
        .def("sample_cost", static_cast<RealVector (Knapsack::*)(Index) const>(&Knapsack::sampleCost), py::arg("seed"))
        ;

    py::class_<Maxcut, Problem, std::shared_ptr<Maxcut>>(m, "Maxcut")
        .def(py::init<Index>(), py::arg("num_vertices"))
        .def_property_readonly("num_vertices", &Maxcut::numVertices)
        .def("get_feasible_set", &Maxcut::getFeasibleSet, py::call_guard<py::gil_scoped_release>())
        .def("get_cost_domain", &Maxcut::getCostDomain)
        .def("sample_cost", static_cast<RealVector (Problem::*)() const>(&Problem::sampleCost))
        .def("sample_cost", static_cast<RealVector (Maxcut::*)(Index) const>(&Maxcut::sampleCost), py::arg("seed"))
        .def("edge_to_index", &Maxcut::edgeToIndex)
        .def("index_to_edge", &Maxcut::indexToEdge);

    py::class_<Tsp, Problem, std::shared_ptr<Tsp>>(m, "Tsp")
        .def(py::init<Index>(), py::arg("num_cities"))
        .def_property_readonly("num_cities", &Tsp::numCities)
        .def("get_feasible_set", &Tsp::getFeasibleSet, py::call_guard<py::gil_scoped_release>())
        .def("get_cost_domain", &Tsp::getCostDomain)
        .def("sample_cost", static_cast<RealVector (Problem::*)() const>(&Problem::sampleCost))
        .def("sample_cost", static_cast<RealVector (Tsp::*)(Index) const>(&Tsp::sampleCost), py::arg("seed"))
        .def("edge_to_index", &Tsp::edgeToIndex)
        .def("index_to_edge", &Tsp::indexToEdge);
}

void bind_module_voronoi(py::module_ &m) {

    py::class_<VoronoiStats>(m, "VoronoiStats")
        .def_readonly("build_time", &VoronoiStats::buildTime)
        .def_readonly("lp_solved", &VoronoiStats::lpSolved)
        ;

    py::class_<Voronoi>(m, "Voronoi")
        .def_property_readonly("dim_points", &Voronoi::dimPoints)
        .def_property_readonly("num_points", &Voronoi::numPoints)
        .def_property_readonly("num_splits", &Voronoi::numSplits)
        .def_property_readonly("num_faces", &Voronoi::numFaces)
        .def_property_readonly("num_edges", &Voronoi::numEdges)
        .def_property_readonly("stats", &Voronoi::stats, py::return_value_policy::reference_internal)
        ;
}

void bind_module_dynprog(py::module_ &m) {
    py::enum_<Exploration>(m, "Exploration")
        .value("greedy", Exploration::GREEDY)
        .value("iterative", Exploration::ITERATIVE)
        .value("exhaustive", Exploration::EXHAUSTIVE)
        ;

    py::enum_<Branching>(m, "Branching")
        .value("ternary", Branching::TERNARY)
        .value("binary", Branching::BINARY)
        ;

    py::enum_<LowerBounding>(m, "LowerBounding")
        .value("fixed", LowerBounding::FIXED)
        .value("backtrack", LowerBounding::BACKTRACK)
        ;

    py::enum_<Positioning>(m, "Positioning")
        .value("online", Positioning::ONLINE)
        .value("precompute", Positioning::PRECOMPUTE)
        ;

    py::enum_<SplitSelection>(m, "SplitSelection")
        .value("all", SplitSelection::ALL)
        .value("sampling", SplitSelection::SAMPLING)
        ;

    py::enum_<SplitScoring>(m, "SplitScoring")
        .value("variance", SplitScoring::VARIANCE)
        .value("entropy", SplitScoring::ENTROPY)
        .value("minmax", SplitScoring::MINMAX)
        .value("none", SplitScoring::NONE)
        .value("random", SplitScoring::RANDOM)
        ;

    py::class_<DynprogStats>(m, "DynprogStats")
        .def_readonly("run_time", &DynprogStats::runTime)
        .def_readonly("num_iters", &DynprogStats::numIters)
        .def_readonly("num_evals", &DynprogStats::numEvals)
        .def_readonly("num_states", &DynprogStats::numStates)
        .def_readonly("num_states_built", &DynprogStats::numStatesBuilt)
        .def_readonly("num_states_closed", &DynprogStats::numStatesClosed)
        .def_readonly("num_states_leafed", &DynprogStats::numStatesLeafed)
        .def_readonly("num_states_pruned", &DynprogStats::numStatesPruned)
        .def_readonly("num_feasibility_checks", &DynprogStats::numFeasibilityChecks)
        .def_readonly("optimal_depth", &DynprogStats::optimalDepth)
        ;

    py::bind_vector<DynprogLogs>(m, "DynprogLogs");
}

void bind_module_tree(py::module_ &m) {

    py::class_<TreeStats>(m, "TreeStats")
        .def_readonly("build_time", &TreeStats::buildTime)
        .def_readonly("dynprog_stats", &TreeStats::dynprogStats, py::return_value_policy::reference_internal)
        .def_readonly("dynprog_logs", &TreeStats::dynprogLogs, py::return_value_policy::reference_internal)
        ;

    py::class_<Tree>(m, "Tree")
        .def_property_readonly("size", &Tree::size)
        .def_property_readonly("width", &Tree::width)
        .def_property_readonly("depth", &Tree::depth)
        .def_property_readonly("stats", &Tree::stats, py::return_value_policy::reference_internal)
        ;
}

void bind_module_ldtree(py::module_ &m) {
    
    py::class_<LDTreeStats>(m, "LDTreeStats")
        .def_readonly("build_time", &LDTreeStats::buildTime)
        ;

    py::class_<LDTree>(m, "LDTree")
        .def(py::init<const std::string&, const std::string&>(),
            py::arg("file_points"),
            py::arg("file_domain") = ""
        )
        .def(py::init<const std::vector<BinaryVector>&, const Domain&>(),
            py::arg("points"),
            py::arg("domain") = Domain()
        )
        .def("build", [](
                LDTree&         self,
                bool            verbose,
                double          logInterval,
                bool            logSave,
                double          timeLimit,
                double          tolerance,
                bool            deduplicate,
                bool            filterChecks,
                Exploration     exploration,
                Branching       branching,
                LowerBounding   lowerBounding,
                Positioning     positioning,
                SplitSelection  splitSelection,
                SplitScoring    splitScoring,
                Index           randomSeed
            ) {
                self.build(
                    verbose,
                    verbose ? &std::cout : nullptr,
                    logInterval,
                    logSave,
                    timeLimit,
                    tolerance,
                    deduplicate,
                    filterChecks,
                    exploration,
                    branching,
                    lowerBounding,
                    positioning,
                    splitSelection,
                    splitScoring,
                    randomSeed
                );
            },
            py::arg("verbose")          = false,
            py::arg("log_interval")     = 5.0,
            py::arg("log_save")         = true,
            py::arg("time_limit")       = std::numeric_limits<double>::infinity(),
            py::arg("tolerance")        = 1e-8,
            py::arg("deduplicate")      = true,
            py::arg("filter_checks")    = true,
            py::arg("exploration")      = Exploration::ITERATIVE,
            py::arg("branching")        = Branching::BINARY,
            py::arg("lower_bounding")   = LowerBounding::BACKTRACK,
            py::arg("positioning")      = Positioning::ONLINE,
            py::arg("split_selection")  = SplitSelection::ALL,
            py::arg("split_scoring")    = SplitScoring::VARIANCE,
            py::arg("random_seed")      = 42,
            py::call_guard<py::gil_scoped_release>())
        .def("query", &LDTree::query,
            py::arg("cost"),
            py::arg("check_domain") = false,
            py::call_guard<py::gil_scoped_release>())
        .def("pprint", [](const LDTree& self, bool tightDisplay) {
                self.pprint(tightDisplay, &std::cout);
            },
            py::arg("tight_display") = false,
            py::call_guard<py::gil_scoped_release>())
        .def("flatten", &LDTree::flatten,
            py::arg("filepath"),
            py::arg("doc") = "",
            py::arg("benchmark_mode") = false,
            py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("domain", &LDTree::domain, py::return_value_policy::reference_internal)
        .def_property_readonly("voronoi", &LDTree::voronoi, py::return_value_policy::reference_internal)
        .def_property_readonly("tree", &LDTree::tree, py::return_value_policy::reference_internal)
        .def_property_readonly("stats", &LDTree::stats, py::return_value_policy::reference_internal)
        ;
}

void bind_module_io(py::module_ &m) {
    m.def("read_points", &readPoints, py::arg("filepath"));
    m.def("write_points", &writePoints, py::arg("filepath"), py::arg("points"));
    m.def("read_domain", &readDomain, py::arg("filepath"));
    m.def("write_domain", &writeDomain, py::arg("filepath"), py::arg("domain"));
    m.def("read_queries", &readQueries, py::arg("filepath"));
    m.def("write_queries", &writeQueries, py::arg("filepath"), py::arg("queries"));
}
