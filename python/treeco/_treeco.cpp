#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>

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

PYBIND11_MODULE(_treeco, m) {
    m.doc() = "Treeco - Linear Decision Tree Policies for Combinatorial Optimization";

    // Dynprog enums
    py::enum_<Exploration>(m, "Exploration")
        .value("greedy", Exploration::GREEDY)
        .value("iterative", Exploration::ITERATIVE)
        .value("exhaustive", Exploration::EXHAUSTIVE)
        .export_values();

    py::enum_<Branching>(m, "Branching")
        .value("ternary", Branching::TERNARY)
        .value("binary", Branching::BINARY)
        .export_values();

    py::enum_<LowerBounding>(m, "LowerBounding")
        .value("fixed", LowerBounding::FIXED)
        .value("backtrack", LowerBounding::BACKTRACK)
        .export_values();

    py::enum_<Positioning>(m, "Positioning")
        .value("online", Positioning::ONLINE)
        .value("precompute", Positioning::PRECOMPUTE)
        .export_values();

    py::enum_<SplitSelection>(m, "SplitSelection")
        .value("all", SplitSelection::ALL)
        .value("sampling", SplitSelection::SAMPLING)
        .export_values();

    py::enum_<SplitScoring>(m, "SplitScoring")
        .value("variance", SplitScoring::VARIANCE)
        .value("entropy", SplitScoring::ENTROPY)
        .value("minmax", SplitScoring::MINMAX)
        .value("none", SplitScoring::NONE)
        .value("random", SplitScoring::RANDOM)
        .export_values();

    py::enum_<Relation>(m, "Relation")
        .value("lt", Relation::LT)
        .value("le", Relation::LE)
        .value("eq", Relation::EQ)
        .value("ge", Relation::GE)
        .value("gt", Relation::GT)
        .value("rt", Relation::RT)
        .value("rf", Relation::RF)
        .export_values();

    // LDTreeStats
    py::class_<LDTreeStats>(m, "LDTreeStats")
        .def(py::init<>())
        .def_readwrite("build_time", &LDTreeStats::buildTime);

    // LDTree
    py::class_<LDTree>(m, "LDTree")
        .def(py::init<const std::string&, const std::string&>(),
             py::arg("file_points"), py::arg("file_domain") = "")
        .def(py::init<const std::vector<BinaryVector>&, const Domain&>(),
             py::arg("points"), py::arg("domain") = Domain())
        .def("build", [](
                LDTree&         self,
                bool            verbose,
                double          logInterval,
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
        .def("domain", &LDTree::domain, py::return_value_policy::reference_internal)
        .def("voronoi", &LDTree::voronoi, py::return_value_policy::reference_internal)
        .def("tree", &LDTree::tree, py::return_value_policy::reference_internal)
        .def("stats", &LDTree::stats, py::return_value_policy::reference_internal);

    // Problem (abstract class)
    py::class_<Problem, std::shared_ptr<Problem>>(m, "Problem")
        .def("dimension", &Problem::dimension)
        .def("get_feasible_set", &Problem::getFeasibleSet)
        .def("get_cost_domain", &Problem::getCostDomain)
        .def("sample_cost", static_cast<RealVector (Problem::*)() const>(&Problem::sampleCost))
        .def("sample_cost", static_cast<RealVector (Problem::*)(Index) const>(&Problem::sampleCost), py::arg("seed"))
        ;

    // Explicit
    py::class_<Explicit, Problem, std::shared_ptr<Explicit>>(m, "Explicit")
        .def(py::init<const std::vector<BinaryVector>&>(), py::arg("feasible_set"))
        .def("get_feasible_set", &Explicit::getFeasibleSet, py::call_guard<py::gil_scoped_release>());

    // Knapsack
    py::class_<Knapsack, Problem, std::shared_ptr<Knapsack>>(m, "Knapsack")
        .def(py::init<Index>(), py::arg("num_items"))
        .def(py::init<Index, const RealVector&, double>(), py::arg("num_items"), py::arg("weights"), py::arg("capacity"))
        .def("num_items", &Knapsack::numItems)
        .def("get_feasible_set", &Knapsack::getFeasibleSet, py::call_guard<py::gil_scoped_release>())
        .def("get_cost_domain", &Knapsack::getCostDomain)
        .def("sample_cost", static_cast<RealVector (Problem::*)() const>(&Problem::sampleCost))
        .def("sample_cost", static_cast<RealVector (Knapsack::*)(Index) const>(&Knapsack::sampleCost), py::arg("seed"))
        .def("weights", &Knapsack::weights)
        .def("capacity", &Knapsack::capacity);

    // Maxcut
    py::class_<Maxcut, Problem, std::shared_ptr<Maxcut>>(m, "Maxcut")
        .def(py::init<Index>(), py::arg("num_vertices"))
        .def("num_vertices", &Maxcut::numVertices)
        .def("get_feasible_set", &Maxcut::getFeasibleSet, py::call_guard<py::gil_scoped_release>())
        .def("get_cost_domain", &Maxcut::getCostDomain)
        .def("sample_cost", static_cast<RealVector (Problem::*)() const>(&Problem::sampleCost))
        .def("sample_cost", static_cast<RealVector (Maxcut::*)(Index) const>(&Maxcut::sampleCost), py::arg("seed"))
        .def("edge_to_index", &Maxcut::edgeToIndex)
        .def("index_to_edge", &Maxcut::indexToEdge);

    // Tsp
    py::class_<Tsp, Problem, std::shared_ptr<Tsp>>(m, "Tsp")
        .def(py::init<Index>(), py::arg("num_cities"))
        .def("num_cities", &Tsp::numCities)
        .def("get_feasible_set", &Tsp::getFeasibleSet, py::call_guard<py::gil_scoped_release>())
        .def("get_cost_domain", &Tsp::getCostDomain)
        .def("sample_cost", static_cast<RealVector (Problem::*)() const>(&Problem::sampleCost))
        .def("sample_cost", static_cast<RealVector (Tsp::*)(Index) const>(&Tsp::sampleCost), py::arg("seed"))
        .def("edge_to_index", &Tsp::edgeToIndex)
        .def("index_to_edge", &Tsp::indexToEdge);

    // IO
    m.def("read_points", &readPoints, py::arg("filepath"));
    m.def("write_points", &writePoints, py::arg("filepath"), py::arg("points"));
    m.def("read_domain", &readDomain, py::arg("filepath"));
    m.def("write_domain", &writeDomain, py::arg("filepath"), py::arg("domain"));
    m.def("read_queries", &readQueries, py::arg("filepath"));
    m.def("write_queries", &writeQueries, py::arg("filepath"), py::arg("queries"));
}
