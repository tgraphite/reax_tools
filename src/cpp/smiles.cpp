#include "smiles.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_set>

namespace reax_smiles {
namespace {

// ---------------------------------------------------------------------------
// small utilities
// ---------------------------------------------------------------------------
std::uint64_t fnv1a64(const std::string& text) {
    std::uint64_t hash = 1469598103934665603ULL;
    for (unsigned char ch : text) {
        hash ^= static_cast<std::uint64_t>(ch);
        hash *= 1099511628211ULL;
    }
    return hash;
}

std::string normalize_symbol(const std::string& element) {
    if (element.empty()) return "*";
    if (element == "X" || element == "x") return "*";
    std::string out = element;
    out[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(out[0])));
    for (std::size_t i = 1; i < out.size(); ++i) {
        out[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(out[i])));
    }
    return out;
}

struct HGraph {
    int n = 0;
    std::vector<std::string> element;
    std::vector<int> charge;
    std::vector<int> hcount;
    std::vector<std::vector<int>> adj;
    std::vector<int> degree;
};

struct Signature {
    int self = 0;
    std::vector<int> neigh;

    bool operator<(const Signature& other) const {
        if (self != other.self) return self < other.self;
        return neigh < other.neigh;
    }
    bool operator==(const Signature& other) const {
        return self == other.self && neigh == other.neigh;
    }
};

std::vector<int> initial_colors(const HGraph& graph) {
    std::vector<std::string> symbols = graph.element;
    std::sort(symbols.begin(), symbols.end());
    symbols.erase(std::unique(symbols.begin(), symbols.end()), symbols.end());
    std::map<std::string, int> element_rank;
    for (std::size_t i = 0; i < symbols.size(); ++i) {
        element_rank[symbols[i]] = static_cast<int>(i);
    }

    // Assign initial colour ids by *sorted* class key, never by atom
    // encounter order.  Otherwise the same molecule with a different input
    // atom order receives different initial colours and hence a different
    // canonical traversal.
    std::set<std::tuple<int, int, int>> keys;
    for (int v = 0; v < graph.n; ++v) {
        keys.emplace(element_rank[graph.element[v]], graph.hcount[v], graph.degree[v]);
    }
    std::map<std::tuple<int, int, int>, int> class_of;
    int next = 0;
    for (const auto& key : keys) class_of[key] = next++;

    std::vector<int> colors(graph.n, 0);
    for (int v = 0; v < graph.n; ++v) {
        auto key = std::make_tuple(element_rank[graph.element[v]], graph.hcount[v], graph.degree[v]);
        colors[v] = class_of[key];
    }
    return colors;
}

// Colour refinement. Returns colours after at most max_rounds; ok is false if
// the partition did not stabilise within the budget.
std::vector<int> refine(const HGraph& graph,
                        std::vector<int> colors,
                        std::size_t max_rounds,
                        bool& ok) {
    ok = true;
    if (graph.n == 0) return colors;

    for (std::size_t round = 0; round < max_rounds; ++round) {
        int old_k = *std::max_element(colors.begin(), colors.end()) + 1;

        std::vector<Signature> signatures(graph.n);
        for (int v = 0; v < graph.n; ++v) {
            signatures[v].self = colors[v];
            signatures[v].neigh.reserve(graph.adj[v].size());
            for (int u : graph.adj[v]) signatures[v].neigh.push_back(colors[u]);
            std::sort(signatures[v].neigh.begin(), signatures[v].neigh.end());
        }

        std::vector<int> order(graph.n);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            return signatures[a] < signatures[b];
        });

        std::vector<int> next(graph.n, 0);
        int new_k = 1;
        next[order[0]] = 0;
        for (std::size_t i = 1; i < order.size(); ++i) {
            if (!(signatures[order[i]] == signatures[order[i - 1]])) ++new_k;
            next[order[i]] = new_k - 1;
        }

        bool stable = (new_k == old_k);
        colors.swap(next);
        if (stable) return colors;
    }

    ok = false;
    return colors;
}

std::string atom_token(const HGraph& graph, int v) {
    const std::string& symbol = graph.element[v];
    std::string token = "[";
    token += symbol;
    if (graph.charge[v] > 0) {
        token += "+";
        if (graph.charge[v] > 1) token += std::to_string(graph.charge[v]);
    } else if (graph.charge[v] < 0) {
        token += "-";
        if (graph.charge[v] < -1) token += std::to_string(-graph.charge[v]);
    }
    token += "]";
    return token;
}

std::string ring_symbol(int number) {
    if (number <= 9) return std::to_string(number);
    return "%" + std::to_string(number);
}

// Write SMILES for one connected component given a canonical vertex order.
std::string write_component(const HGraph& graph, const std::vector<int>& order) {
    const int n = graph.n;
    std::vector<int> rank(n, 0);
    for (int i = 0; i < n; ++i) rank[order[i]] = i;

    std::vector<std::vector<int>> adj_by_rank(n);
    for (int v = 0; v < n; ++v) {
        adj_by_rank[v] = graph.adj[v];
        std::sort(adj_by_rank[v].begin(), adj_by_rank[v].end(),
                  [&](int a, int b) { return rank[a] < rank[b]; });
    }

    std::vector<int> parent(n, -1);
    std::vector<char> visited(n, 0);
    std::vector<std::vector<int>> children(n);
    std::vector<std::pair<int, std::size_t>> stack;
    int root = order[0];
    visited[root] = 1;
    stack.push_back({root, 0});
    while (!stack.empty()) {
        int v = stack.back().first;
        std::size_t& next = stack.back().second;
        if (next < adj_by_rank[v].size()) {
            int u = adj_by_rank[v][next++];
            if (!visited[u]) {
                visited[u] = 1;
                parent[u] = v;
                children[v].push_back(u);
                stack.push_back({u, 0});
            }
        } else {
            stack.pop_back();
        }
    }

    // Non-tree edges are ring closures. Record each edge once (lower rank first).
    std::vector<std::pair<int, int>> ring_edges;
    for (int v = 0; v < n; ++v) {
        for (int u : graph.adj[v]) {
            if (rank[v] < rank[u] && parent[u] != v && parent[v] != u) {
                ring_edges.emplace_back(v, u);
            }
        }
    }
    std::sort(ring_edges.begin(), ring_edges.end(), [&](const auto& a, const auto& b) {
        if (rank[a.first] != rank[b.first]) return rank[a.first] < rank[b.first];
        return rank[a.second] < rank[b.second];
    });

    std::vector<std::vector<int>> ring_at(n);
    for (std::size_t i = 0; i < ring_edges.size(); ++i) {
        int number = static_cast<int>(i) + 1;
        ring_at[ring_edges[i].first].push_back(number);
        ring_at[ring_edges[i].second].push_back(number);
    }
    for (int v = 0; v < n; ++v) {
        std::sort(ring_at[v].begin(), ring_at[v].end());
    }

    std::string out;
    std::function<void(int)> emit = [&](int v) {
        out += atom_token(graph, v);
        for (int number : ring_at[v]) out += ring_symbol(number);
        for (int i = 0; i < graph.hcount[v]; ++i) out += "([H])";
        // In SMILES a branch written after an inline child would attach to
        // the last atom of that child's subtree.  Emit all but the final
        // child as branches of v, then continue with the final child.
        for (std::size_t i = 0; i + 1 < children[v].size(); ++i) {
            out += "(";
            emit(children[v][i]);
            out += ")";
        }
        if (!children[v].empty()) emit(children[v].back());
    };
    emit(root);
    return out;
}

struct SearchContext {
    std::size_t leaves = 0;
    bool ok = true;
    std::string error;
};

std::string individualise(const HGraph& graph,
                          std::vector<int> colors,
                          const Options& options,
                          SearchContext& context) {
    if (!context.ok) return "";

    int k = *std::max_element(colors.begin(), colors.end()) + 1;
    if (k == graph.n) {
        if (++context.leaves > options.max_ir_leaves) {
            context.ok = false;
            context.error = "individualisation budget exceeded";
            return "";
        }
        std::vector<int> order(graph.n, 0);
        for (int v = 0; v < graph.n; ++v) order[colors[v]] = v;
        return write_component(graph, order);
    }

    std::map<int, std::vector<int>> cells;
    for (int v = 0; v < graph.n; ++v) cells[colors[v]].push_back(v);

    int chosen_color = -1;
    std::size_t chosen_size = std::numeric_limits<std::size_t>::max();
    for (const auto& [color, members] : cells) {
        if (members.size() <= 1) continue;
        if (members.size() < chosen_size ||
            (members.size() == chosen_size && (chosen_color < 0 || color < chosen_color))) {
            chosen_color = color;
            chosen_size = members.size();
        }
    }
    if (chosen_color < 0) {
        // Should not happen: k < n implies a non-singleton cell exists.
        context.ok = false;
        context.error = "internal partition error";
        return "";
    }
    if (chosen_size > options.max_cell) {
        context.ok = false;
        context.error = "individualisation cell too large";
        return "";
    }

    std::string best;
    for (int v : cells[chosen_color]) {
        std::vector<int> next_colors = colors;
        next_colors[v] = k;
        bool refined_ok = true;
        next_colors = refine(graph, next_colors, options.max_refine_rounds, refined_ok);
        if (!refined_ok) {
            context.ok = false;
            context.error = "colour refinement budget exceeded";
            return "";
        }
        std::string candidate = individualise(graph, next_colors, options, context);
        if (!context.ok) return "";
        if (best.empty() || candidate < best) best = candidate;
    }
    return best;
}

// Split the H-suppressed graph into connected components and canonicalise each.
std::vector<int> component_of(const HGraph& graph) {
    std::vector<int> component(graph.n, -1);
    int next = 0;
    for (int start = 0; start < graph.n; ++start) {
        if (component[start] != -1) continue;
        std::vector<int> stack = {start};
        component[start] = next;
        while (!stack.empty()) {
            int v = stack.back();
            stack.pop_back();
            for (int u : graph.adj[v]) {
                if (component[u] == -1) {
                    component[u] = next;
                    stack.push_back(u);
                }
            }
        }
        ++next;
    }
    return component;
}

HGraph subgraph(const HGraph& graph, const std::vector<int>& vertices) {
    HGraph out;
    out.n = static_cast<int>(vertices.size());
    std::map<int, int> local;
    for (int i = 0; i < out.n; ++i) local[vertices[i]] = i;
    out.element.resize(out.n);
    out.charge.resize(out.n);
    out.hcount.resize(out.n);
    out.adj.resize(out.n);
    out.degree.resize(out.n);
    for (int i = 0; i < out.n; ++i) {
        int v = vertices[i];
        out.element[i] = graph.element[v];
        out.charge[i] = graph.charge[v];
        out.hcount[i] = graph.hcount[v];
    }
    for (int i = 0; i < out.n; ++i) {
        int v = vertices[i];
        for (int u : graph.adj[v]) {
            auto it = local.find(u);
            if (it != local.end()) out.adj[i].push_back(it->second);
        }
        std::sort(out.adj[i].begin(), out.adj[i].end());
        out.adj[i].erase(std::unique(out.adj[i].begin(), out.adj[i].end()), out.adj[i].end());
        out.degree[i] = static_cast<int>(out.adj[i].size());
    }
    return out;
}

bool build_h_suppressed(const std::vector<Atom>& atoms,
                        const std::vector<Bond>& bonds,
                        HGraph& graph,
                        std::string& error) {
    const int n_full = static_cast<int>(atoms.size());
    if (n_full <= 0) {
        error = "empty atom list";
        return false;
    }

    std::vector<std::vector<int>> full_adj(n_full);
    std::set<std::pair<int, int>> seen;
    for (const auto& bond : bonds) {
        if (bond.i < 0 || bond.j < 0 || bond.i >= n_full || bond.j >= n_full) {
            error = "bond index out of range";
            return false;
        }
        if (bond.i == bond.j) {
            error = "self bond";
            return false;
        }
        int a = std::min(bond.i, bond.j);
        int b = std::max(bond.i, bond.j);
        if (!seen.emplace(a, b).second) continue;
        full_adj[a].push_back(b);
        full_adj[b].push_back(a);
    }

    std::vector<std::string> symbol(n_full);
    for (int i = 0; i < n_full; ++i) symbol[i] = normalize_symbol(atoms[i].element);

    std::vector<char> keep(n_full, 1);
    std::vector<int> hcount(n_full, 0);
    for (int i = 0; i < n_full; ++i) {
        if (symbol[i] != "H") continue;
        if (full_adj[i].size() == 1) {
            int neighbour = full_adj[i][0];
            if (symbol[neighbour] != "H") {
                keep[i] = 0;
                hcount[neighbour] += 1;
            }
        }
    }

    std::vector<int> local(n_full, -1);
    std::vector<int> vertices;
    for (int i = 0; i < n_full; ++i) {
        if (!keep[i]) continue;
        local[i] = static_cast<int>(vertices.size());
        vertices.push_back(i);
    }
    if (vertices.empty()) {
        error = "no atoms after hydrogen suppression";
        return false;
    }

    graph.n = static_cast<int>(vertices.size());
    graph.element.resize(graph.n);
    graph.charge.resize(graph.n);
    graph.hcount.resize(graph.n);
    graph.adj.resize(graph.n);
    graph.degree.resize(graph.n);
    for (int i = 0; i < graph.n; ++i) {
        int old = vertices[i];
        graph.element[i] = symbol[old];
        graph.charge[i] = atoms[old].charge;
        graph.hcount[i] = hcount[old];
    }
    for (int i = 0; i < graph.n; ++i) {
        int old = vertices[i];
        for (int nb : full_adj[old]) {
            int j = local[nb];
            if (j >= 0) graph.adj[i].push_back(j);
        }
        std::sort(graph.adj[i].begin(), graph.adj[i].end());
        graph.adj[i].erase(std::unique(graph.adj[i].begin(), graph.adj[i].end()), graph.adj[i].end());
        graph.degree[i] = static_cast<int>(graph.adj[i].size());
    }
    return true;
}

Result make_fallback(const std::string& error) {
    Result result;
    result.fallback = true;
    result.error = error;
    return result;
}

}  // namespace

Result canonical_smiles(const std::vector<Atom>& atoms,
                        const std::vector<Bond>& bonds,
                        const Options& options) {
    Result result;
    if (atoms.empty()) {
        result.error = "empty atom list";
        return result;
    }
    if (atoms.size() > options.max_atoms) {
        return make_fallback("atom count exceeds max_atoms");
    }

    HGraph graph;
    std::string error;
    if (!build_h_suppressed(atoms, bonds, graph, error)) {
        result.error = error;
        return result;
    }

    std::vector<int> component = component_of(graph);
    int component_count = 0;
    for (int value : component) component_count = std::max(component_count, value + 1);

    std::vector<std::string> pieces;
    SearchContext context;
    for (int c = 0; c < component_count; ++c) {
        std::vector<int> vertices;
        for (int v = 0; v < graph.n; ++v) {
            if (component[v] == c) vertices.push_back(v);
        }
        HGraph part = subgraph(graph, vertices);
        std::vector<int> colors = initial_colors(part);
        bool ok = true;
        colors = refine(part, colors, options.max_refine_rounds, ok);
        if (!ok) {
            return make_fallback("colour refinement budget exceeded");
        }
        std::string piece = individualise(part, colors, options, context);
        if (!context.ok) {
            return make_fallback(context.error);
        }
        pieces.push_back(piece);
    }

    std::sort(pieces.begin(), pieces.end());
    std::string joined;
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        if (i > 0) joined += ".";
        joined += pieces[i];
    }
    if (joined.empty()) {
        return make_fallback("empty SMILES");
    }
    result.smiles = joined;
    return result;
}

std::uint64_t smiles_id(const std::string& smiles) {
    return fnv1a64("smiles|" + smiles);
}

std::uint64_t formula_id(const std::string& formula) {
    return fnv1a64("formula|" + formula);
}

}  // namespace reax_smiles
