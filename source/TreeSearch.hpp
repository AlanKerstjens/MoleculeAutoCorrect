#ifndef _TREE_SEARCH_HPP_
#define _TREE_SEARCH_HPP_

#include <boost/graph/adjacency_list.hpp>
#include <boost/dynamic_bitset.hpp>
#include <concepts>
#include <random>
#include <queue>
#include <cmath>

template<class T>
concept Hashable = requires(T x) {
  { std::hash<T>{}(x) } -> std::convertible_to<std::size_t>;
};

template<class Vertex>
requires std::default_initializable<Vertex> && Hashable<Vertex>
class TreeSearch {
public:
  typedef boost::adjacency_list<
    boost::vecS, boost::vecS, boost::bidirectionalS, Vertex
  > Tree;

  // Function returning a pair, where first is a Vertex to expand from the Tree.
  // The policy may signal the end of the search by returning second false.
  typedef std::function<
    std::pair<typename Tree::vertex_descriptor, bool>(const TreeSearch&)
  > SelectionPolicy;

  // Function returning a new Vertex from an existing Vertex. The function is
  // allowed to fail and return a std::nullopt. Note that the call may modify
  // the state of the functor and/or Vertex. This allows the user to control the
  // output of future calls of Expansion.
  typedef std::function<
    std::optional<Vertex>(Vertex&)
  > Expansion;

  // Function returning a reward associated with a Vertex.
  // Required by certain expansion policies.
  typedef std::function<
    double(const Vertex&)
  > RewardFunction;

  // Function signalling whether the search should end or not.
  typedef std::function<
    bool(const TreeSearch&)
  > TerminationPolicy;

  // Map from Vertex hash to Tree::vertex_descriptor. Useful to check if a
  // Vertex resulting from Expansion has been generated before.
  typedef std::unordered_map<std::size_t, typename Tree::vertex_descriptor>
    VertexHashMap;

private:
  Tree tree;
  Expansion expansion;
  RewardFunction reward_function;
  std::hash<Vertex> vertex_hasher;
  VertexHashMap vertex_hash_map;
  std::vector<std::size_t> vertex_depths;
  std::vector<std::size_t> vertex_visits;
  std::vector<double> vertex_rewards;
  std::vector<bool> vertex_expandable;
  std::size_t max_n_vertices = 0, max_depth = 0;
  bool forward_depth = false;

public:
  TreeSearch(
    const Vertex& start_vertex,
    const Expansion& expansion,
    const RewardFunction& reward_function = nullptr,
    std::size_t max_n_vertices = 10000, std::size_t max_depth = 20,
    bool forward_depth = false) :
    expansion(expansion),
    reward_function(reward_function),
    max_n_vertices(max_n_vertices), max_depth(max_depth),
    forward_depth(forward_depth) {
    AddVertex(start_vertex, 0);
  };

private:
  std::pair<typename Tree::vertex_descriptor, bool> AddVertex(
    const Vertex& vertex,
    std::size_t depth) {
    std::size_t vertex_hash = vertex_hasher(vertex);
    typename VertexHashMap::const_iterator it = vertex_hash_map.find(vertex_hash);
    if (it == vertex_hash_map.cend()) {
      typename Tree::vertex_descriptor v = boost::add_vertex(vertex, tree);
      vertex_hash_map.emplace(vertex_hash, v);
      vertex_depths.push_back(depth);
      vertex_visits.push_back(0);
      vertex_rewards.push_back(0);
      vertex_expandable.push_back(depth < max_depth);
      return {v, true};
    } else {
      return {it->second, false};
    };
  };

  std::pair<typename Tree::vertex_descriptor, bool> Expand(
    Tree::vertex_descriptor v1) {
    if (!vertex_expandable[v1]) {
      return {v1, false};
    };
    std::optional<Vertex> vertex2 = expansion(tree[v1]);
    if (!vertex2) {
      vertex_expandable[v1] = false;
      return {v1, false};
    };
    auto [v2, is_new_vertex] = AddVertex(
      std::move(*vertex2), vertex_depths[v1] + 1);
    if (is_new_vertex) {
      boost::add_edge(v1, v2, tree);
      if (reward_function) {
        boost::dynamic_bitset<> rewarded_vertices (boost::num_vertices(tree));
        Backpropagate(v2, reward_function(tree[v2]), rewarded_vertices);
      };
    } else {
      // Checking whether the edge already exists or not is necessary since
      // different expansions of v1 could yield the same v2.
      auto [edge, edge_exists] = boost::edge(v1, v2, tree);
      if (!edge_exists) {
        boost::add_edge(v1, v2, tree);
        if (forward_depth) {
          ForwardDepth(v2, vertex_depths[v1] + 1);
        };
        vertex_visits[v1] += vertex_visits[v2];
        vertex_rewards[v1] += vertex_rewards[v2];
      };
    };
    return {v2, is_new_vertex};
  };

  void ForwardDepth(
    Tree::vertex_descriptor v,
    std::size_t new_depth) {
    std::size_t old_depth = vertex_depths[v];
    if (new_depth >= old_depth) {
      return;
    };
    for (typename Tree::edge_descriptor e :
      boost::make_iterator_range(boost::out_edges(v, tree))) {
      ForwardDepth(boost::target(e, tree), new_depth + 1);
    };
  };

  void Backpropagate(
    Tree::vertex_descriptor v,
    double reward,
    boost::dynamic_bitset<>& rewarded) {
    if (rewarded[v]) {
      return;
    };
    ++vertex_visits[v];
    vertex_rewards[v] += reward;
    rewarded.set(v);
    for (typename Tree::edge_descriptor e :
      boost::make_iterator_range(boost::in_edges(v, tree))) {
      Backpropagate(boost::source(e, tree), reward, rewarded);
    };
  };

public:
  std::size_t Search(
    const SelectionPolicy& selection_policy,
    const TerminationPolicy& terminator) {
    while (boost::num_vertices(tree) < max_n_vertices) {
      auto [parent, selected] = selection_policy(*this);
      if (!selected) {
        break;
      };
      auto [child, is_new_vertex] = Expand(parent);
      if (is_new_vertex && terminator(*this)) {
        break;
      };
    };
    return boost::num_vertices(tree);
  };

  std::vector<typename Tree::vertex_descriptor> ShortestPathToRoot(
    Tree::vertex_descriptor v) const {
    // There may be more than one equally short path to the root vertex.
    // Only one of them is returned here.
    std::size_t depth = vertex_depths[v];
    std::vector<typename Tree::vertex_descriptor> path;
    path.reserve(depth + 1);
    path.push_back(v);
    while (depth) {
      typename Tree::vertex_descriptor shallowest_parent;
      std::size_t shallowest_depth = depth;
      for (typename Tree::edge_descriptor e :
        boost::make_iterator_range(boost::in_edges(v, tree))) {
        typename Tree::vertex_descriptor parent = boost::source(e, tree);
        if (vertex_depths[parent] < shallowest_depth) {
          shallowest_parent = parent;
          shallowest_depth = vertex_depths[parent];
        };
      };
      path.push_back(shallowest_parent);
      v = shallowest_parent;
      depth = shallowest_depth;
    };
    return path;
  };

  std::vector<std::pair<typename Tree::vertex_descriptor, double>> TopVertices(
    const std::function<double(const Vertex&)>& scoring_function,
    std::size_t n) const {
    std::size_t n_vertices = boost::num_vertices(tree);
    std::vector<std::pair<typename Tree::vertex_descriptor, double>> top_vertices;
    top_vertices.reserve(n_vertices);
    for (typename Tree::vertex_descriptor v : boost::make_iterator_range(
      boost::vertices(tree))) {
      top_vertices.emplace_back(v, scoring_function(tree[v]));
    };
    std::sort(top_vertices.begin(), top_vertices.end(),
      [] (const auto& p1, const auto& p2) {
        if (p1.second != p2.second) {
          return p1.second > p2.second;
        };
        return p1.first < p2.first;
      }
    );
    std::size_t m = n < n_vertices ? n : n_vertices;
    return std::vector<std::pair<typename Tree::vertex_descriptor, double>> (
      top_vertices.begin(), top_vertices.begin() + m);
  };

  const Tree& GetTree() const {
    return tree;
  };

  const Vertex& GetVertex(typename Tree::vertex_descriptor v) const {
    return tree[v];
  };

  std::size_t GetVertexDepth(typename Tree::vertex_descriptor v) const {
    return vertex_depths[v];
  };

  Tree::vertex_descriptor GetVertexDescriptor(const Vertex& vertex) const {
    return vertex_hash_map.at(vertex_hasher(vertex));
  };

  Tree::vertex_descriptor GetLastVertexDescriptor() const {
    return boost::num_vertices(tree) - 1;
  };

  const Vertex& GetLastVertex() const {
    return tree[GetLastVertexDescriptor()];
  };

  const std::vector<std::size_t>& GetVertexDepths() const {
    return vertex_depths;
  };

  const std::vector<std::size_t>& GetVertexVisits() const {
    return vertex_visits;
  };

  const std::vector<double>& GetVertexRewards() const {
    return vertex_rewards;
  };

  const std::vector<bool>& GetExpandableVerticesMask() const {
    return vertex_expandable;
  };

  std::size_t GetMaxNVertices() const {
    return max_n_vertices;
  };

  std::size_t GetMaxDepth() const {
    return max_depth;
  };
};

template<class Vertex>
class GreedyPolicy {
  typedef TreeSearch<Vertex>::Tree Tree;

  // Negative priority values are reserved for unexpandable vertices.
  typedef std::function<
    double(const TreeSearch<Vertex>&, typename Tree::vertex_descriptor)
  > PriorityFunction;

  typedef std::pair<typename Tree::vertex_descriptor, double> QV; // Q'd Vertex
  struct QVComp {
    bool operator()(const QV& qv1, const QV& qv2) const {
      return qv1.second < qv2.second;
    };
  };

  PriorityFunction priority_function;
  std::priority_queue<QV, std::vector<QV>, QVComp> queue;
  std::size_t n_known_vertices = 0;

public:
  GreedyPolicy(const PriorityFunction& priority_function) :
    priority_function(priority_function) {};

  std::pair<typename Tree::vertex_descriptor, bool> operator()(
    const TreeSearch<Vertex>& tree_search) {
    // Calculate the priority of vertices that were added to the Tree since
    // the last operator() call.
    const Tree& tree = tree_search.GetTree();
    std::size_t n_vertices = boost::num_vertices(tree);
    const auto& vertex_expandable = tree_search.GetExpandableVerticesMask();
    for (std::size_t v = n_known_vertices; v < n_vertices; ++v) {
      if (vertex_expandable[v]) {
        double priority = priority_function(tree_search, v);
        if (priority >= 0.0) {
          queue.emplace(v, priority);
        };
      };
    };
    n_known_vertices = n_vertices;
    if (queue.empty()) {
      return {0, false};
    };
    // Return the Vertex with the highest priority. Some queued vertices may no 
    // longer be expandable. If that's the case remove them from the queue.
    auto [v, priority] = queue.top();
    while (!vertex_expandable[v]) {
      queue.pop();
      if (queue.empty()) {
        return {0, false};
      };
      std::tie(v, priority) = queue.top();
    };
    return {v, true};
  };
};

template <class Vertex>
class UpperConfidenceTree {
  typedef TreeSearch<Vertex>::Tree Tree;
  double c = 1.0;

private:
  double UCB1(
    const TreeSearch<Vertex>& tree_search,
    typename Tree::vertex_descriptor child,
    typename Tree::vertex_descriptor parent) const {
    const std::vector<double>& rewards = tree_search.GetVertexRewards();
    const std::vector<std::size_t>& visits = tree_search.GetVertexVisits();
    double xplt = rewards[child] / visits[child];
    double xplr = sqrt(log(visits[parent])/visits[child]);
    return xplt + c * xplr;
  };

public:
  UpperConfidenceTree(double c = sqrt(2)) : c(c) {};

  std::pair<typename Tree::vertex_descriptor, bool> operator()(
    const TreeSearch<Vertex>& tree_search) const {
    const Tree& tree = tree_search.GetTree();
    const auto& vertex_expandable = tree_search.GetExpandableVerticesMask();
    typename Tree::vertex_descriptor v = 0; // Root
    while (!vertex_expandable[v]) {
      if (!boost::out_degree(v, tree)) {
        break;
      };
      typename Tree::vertex_descriptor best_child;
      double best_ucb = -std::numeric_limits<double>::infinity();
      for (typename Tree::edge_descriptor e : boost::make_iterator_range(
        boost::out_edges(v, tree))) {
        typename Tree::vertex_descriptor child = boost::target(e, tree);
        double ucb = UCB1(tree_search, child, v);
        if (ucb > best_ucb) {
          best_child = child;
          best_ucb = ucb;
        };
      };
      v = best_child;
    };
    return {v, vertex_expandable[v]};
  };
};

template <class Vertex>
class RandomPolicy {
  typedef TreeSearch<Vertex>::Tree Tree;
  std::mt19937 prng;

public:
  RandomPolicy() {
    std::random_device rd;
    prng.seed(rd());
  };

  RandomPolicy(std::mt19937::result_type seed) {
    prng.seed(seed);
  };

  std::pair<typename Tree::vertex_descriptor, bool> operator()(
    const TreeSearch<Vertex>& tree_search) {
    const Tree& tree = tree_search.GetTree();
    const auto& vertex_expandable = tree_search.GetExpandableVerticesMask();
    std::size_t n_vertices = boost::num_vertices(tree);
    if (n_vertices < 2) {
      return {0, n_vertices ? vertex_expandable[0] : false};
    };
    std::uniform_int_distribution<std::size_t> distribution (0, n_vertices - 1);
    std::size_t vertex_idx = distribution(prng);
    for (std::size_t i = 0; i < n_vertices; ++i) {
      if (vertex_expandable[vertex_idx]) {
        return {vertex_idx, true};
      };
      if (++vertex_idx >= n_vertices) {
        vertex_idx = 0;
      };
    };
    return {0, false};
  };
};

#endif // !_TREE_SEARCH_HPP_
