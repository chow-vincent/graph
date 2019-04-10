#pragma once

// TODO:
// - find out all the package requirements needed to run these files
// - eliminate dependency on cme212 files

#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/counting_iterator.h>

#include <vector>
#include <set>
#include <cassert>

#include "utils/util.hpp"
// #include "CME212/Point.hpp"
#include "utils/point.hpp"

template <typename V, typename E>
class Graph {
  private:

    // forward declarations of underlying data
    struct NodeInternal;
    struct EdgeInternal;

  public:
    // forward declaration of Node and Edge
    class Node;
    using node_type = Node;

    class Edge;
    using edge_type = Edge;

    using NodeValue = V;
    using EdgeValue = E;

    using NodeID = std::size_t;
    using EdgeID = std::size_t;

    // aliases for underlying data structures
    using NodeVector = std::vector<NodeInternal>;
    using EdgeVector = std::vector<EdgeInternal>;
    using AdjEdgesSet = std::set<EdgeID>;

    // forward declarations of iterators
    class IncidentIterator;

    // constructor for empty graph
    Graph() {
    }

    // default destructor
    ~Graph() = default;

    // Node class. Proxy for NodeInternal.
    class Node : private totally_ordered<Node> {
      public:

        // constructor for invalid node
        Node() {
          node_id_ = -1;
          graph_ptr_ = nullptr;
        }

        Node(NodeID node_id, const Graph *graph_ptr)
          : node_id_(node_id), graph_ptr_(const_cast<Graph*>(graph_ptr)) {
        }

        // modifiable, non-const ref to position
        // precondition assumes nodes are valid for these fxns to be called
        Point &position() {
          return (graph_ptr_->nodes_)[node_id_].position;
        }

        const Point &position() const {
          return (graph_ptr_->nodes_)[node_id_].position;
        }

        // modifiable, non-const ref to value
        NodeValue &value() {
          return (graph_ptr_->nodes_)[node_id_].node_value;
        }

        const NodeValue &value() const {
          return (graph_ptr_->nodes_)[node_id_].node_value;
        }

        std::size_t degree() const {
          const auto &adj_edges = (graph_ptr_->nodes_)[node_id_].adjacent_edges;
          return adj_edges.size();
        }

        IncidentIterator edge_begin() const {
          const NodeInternal &node = (graph_ptr_->nodes_)[node_id_];
          return IncidentIterator((node.adjacent_edges).begin(), node_id_, graph_ptr_);
        }

        IncidentIterator edge_end() const {
          const NodeInternal &node = (graph_ptr_->nodes_)[node_id_];
          return IncidentIterator((node.adjacent_edges).end(), node_id_, graph_ptr_);
        }

        // comparators
        bool operator==(const Node &n) const {
          return (graph_ptr_ == n.get_graph_pointer()) &&
                  (node_id_ == n.get_node_id());
        }

        bool operator<(const Node &n) const {
          return node_id_ < n.get_node_id();
        }

        bool operator>(const Node &n) const {
          return node_id_ > n.get_node_id();
        }

        // ID, synonymous with index
        NodeID get_node_id() const {
          return node_id_;
        }

        // index, synonymous with ID
        NodeID index() const {
          return node_id_;
        }

        Graph *get_graph_pointer() const {
          return graph_ptr_;
        }

      private:
        friend class Graph;

        // ID of node proxy (index into NodeVector for this implementation)
        NodeID node_id_;

        // pointer to the graph to which this node belongs
        Graph *graph_ptr_;
    };

    // size of the graph
    std::size_t size() const {
      return nodes_.size();
    }

    std::size_t num_nodes() const {
      return size();
    }

    Node add_node(const Point &position, const NodeValue &node_value = NodeValue()) {
      nodes_.emplace_back(position, node_value);
      return Node(nodes_.size() - 1, this);
    }

    std::size_t remove_node(const Node &n) {

      // case where node not found
      if (!has_node(n)) {
        return 0;
      }

      const NodeID last_node_id = nodes_.size() - 1;
      const NodeID node_id = n.get_node_id(); // node to be deleted

      // remove edges incident to n
      AdjEdgesSet adj_edges_copy = nodes_[node_id].adjacent_edges;
      for (EdgeID edge_id : adj_edges_copy) {
        remove_edge(edge(edge_id));
      }

      update_adjacent_edges(nodes_.back().adjacent_edges, node_id, last_node_id);

      // swap and pop
      std::swap(nodes_[node_id], nodes_[last_node_id]);
      nodes_.pop_back();
      return 1;
    }

    // helper to update the ID of the node swapped in for the deleted one
    void update_adjacent_edges(AdjEdgesSet &adj_edges, const NodeID node_id, const NodeID last_node_id) {
      // if we are removing last node, do nothing
      if (node_id == last_node_id) {
        return;
      }

      for (auto it = adj_edges.begin(); it != adj_edges.end(); ++it) {
        EdgeInternal &edge = edges_[*it];
        if (edge.node1_id == last_node_id) {
          edge.node1_id = node_id;
        } else {
          edge.node2_id = node_id;
        }
      }
    }

    // TODO?:
    // NodeIterator remove_node(NodeIterator n_it) {
    // }

    bool has_node(const Node &n) const {
      return n.get_node_id() < nodes_.size() && n.get_graph_pointer() == this;
    }

    // precondition that id in range
    Node node(NodeID id) const {
      return Node(id, this);
    }

    // Edge class. Proxy for EdgeInternal.
    class Edge : private totally_ordered<Edge> {
      public:
        // constructor for invalid edge
        Edge() {
          edge_id_ = -1;
          node1_id_ = -1;
          node2_id_ = -1;
          graph_ptr_ = nullptr;
        }

        Edge(EdgeID edge_id, NodeID node1_id, NodeID node2_id, const Graph *graph_ptr)
          : edge_id_(edge_id), node1_id_(node1_id), node2_id_(node2_id),
            graph_ptr_(const_cast<Graph*>(graph_ptr)) {
        }

        Node node1() const {
          return Node(node1_id_, graph_ptr_);
        }

        Node node2() const {
          return Node(node2_id_, graph_ptr_);
        }

        EdgeValue &value() {
          return (graph_ptr_->edges_)[edge_id_].edge_value;
        }

        const EdgeValue &value() const {
          return (graph_ptr_->edges_)[edge_id_].edge_value;
        }

        // euclidean length of edge
        double length() const {
          const Point pos1 = (graph_ptr_->nodes_)[node1_id_].position;
          const Point pos2 = (graph_ptr_->nodes_)[node2_id_].position;
          return norm_2(pos1 - pos2);
        }

        bool operator==(const Edge &e) const {
          return (graph_ptr_ == e.get_graph_pointer()) &&
                  (edge_id_ == e.get_edge_id());
        }

        bool operator<(const Edge &e) const {
          return (edge_id_ < e.get_edge_id() || graph_ptr_ < e.get_graph_pointer());
        }

        EdgeID get_edge_id() const {
          return edge_id_;
        }

        Graph *get_graph_pointer() const {
          return graph_ptr_;
        }

      private:
        friend class graph;

        // ID of edge proxy (index into EdgeVector for this implementation)
        EdgeID edge_id_;

        // NOTE: not necessarily same order as stored in NodeInternal
        // the node spawning the incident iterator is always set to be node1
        NodeID node1_id_;
        NodeID node2_id_;

        // pointer to the graph to which this node belongs
        Graph *graph_ptr_;
    };

    std::size_t num_edges() const {
      return edges_.size();
    }

    Edge edge(EdgeID id) const {
      const EdgeInternal &edge = edges_[id];
      return Edge(id, edge.node1_id, edge.node2_id, this);
    }

    std::size_t remove_edge(const Edge &e) {

      // case where edge not found
      if (!has_edge(e.node1(), e.node2())) {
        return 0;
      }

      const EdgeID last_edge_id = edges_.size() - 1;
      const EdgeID edge_id = e.get_edge_id();

      // retrieve edge to be deleted
      const EdgeInternal &edge = edges_[edge_id];

      // delete edge_id from adj edges of both nodes
      NodeInternal &node_a = nodes_[edge.node1_id];
      NodeInternal &node_b = nodes_[edge.node2_id];
      (node_a.adjacent_edges).erase(edge_id);
      (node_b.adjacent_edges).erase(edge_id);

      // update the edge ID of the one to be swapped for deleted edge
      if (last_edge_id != edge_id) {
        EdgeInternal &last_edge = edges_[last_edge_id];
        AdjEdgesSet &adj_edges1 = nodes_[last_edge.node1_id].adjacent_edges;
        AdjEdgesSet &adj_edges2 = nodes_[last_edge.node2_id].adjacent_edges;
        adj_edges1.erase(last_edge_id);
        adj_edges1.insert(edge_id);
        adj_edges2.erase(last_edge_id);
        adj_edges2.insert(edge_id);
      }

      std::swap(edges_[edge_id], edges_[last_edge_id]);
      edges_.pop_back();
      return 1;
    }

    std::size_t remove_edge(const Node &a, const Node &b) {
      Edge edge = find_connecting_edge(a, b);
      return remove_edge(edge);
    }

    // TODO?:
    // EdgeIterator remove_edge(EdgeIterator e_it) {
    // }

    bool has_edge(const Node &a, const Node &b) const {
      Edge edge = find_connecting_edge(a, b);
      return edge.get_graph_pointer(); // nullptr if edge invalid
    }

    // helper function to find edge connecting two nodes
    Edge find_connecting_edge(const Node &a, const Node &b) const {

      // invalid edge if either node is invalid
      if (!has_node(a) || !has_node(b)) {
        return Edge();
      }

      // invalid edge if a and bare the same node
      if (a.get_node_id() == b.get_node_id()) {
        return Edge();
      }

      // iterate through incident edges to find the right one
      const NodeInternal &node = nodes_[a.get_node_id()];
      for (const EdgeID edge_id : node.adjacent_edges) {
        const EdgeInternal &edge = edges_[edge_id];
        if (edge.node1_id == b.get_node_id() || edge.node2_id == b.get_node_id()) {
          return Edge(edge_id, edge.node1_id, edge.node2_id, this);
        }
      }

      return Edge(); // invalid if could not find
    }

    Edge add_edge(const Node& a, const Node& b, const EdgeValue &edge_value = EdgeValue()) {
      Edge edge = find_connecting_edge(a, b);

      // return Edge if already found
      if (edge.get_graph_pointer()) {
        return edge;
      }

      const EdgeID new_edge_id = edges_.size();

      // add EdgeInternal to graph
      NodeID id_a = a.get_node_id();
      NodeID id_b = b.get_node_id();
      edges_.emplace_back(id_a, id_b, edge_value);

      // add to adjacent edges of Nodes a and b
      NodeInternal &node_a = nodes_[id_a];
      NodeInternal &node_b = nodes_[id_b];
      (node_a.adjacent_edges).insert(new_edge_id);
      (node_b.adjacent_edges).insert(new_edge_id);

      return Edge(new_edge_id, id_a, id_b, this);
    }

    void clear() {
      nodes_.clear();
      edges_.clear();
    }

    // functor for NodeIterator class
    struct GetNode {
      Graph *graph_ptr_;

      Node operator()(const NodeID id) {
        return Node(id, graph_ptr_);
      }

      GetNode(const Graph *graph_ptr) : graph_ptr_(const_cast<Graph*>(graph_ptr)) {
      }
    };

    // class for iterating over graph's nodes
    class NodeIterator : public thrust::transform_iterator<GetNode, thrust::counting_iterator<NodeID>, Node> {
      public:
        using super_t = thrust::transform_iterator<GetNode, thrust::counting_iterator<NodeID>, Node>;
        NodeIterator() {
        }

        NodeIterator(const super_t &ti) : super_t{ti} {
        }

      private:
        friend class Graph;
        NodeIterator(const thrust::counting_iterator<NodeID> &it, const Graph *g)
          : super_t{it, GetNode(g)} {
        }
    };

    NodeIterator node_begin() const {
      return NodeIterator(thrust::counting_iterator<NodeID>(0), this);
    }

    NodeIterator node_end() const {
      return NodeIterator(thrust::counting_iterator<NodeID>(num_nodes()), this);
    }

    class IncidentIterator {
      public:
        // aliases to use STL iterator_traits
        using value_type = Edge;
        using pointer    = Edge*;
        using reference  = Edge&;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;

        // constructor for invalid incident iterator
        IncidentIterator() {
        }

        IncidentIterator(typename AdjEdgesSet::iterator adj_edges_iterator,
                          const NodeID node_id, const Graph *graph_ptr)
          : adj_edges_iterator_(adj_edges_iterator), node_id_(node_id),
            graph_ptr_(const_cast<Graph*>(graph_ptr)) {
        }

        Edge operator*() const {
          const EdgeID edge_id = *adj_edges_iterator_;
          const EdgeInternal &edge = (graph_ptr_->edges_)[edge_id];

          // manipulate node1 and node2 so that node1 is node_id_
          if (node_id_ == edge.node1_id) {
            return Edge(edge_id, edge.node1_id, edge.node2_id, graph_ptr_);
          } else {
            return Edge(edge_id, edge.node2_id, edge.node1_id, graph_ptr_);
          }
        }

        IncidentIterator &operator++() {
          ++adj_edges_iterator_;
          return *this;
        }

        bool operator==(const IncidentIterator &iit) const {
          return (adj_edges_iterator_ == iit.get_iterator()) &&
                  (node_id_ == iit.get_node_id()) &&
                  (graph_ptr_ == iit.get_graph_pointer());
        }

        bool operator!=(const IncidentIterator &iit) const {
          return !(*this == iit);
        }

        typename AdjEdgesSet::iterator get_iterator() const {
          return adj_edges_iterator_;
        }

        NodeID get_node_id() const {
          return node_id_;
        }

        Graph *get_graph_pointer() const {
          return graph_ptr_;
        }

      private:
        friend class Graph;

        typename AdjEdgesSet::iterator adj_edges_iterator_;

        // ID of the node spawning this iterator (to ensure proper orientation)
        NodeID node_id_;

        Graph *graph_ptr_;
    };

    // functor for EdgeIterator class
    struct GetEdge {
      Graph *graph_ptr_;

      Edge operator()(const EdgeID id) {
        const EdgeInternal &edge = (graph_ptr_->edges_)[id];
        return Edge(id, edge.node1_id, edge.node2_id, graph_ptr_);
      }

      GetEdge(const Graph *graph_ptr) : graph_ptr_(const_cast<Graph*>(graph_ptr)) {
      }
    };

    // class for iterating over graph's edges
    class EdgeIterator : public thrust::transform_iterator<GetEdge, thrust::counting_iterator<EdgeID>, Edge> {
      public:
        using super_t = thrust::transform_iterator<GetEdge, thrust::counting_iterator<EdgeID>, Edge>;
        EdgeIterator() {
        }

        EdgeIterator(const super_t &ti) : super_t{ti} {
        }

      private:
        friend class Graph;
        EdgeIterator(const thrust::counting_iterator<NodeID> &it, const Graph *g)
          : super_t{it, GetEdge(g)} {
        }
    };

    EdgeIterator edge_begin() const {
      return EdgeIterator(thrust::counting_iterator<EdgeID>(0), this);
    }

    EdgeIterator edge_end() const {
      return EdgeIterator(thrust::counting_iterator<EdgeID>(num_edges()), this);
    }

  private:

    // true, underlying node
    struct NodeInternal {
      Point position;
      NodeValue node_value;
      AdjEdgesSet adjacent_edges;
      NodeInternal(const Point &position_, const NodeValue &node_value_)
        : position(position_), node_value(node_value_) {
      }
    };

    // true, underlying edge
    struct EdgeInternal {
      NodeID node1_id;
      NodeID node2_id;
      EdgeValue edge_value;
      EdgeInternal(const NodeID &node1_id_, const NodeID &node2_id_, const EdgeValue &edge_value_)
        : node1_id(node1_id_), node2_id(node2_id_), edge_value(edge_value_) {
      }
    };

    // underlying data structures storing graph information
    NodeVector nodes_;
    EdgeVector edges_;

};
