#pragma once

#include <vector>
#include <cassert>

template <typename V, typename E>
class Graph {
  private:

    // forward declarations of underlying data
    struct NodeInternal;
    struct EdgeInternal;

  public:
    using NodeValue = V;
    using EdgeValue = E;

    using NodeID = std::size_t;
    using EdgeID = std::size_t;

    // aliases for underlying data structures
    using NodeVector = std::vector<NodeInternal>;
    using EdgeVector = std::vector<EdgeInternal>;
    using AdjEdgesVector = std::vector<EdgeID>;

    // forward declarations of iterators
    class IncidentIterator;

    // constructor for empty graph
    Graph() {
    }

    // default destructor
    ~Graph() = default;

    // Node class. Proxy for NodeInternal.
    // TODO:
    // - without totally_ordered, will need to define != operator myself
    class Node {
      public:
        Node() {
          node_id_ = -1;
          graph_ptr_ = nullptr;
        }

        Node(NodeID node_id, const Graph *graph_ptr)
          : node_id_(node_id), graph_ptr_(const_cast<Graph*>(graph_ptr)) {
        }

        // modifiable, non-const ref to position
        // TODO: precondition is assuming node is valid for now
        // may need asserts later though for debugging
        Point &position() {
          return (graph_ptr_->nodes_)[node_id_].position;
        }

        const Point &position() const {
          return (graph_ptr_->nodes_)[node_id_].position;
        }

        // ID, synonymous with index into NodeVector
        NodeID get_node_id() const {
          return node_id_;
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

        // TODO: come back to this after incident iterator implemented
        IncidentIterator edge_begin() const {
        }

        IncidentIterator edge_end() const {
        }

        // comparators
        bool operator==(const Node &n) const {
          return (graph_ptr_ == n.get_graph_pointer()) &&
                  (node_id_ == n.get_node_id());
        }

        // TODO: see if this works
        bool operator!=(const Node &n) const {
          return !(*this == n);
        }

        bool operator<(const Node &n) const {
          return node_id_ < n.get_node_id();
        }

        Graph *get_graph_pointer() const {
          return graph_ptr_;
        }

      private:
        friend class Graph;

        // ID of node proxy.
        // Equivalent to index into NodeVector for this implementation.
        NodeID node_id_;

        // pointer to the graph to which this node belongs
        Graph *graph_ptr_;
    };

    std::size_t size() const {
      return nodes_.size();
    }

    std::size_t num_nodes() const {
      return size();
    }

    Node add_node(const Point &position,
      const NodeValue &node_value = NodeValue()) {
      nodes_.emplace_back(position, node_value);
      return Node(nodes_.size(), this);
    }

    // TODO:
    // - change this to return bool?
    std::size_t remove_node(const Node &n) {
    }

    // TODO?:
    // NodeIterator remove_node(NodeIterator n_it) {
    // }

    // TODO:
    bool has_node(const Node &n) const {
    }

    // TODO:
    // - maybe change size_t to NodeID and i to id
    Node node(std::size_t i) const {
    }

    // TODO:
    // - without totally_ordered, will need to define != operator myself
    class Edge {
    };

    // TODO:
    std::size_t num_edges() const {
    }

    // TODO:
    // - maybe change size_t to EdgeID and i to id
    Edge edge(std::size_t i) const {
    }

    // TODO:
    std::size_t remove_edge(const Edge &e) {
    }

    // TODO:
    std::size_t remove_edge(const Node &a, const Node &b) {
    }

    // TODO?:
    // EdgeIterator remove_edge(EdgeIterator e_it) {
    // }

    // TODO:
    bool has_edge(const Node &a, const Node &b) const {
    }

    // TODO:
    Edge add_edge(const Node& a, const Node& b,
      const edge_value_type &edge_value = edge_value_type()) {
    }

    // TODO:
    void clear() {
    }

    // TODO:
    class NodeIterator {
    };

    // TODO:
    NodeIterator node_begin() const {
    }

    // TODO:
    NodeIterator node_end() const {
    }

    // TODO:
    class IncidentIterator {
    };

    // TODO:
    class EdgeIterator {
    };

    // TODO:
    EdgeIterator edge_begin() const {
    }

    EdgeIterator edge_end() const {
    }

  private:

    // TODO:
    struct NodeInternal {
      // TODO: do I even need to store the ID, now that I enforce ID == index?
      // NodeID node_id;
      Point position;
      NodeValue node_value;
      AdjEdgesVector adjacent_edges;
      NodeInternal(const NodeID node_id_, const Point &position_,
        const NodeValue &node_value_)
        : node_id(node_id_), position(position_), node_value(node_value_) {
      }
    };

    // TODO:
    struct EdgeInternal {
      // EdgeID edge_id;
      NodeID node1_id;
      NodeID node2_id;
      EdgeValue edge_value;
      EdgeInternal(const EdgeID edge_id_, const NodeID &node1_id_,
        const NodeID &node2_id_, const EdgeVAlue &edge_value_)
        : edge_id(edge_id_), node1_id(node1_id_),
          node2_id(node2_id_), edge_value(edge_value_) {
      }
    };

    NodeVector nodes_;
    EdgeVector edges_;

}
