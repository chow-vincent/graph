/**
 * Implementation of mass-spring system using Graph class.
 * Can be used for simple physics simulations of objects represented as Graphs.
 */

#include <omp.h>
#include <thrust/system/omp/execution_policy.h>
#include <thrust/for_each.h>

#include <fstream>
#include <chrono>
#include <thread>
#include <cmath>

#include "utils/sfml_viewer.hpp"
#include "utils/util.hpp"
#include "utils/color.hpp"
#include "utils/point.hpp"

#include "graph.hpp"

// gravity in meters/sec^2
static constexpr double grav = 9.81;

// spring constant
static constexpr double kSpring = 100;

// rest length (default uniform rest length for Prob1).
static constexpr double kRestLength = 0.2;

// custom "value" to store in graph Nodes
struct NodeData {
  Point vel;   // Node "velocity"
  double mass; // Node "mass"
  NodeData() : vel(0), mass(1) {}
};

// custom "value" to store rest length of graph Edges
struct EdgeData {
  double length;
  EdgeData() : length(0) {}
};

// define aliases
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::Node;
using Edge = typename GraphType::Edge;

// update position functor for thrust::for_each in symp_euler_step
struct UpdatePosition {
  const double dt_;

  void operator()(Node n) {
    n.position() += n.value().vel * dt_;
  }

  UpdatePosition(const double dt) : dt_(dt) {
  }
};

// update velocity functor for thrust::for_each in symp_euler_step
template <typename F>
struct UpdateVelocity {
  double t_;
  double dt_;
  F force_;

  void operator()(Node n) {
    n.value().vel += force_(n, t_) * (dt_ / n.value().mass);
  }

  UpdateVelocity(double t, double dt, F force)
    : t_(t), dt_(dt), force_(force) {
  }
};

/** Modify a graph's nodes using symplectic Euler method with given Node
 *  forces and constraints. Operations that can be parallelized are parallized
 *  with thrust and OpenMP.
 *
 * @param[in,out] g      Graph
 * @param[in]     t      Current time (for time-dependent forces/constraints)
 * @param[in]     dt     Time step
 * @param[in]     force  Functor defining force on Node
 * @return        next time step (@a t + @a dt)
 *
 * @tparam G graph where Nodes store mass and velocities, and Edges store lengths
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C functor called to apply constraints on Graph Nodes.
 *
 */
template <typename G, typename F, typename C>
double symp_euler_step(G &g, double t, double dt, F force, C constraint) {

  // update positions
  thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), UpdatePosition(dt));

  // apply constraints
  constraint(t);

  // update velocities
  thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), UpdateVelocity<F>(t, dt, force));

  return t + dt;
}

// gravitational force functor
struct GravityForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return n.value().mass*Point(0.0, 0.0, -grav);
  }
};

// spring force functor
struct MassSpringForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;

    Point total_force = Point(0.0, 0.0, 0.0);
    Point pos_i = n.position();

    // add up spring forces
    for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      Point pos_j = (*iter).node2().position();

      const double length = (*iter).length();
      total_force += -kSpring*((pos_i - pos_j)/length)*(length - (*iter).value().length);
    }
    return total_force;
  }
};

// damping force functor
struct DampingForce {
  const double damping_coeff;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    return -damping_coeff*n.value().vel;
  }
  DampingForce(const double damping_coeff_) : damping_coeff(damping_coeff_) {
  }
};

// functor for combining three forces
template <typename F1, typename F2, typename F3>
struct MakeCombinedForce {
  F1 f1;
  F2 f2;
  F3 f3;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n, t) + f2(n, t) + f3(n, t);
  }
  MakeCombinedForce(F1 f1_, F2 f2_, F3 f3_)
    : f1(f1_), f2(f2_), f3(f3_) {
  }
};

// functor for combining two forces, (template specialization used)
template <typename F1, typename F2>
struct MakeCombinedForce<F1, F2, void> {
  F1 f1;
  F2 f2;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1(n, t) + f2(n, t);
  }
  MakeCombinedForce(F1 f1_, F2 f2_)
    : f1(f1_), f2(f2_){
  }
};

// for keeping track of pins
struct Pin {
  Point &pos;
  Point ref_pos;
  Point &vel;

  Pin(Point &pos_, Point ref_pos_, Point &vel_)
    : pos(pos_), ref_pos(ref_pos_), vel(vel_) {
    }
};

// pin constraint (keep pins at same position and zero vel)
struct PinConstraint {
  PinConstraint(std::vector<Pin> pins) {
    for (auto elt : pins) {
      elt.pos = elt.ref_pos;          // reset position
      elt.vel = Point(0.0, 0.0, 0.0); // reset velocity
    }
  }
};

// plane constraint functor
struct PlaneConstraint {
  const double kPlaneConstraint = -0.75;

  template <typename NODE>
  void operator()(NODE n, double t) {
    (void) t;
    if (n.position().z < kPlaneConstraint) {
      n.position().z = kPlaneConstraint; // set to nearest point on the plane
      n.value().vel.z = 0.0;             // set z comp to zero
    }
  }
};

// helper functions
Point project_onto_sphere(const Point &p, const Point &center, const double radius) {
  Point p_sphere = p - center;
  p_sphere = (p_sphere/norm_2(p_sphere))*(radius);
  return p_sphere + center;
}

Point compute_orthogonal_vel(const Point &vel, const Point &center) {
  Point R = (vel - center)/norm_2(vel - center);
  return inner_prod(vel, R)*R;
}

// functor for sphere constraint
struct SphereConstraint {
  const Point kCenter = Point(0.5, 0.5, -0.5);
  const double kRadius = 0.15;

  template <typename NODE>
  void operator()(NODE n, double t) {
    (void) t;
    if (norm_2(n.position() - kCenter) < kRadius) {
      n.position()  = project_onto_sphere(n.position(), kCenter, kRadius);
      n.value().vel -= compute_orthogonal_vel(n.value().vel, kCenter);
    }
  }
};

// functor for transforming node to point iterator
struct NodeToPoint {
  template <typename NODE>
  Point operator()(NODE n) {
    return n.position();
  }
};

// functor for transforming incident edge to node2 point iterator
struct EdgeToPoint {
  template <typename EDGE>
  Point operator()(EDGE e) {
    return e.node2().position();
  }
};

template <typename C1, typename C2>
struct ApplyIterativeConstraint {
  C1 c1;
  C2 c2;
  double t;

  template <typename NODE>
  void operator()(NODE n) {
    c1(n, t); c2(n, t);
  }

  ApplyIterativeConstraint(C1 c1_, C2 c2_, double t_) : c1(c1_), c2(c2_), t(t_) {
  }
};

template <typename C1, typename C2, typename C3>
struct MakeCombinedConstraint {
  C1 c1;
  C2 c2;
  C3 c3;
  const GraphType &graph;

  void operator()(double t) {

    // c1 (pin) constraint applied when its constructor called

    // apply c2, c3 constraints
    thrust::for_each(thrust::omp::par, graph.node_begin(), graph.node_end(),
                      ApplyIterativeConstraint<C2, C3>(c2, c3, t));
  }

  MakeCombinedConstraint(C1 c1_, C2 c2_, C3 c3_, const GraphType &graph_)
    : c1(c1_), c2(c2_), c3(c3_), graph(graph_) {
  }
};

int main(int argc, char** argv)
{
  // check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // construct graph
  GraphType graph;
  std::vector<GraphType::Node> nodes;

  // nodes_file from first input arg
  std::ifstream nodes_file(argv[1]);

  // each line of node_file is a 3D point, add to graph
  Point p;
  while (GraphUtil::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // tets_file from second input arg
  std::ifstream tets_file(argv[2]);

  // each line of tets_file as four ints referring to node IDs
  std::array<int,4> t;
  while (GraphUtil::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // set initial conditions for Nodes
  std::vector<Pin> pins;
  const std::size_t num_nodes = graph.size();
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    Node node = *it;
    node.value().mass = 1.0/num_nodes;
    if (node.position() == Point(0.0, 0.0, 0.0) || node.position() == Point(1.0, 0.0, 0.0)) {
      pins.emplace_back(node.position(), node.position(), node.value().vel);
    }
  }

  // set initial condition for edges (rest lengths)
  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it) {
    (*it).value().length = (*it).length();
  }

  // print num of nodes and edges
  std::cout << "# Nodes: " << graph.num_nodes() << std::endl;
  std::cout << "# Edges: " << graph.num_edges() << std::endl;

  // launch viewer
  GraphUtil::SFML_Viewer viewer;

  // load graph into viewer
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  // center view
  viewer.center_view();

  // setup thread to enable simultaneous interaction and simulation
  bool interrupt_sim_thread = false;
  auto sim_thread = std::thread([&](){

    // setup simulation
    double dt = 0.0005;
    double t_start = 0;
    double t_end = 5.0;

    // begin simulation
    for (double t = t_start; t < t_end && !interrupt_sim_thread; t += dt) {

      symp_euler_step(graph, t, dt,
        MakeCombinedForce<GravityForce, MassSpringForce, DampingForce>(
          GravityForce(), MassSpringForce(), DampingForce(1.0/graph.size())),
        MakeCombinedConstraint<PinConstraint, PlaneConstraint, SphereConstraint>(
          PinConstraint(pins), PlaneConstraint(), SphereConstraint(), graph));

      viewer.clear();
      node_map.clear();

      // update viewer with updated Graph Nodes and Edges
      viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
      viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
      viewer.set_label(t);
    }

  });  // simulation thread

  viewer.event_loop();

  // killing window means returned from event loop
  interrupt_sim_thread = true;
  sim_thread.join();

  return 0;
}
