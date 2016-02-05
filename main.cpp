#include "mfem.hpp"

#include <algorithm>
#include <cfloat>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

//#define DEBUG

using namespace std;
using namespace mfem;

const int DIM = 2; // working dimension

const double SAME_POINT_TOLERANCE = 1e-8;

const int SURFACE_ID              = 4;
const int COARSE_EDGE_ID          = 5;
const int NEW_LINE_ID             = 7;
const int BOUNDARY_COARSE_EDGE_ID = 9;

/**
 * Distance between two vertices
 */
double distance_between(const Vertex &a, const Vertex &b)
{
  double dist = 0.0;
  for (int d = 0; d < DIM; ++d)
    dist += (a(d)-b(d)) * (a(d)-b(d));
  return sqrt(dist);
}

/**
 * Special structure for sorting intermediate points inside edges - the points
 * of intersections with other edges
 */
struct CompareInterPoints
{
  int point_beg; ///< index of beginning point of an edge
  
  ///< all points of the mesh (vertices and intersections)
  const vector<Vertex>& points; 

  CompareInterPoints(int point_beg_, const vector<Vertex> &points_)
    : point_beg(point_beg_)
    , points(points_)
  { }

  /**
   * Actual comparison between two intermediate points within an edge structure
   */
  bool operator () (int a, int b)
  {
    // which point is closer to the beginning of the edge comes first when 
    // sorting
    const double dist_a = distance_between(points[point_beg], points[a]);
    const double dist_b = distance_between(points[point_beg], points[b]);
    return dist_a < dist_b;
  }
};

/**
 * Check if the given vector contains a given number.
 */
bool vector_has_number(const vector<int>& thevector, int number)
{
  return (find(thevector.begin(), thevector.end(), number) != thevector.end());
}

/**
 * Any edge in the geometry - it describes the edges of an initial coarse grid
 * and any new lines that are inserted
 */
struct Edge
{
  Edge(int beg = -1, int end = -1, int id = 1)
    : point_beg(beg)
    , point_end(end)
    , intermediate_points()
    , ID(id)
  { }
  
  Edge(const Edge &e)
    : point_beg(e.point_beg)
    , point_end(e.point_end)
    , intermediate_points(e.intermediate_points)
    , ID(e.ID)
  { }
  
  /**
   * Add an point index to a list of intermediate points - if it's not there yet
   */
  void add_index(int ind)
  {
    if (ind == point_beg || ind == point_end) return;
    if (find(intermediate_points.begin(), intermediate_points.end(), ind) ==
        intermediate_points.end())
      intermediate_points.push_back(ind);
  }

  /**
   * Sort the intermediate points in such a way that they go from the beginning
   * of the edge to the end point
   */
  void sort_intermediate_points(const vector<Vertex> &points)
  {
    if (intermediate_points.empty() ||
        intermediate_points.size() == 1)
      return; // there is nothing to sort

    // special sorting object
    CompareInterPoints compare_inter_points(point_beg, points);

    sort(intermediate_points.begin(),
         intermediate_points.end(),
         compare_inter_points);
  }

  /**
   * Check if the edge lies on a boundary
   */
  bool on_boundary(const vector<int>& left_boundary_points,
                   const vector<int>& right_boundary_points,
                   const vector<int>& bottom_boundary_points,
                   const vector<int>& top_boundary_points) const
  {
    if (vector_has_number(left_boundary_points, point_beg) &&
        vector_has_number(left_boundary_points, point_end))
      return true;
    if (vector_has_number(right_boundary_points, point_beg) &&
        vector_has_number(right_boundary_points, point_end))
      return true;
    if (vector_has_number(bottom_boundary_points, point_beg) &&
        vector_has_number(bottom_boundary_points, point_end))
      return true;
    if (vector_has_number(top_boundary_points, point_beg) &&
        vector_has_number(top_boundary_points, point_end))
      return true;
    return false;
  }
  
  int point_beg, point_end; ///< beginning and end points (indices)
  
  /// < indices of intermediate points (intersections)
  vector<int> intermediate_points;

  int ID; ///< Identification number for Physical Line in Gmsh .geo file
};

/**
 * Get an index of a point in the list of the existing points. If there is no
 * such point in the list, we add it to the list, and return the new index.
 */
int get_point_index(vector<Vertex> &points, const Vertex &point_exam)
{
  for (size_t pi = 0; pi < points.size(); ++pi)
  {
    if (distance_between(points[pi], point_exam) < SAME_POINT_TOLERANCE)
      return pi; // index of the existing point
  }

  points.push_back(point_exam);
  return points.size()-1; // index of the new point
}

/**
 * Almost the same as 'get_point_index', but if the point can't be found in the
 * list, an exception is thrown, because here we assume that the point must be
 * there.
 */
int find_point_index(const vector<Vertex> &points, const Vertex &point_exam)
{
  for (size_t pi = 0; pi < points.size(); ++pi)
  {
    if (distance_between(points[pi], point_exam) < SAME_POINT_TOLERANCE)
      return pi;
  }

  throw runtime_error("find_point_index: not found");
  return -1;
}

/**
 * Check if the number 'x' is between 'a' and 'b' with respect to some tolerance
 */
bool in_range(double a, double b, double x)
{
  const double tol = SAME_POINT_TOLERANCE;
  double min_ = min(a, b);
  double max_ = max(a, b);
  if (x+tol >= min_ && x-tol <= max_)
    return true;
  return false;
}

/**
 * Compute the coordinates of the intersection points between two segments
 * defined by vertices: segment1 = beg1 -> end1, segment2 = beg2 -> end2.
 * The function returns a number: 0 if the segments intersect (in this
 * case the point of intersection is (x_intersect, y_intersect)), 1
 * if the segments do not intersect (in this case x_intersect and y_intersect
 * may contains garbage), and 2 if they coincide (in this case the new edge
 * shouldn't be added).
 */
int compute_intersection(const Vertex& beg1, const Vertex& end1,
                         const Vertex& beg2, const Vertex& end2,
                         double &x_intersect, double &y_intersect)
{
  const double tol = SAME_POINT_TOLERANCE;

  if (fabs(end2(0) - beg2(0)) < tol)
    throw runtime_error("The new line is vertical - not expected");

  const double a2 = (end2(1) - beg2(1)) / (end2(0) - beg2(0));
  const double b2 = beg2(1) - a2 * beg2(0);

  if (fabs(end1(0) - beg1(0)) < tol) // vertical line
  {
    x_intersect = beg1(0); // only possible intersection
    if (!in_range(beg2(0), end2(0), x_intersect))
      return 1; // no intersection
  }
  else
  {
    const double a1 = (end1(1) - beg1(1)) / (end1(0) - beg1(0));
    const double b1 = beg1(1) - a1 * beg1(0);

#if defined(DEBUG)
    cout << "a1 = " << a1 << " b1 = " << b1 << " a2 = " << a2 << " b2 = " << b2 << " ";
#endif

    if (fabs(a1 - a2) < tol)
    {
      if (fabs(b1 - b2) < tol)
        return 2; // the segments coincide
      else
        return 1; // no intersection - segments are parallel
    }

#if defined(DEBUG)
    cout << "beg1(0) = " << beg1(0) << " beg1(1) = " << beg1(1) << " ";
    cout << "end1(0) = " << end1(0) << " end1(1) = " << end1(1) << " ";
    cout << "beg2(0) = " << beg2(0) << " beg2(1) = " << beg2(1) << " ";
    cout << "end2(0) = " << end2(0) << " end2(1) = " << end2(1) << " ";
#endif
  
    x_intersect = (b2 - b1) / (a1 - a2);
    if (!in_range(beg1(0), end1(0), x_intersect) ||
        !in_range(beg2(0), end2(0), x_intersect))
      return 1; // no intersection
  }
    
  y_intersect = a2 * x_intersect + b2;

#if defined(DEBUG)
  cout << "x_inter = " << x_intersect << " y_inter = " << y_intersect << endl;
#endif

  if (!in_range(beg1(1), end1(1), y_intersect) ||
      !in_range(beg2(1), end2(1), y_intersect))
    return 1; // no intersection
    
  return 0; // there is an intersection at the point (x_intersect, y_intersect)
}

/**
 * Find an intersection points between all the existing edges and a new one.
 * Then the new edge is appended to the list of existing edges.
 */
void find_intersection_points(vector<Edge> &edges, 
                              vector<Vertex> &points, 
                              Edge &new_edge)
{
  const Vertex& beg2 = points[new_edge.point_beg];
  const Vertex& end2 = points[new_edge.point_end];

  for (size_t ei = 0; ei < edges.size(); ++ei)
  {
    const Vertex& beg1 = points[edges[ei].point_beg];
    const Vertex& end1 = points[edges[ei].point_end];
    
    double x_intersect, y_intersect;
    int intersect = compute_intersection(beg1, end1, beg2, end2,
                                         x_intersect, y_intersect);
    if (intersect == 0)
    {
      const int ind = get_point_index(points, Vertex(x_intersect, y_intersect));
      edges[ei].add_index(ind);
      new_edge.add_index(ind);
    }
  }
  
  edges.push_back(new_edge);
}

/**
 * Get a list of indices of points that have the same x-coordinate as a given
 * one.
 */
void get_list_of_point_indices_x(const vector<Vertex>& points,
                                 double x_coord,
                                 vector<int>& points_indices)
{
  points_indices.clear();
  for (size_t pi = 0; pi < points.size(); ++pi)
    if (fabs(points[pi](0) - x_coord) < SAME_POINT_TOLERANCE)
      points_indices.push_back(pi);
}

/**
 * Get a list of indices of points that have the same y-coordinate as a given
 * one.
 */
void get_list_of_point_indices_y(const vector<Vertex>& points,
                                 double y_coord,
                                 vector<int>& points_indices)
{
  points_indices.clear();
  for (size_t pi = 0; pi < points.size(); ++pi)
    if (fabs(points[pi](1) - y_coord) < SAME_POINT_TOLERANCE)
      points_indices.push_back(pi);
}

/**
 * Create a list of well-connected lines lying on a left boundary, and insert
 * their number in the array which is going to be used in Gmsh .geo file.
 */
void make_left_boundary_line(ostream& out,
                             const vector<Vertex>& points,
                             double x_min,
                             int line_beg,
                             int line_end,
                             vector<int>& points_indices)
{
  get_list_of_point_indices_x(points, x_min, points_indices);
  CompareInterPoints compare_inter_points(line_beg, points);
  sort(points_indices.begin(), points_indices.end(), compare_inter_points);
  if (points_indices.front() != line_beg || points_indices.back() != line_end)
    throw runtime_error("Incorrect list of left boundary nodes");
  for (size_t i = 0; i < points_indices.size()-1; ++i)
  {
    out << "l=newl; Line(l)={" << points_indices[i]+1
        << "," << points_indices[i+1]+1 << "}; "
        << "left_boundary[" << i << "]=l;\n";
  }
}

/**
 * Create a list of well-connected lines lying on a right boundary, and insert
 * their number in the array which is going to be used in Gmsh .geo file.
 */
void make_right_boundary_line(ostream& out,
                              const vector<Vertex>& points,
                              double x_max,
                              int line_beg,
                              int line_end,
                              vector<int>& points_indices)
{
  get_list_of_point_indices_x(points, x_max, points_indices);
  CompareInterPoints compare_inter_points(line_beg, points);
  sort(points_indices.begin(), points_indices.end(), compare_inter_points);
  if (points_indices.front() != line_beg || points_indices.back() != line_end)
    throw runtime_error("Incorrect list of left boundary nodes");
  for (size_t i = 0; i < points_indices.size()-1; ++i)
  {
    out << "l=newl; Line(l)={" << points_indices[i]+1
        << "," << points_indices[i+1]+1 << "}; "
        << "right_boundary[" << i << "]=l;\n";
  }
}

/**
 * Create a list of well-connected lines lying on a bottom boundary, and insert
 * their number in the array which is going to be used in Gmsh .geo file.
 */
void make_bottom_boundary_line(ostream& out,
                               const vector<Vertex>& points,
                               double y_min,
                               int line_beg,
                               int line_end,
                               vector<int>& points_indices)
{
  get_list_of_point_indices_y(points, y_min, points_indices);
  CompareInterPoints compare_inter_points(line_beg, points);
  sort(points_indices.begin(), points_indices.end(), compare_inter_points);
  if (points_indices.front() != line_beg || points_indices.back() != line_end)
    throw runtime_error("Incorrect list of bottom boundary nodes");
  for (size_t i = 0; i < points_indices.size()-1; ++i)
  {
    out << "l=newl; Line(l)={" << points_indices[i]+1
        << "," << points_indices[i+1]+1 << "}; "
        << "bottom_boundary[" << i << "]=l;\n";
  }
}

/**
 * Create a list of well-connected lines lying on a top boundary, and insert
 * their number in the array which is going to be used in Gmsh .geo file.
 */
void make_top_boundary_line(ostream& out,
                            const vector<Vertex>& points,
                            double y_max,
                            int line_beg,
                            int line_end,
                            vector<int>& points_indices)
{
  get_list_of_point_indices_y(points, y_max, points_indices);
  CompareInterPoints compare_inter_points(line_beg, points);
  sort(points_indices.begin(), points_indices.end(), compare_inter_points);
  if (points_indices.front() != line_beg || points_indices.back() != line_end)
    throw runtime_error("Incorrect list of top boundary nodes");
  for (size_t i = 0; i < points_indices.size()-1; ++i)
  {
    out << "l=newl; Line(l)={" << points_indices[i]+1
        << "," << points_indices[i+1]+1 << "}; "
        << "top_boundary[" << i << "]=l;\n";
  }
}

/**
 * Write down the geometry of intersected lines to build a mesh later on.
 */
void make_geo_file(ostream &out,
                   const Vertex &min_coord,
                   const Vertex &max_coord,
                   const vector<Edge> &edges,
                   const vector<Vertex> &points)
{
  out << "cl = 10;\n"; // that can be changed manually later on
  
  for (size_t pi = 0; pi < points.size(); ++pi)
    out << "Point(" << pi+1 << ")={"
        << points[pi](0) << "," << points[pi](1) << ",0,cl};\n";

  int bdr_i1 = find_point_index(points, Vertex(min_coord(0), min_coord(1)));
  int bdr_i2 = find_point_index(points, Vertex(max_coord(0), min_coord(1)));
  int bdr_i3 = find_point_index(points, Vertex(max_coord(0), max_coord(1)));
  int bdr_i4 = find_point_index(points, Vertex(min_coord(0), max_coord(1)));

  vector<int> left_points, right_points, bottom_points, top_points;

  make_bottom_boundary_line(out, points, min_coord(1), bdr_i1, bdr_i2, bottom_points);
  make_right_boundary_line (out, points, max_coord(0), bdr_i2, bdr_i3, right_points);
  make_top_boundary_line   (out, points, max_coord(1), bdr_i3, bdr_i4, top_points);
  make_left_boundary_line  (out, points, min_coord(0), bdr_i4, bdr_i1, left_points);


  out << "Line Loop(1)={bottom_boundary[],right_boundary[],top_boundary[],left_boundary[]};\n";
  out << "Plane Surface(1)={1};\n";
  out << "Physical Surface(" << SURFACE_ID << ")={1};\n";
  out << "Physical Line(" << BOUNDARY_COARSE_EDGE_ID
      << ")={bottom_boundary[],right_boundary[],top_boundary[],left_boundary[]};\n";

  int n_coarse_edge = 0, n_new_line = 0;
  for (size_t ei = 0; ei < edges.size(); ++ei)
  {
    const Edge &edge = edges[ei];
    if (edge.on_boundary(left_points, right_points, bottom_points, top_points))
      continue;
    out << "l=newl; Line(l)={" << edge.point_beg+1 << ",";
    for (size_t ip = 0; ip < edge.intermediate_points.size(); ++ip)
    {
      out << edge.intermediate_points[ip]+1 << "}; Line{l} In Surface{1}; ";
      switch (edge.ID)
      {
        case COARSE_EDGE_ID:
        {
          out << "coarse_edges[" << n_coarse_edge++ << "]=l;\n";
          break;
        }
        case NEW_LINE_ID:
        {
          out << "new_lines[" << n_new_line++ << "]=l;\n";
          break;
        }
        default:
          throw runtime_error("Unexpected edge ID");
      }
      out << "l=newl; Line(l)={" << edge.intermediate_points[ip]+1 << ",";
    }
    out << edge.point_end+1 << "}; Line{l} In Surface{1}; ";
    switch (edge.ID)
    {
      case COARSE_EDGE_ID:
      {
        out << "coarse_edges[" << n_coarse_edge++ << "]=l;\n";
        break;
      }
      case NEW_LINE_ID:
      {
        out << "new_lines[" << n_new_line++ << "]=l;\n";
        break;
      }
      default:
        throw runtime_error("Unexpected edge ID");
    }
  }
  out << "Physical Line(" << COARSE_EDGE_ID << ")={coarse_edges[]};\n";
  out << "Physical Line(" << NEW_LINE_ID    << ")={new_lines[]};\n";
}

/**
 * Show the points
 */
void print_points(const vector<Vertex> &points)
{
  for (size_t pi = 0; pi < points.size(); ++pi)
    cout << pi+1 << " " << points[pi](0) << " " << points[pi](1) << "\n";
  cout << "======================\n";
}

/**
 * Show the edges
 */
void print_edges(const vector<Edge> &edges)
{
  for (size_t ei = 0; ei < edges.size(); ++ei)
  {
    const Edge &edge = edges[ei];
    cout << ei+1 << " : " << edge.point_beg+1 << " -> " << edge.point_end+1 << " : ";
    for (size_t ip = 0; ip < edge.intermediate_points.size(); ++ip)
      cout << edge.intermediate_points[ip]+1 << " ";
    cout << "\n";
  }
  cout << "======================\n";
}



//==============================================================================
//
// Main function
//
//==============================================================================
int main(int argc, char *argv[])
{
   const char *coarsemesh_file = "undefined-file.msh";
   const char *lines_file      = "undefined-file.txt";
   const char *result_file     = "result.geo";

   OptionsParser args(argc, argv);
   args.AddOption(&coarsemesh_file, "-m", "--mesh-file", "Coarse mesh file");
   args.AddOption(&lines_file, "-l", "--lines-file",
                  "File containing descriptions of lines to be inserted");
   args.AddOption(&result_file, "-r", "--result-file",
                  "File containing resulting geometry (optional: result.geo by default)");
   args.Parse();
   if (argc == 1 || !args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   ifstream inmesh(coarsemesh_file);
   if (!inmesh)
   {
      cerr << "\nCan not open coarse mesh file: " << coarsemesh_file << endl;
      return 2;
   }
   const bool generate_edges  = 1;
   const bool refine          = 0;
   const bool fix_orientation = 0;
   Mesh mesh(inmesh, generate_edges, refine, fix_orientation);
   inmesh.close();

   if (mesh.Dimension() != DIM)
   {
     cerr << "\nCoarse mesh must be " << DIM << "D\n" << endl;
     return 2;
   }

   // copy all the mesh vertices in the list of points (it will be extended by
   // the intersection points)
   vector<Vertex> points(mesh.GetNV());

   Vertex min_coord(DBL_MAX, DBL_MAX);
   Vertex max_coord(DBL_MIN, DBL_MIN);
   for (int vi = 0; vi < mesh.GetNV(); ++vi)
   {
     double *v_coord = mesh.GetVertex(vi);
     for (int d = 0; d < DIM; ++d)
     {
       min_coord(d) = min(min_coord(d), v_coord[d]);
       max_coord(d) = max(max_coord(d), v_coord[d]);
       points[vi](d)= v_coord[d];
     }
   }
    
   cout << "min_coord = ( " << min_coord(0) << ", " << min_coord(1) << " )\n";
   cout << "max_coord = ( " << max_coord(0) << ", " << max_coord(1) << " )\n";

   // initial mesh edges - will be extended by the lines from the lines_file
   cout << "mesh edges: " << mesh.GetNEdges() << "\n";
   vector<Edge> edges(mesh.GetNEdges()); 
   for (int ei = 0; ei < mesh.GetNEdges(); ++ei)
   {
     Array<int> edge_vertices;
     mesh.GetEdgeVertices(ei, edge_vertices);
     edges[ei] = Edge(edge_vertices[0], edge_vertices[1], COARSE_EDGE_ID);
   }

   ifstream inlines(lines_file);
   if (!inlines)
   {
     cerr << "\nCan not open file with lines: " << lines_file << "\n" << endl;
     return 2;
   }

#if defined(DEBUG)
   print_edges(edges);
#endif

   string line;
   int n_inserted_lines = 0;
   int line_index = 0;
   while (getline(inlines, line))
   {
     ++line_index;
     if (line[0] == '#') continue;

     double x0, y0, x1, y1;
     istringstream istr(line);
     istr >> x0 >> y0 >> x1 >> y1;

     if (x0+SAME_POINT_TOLERANCE < min_coord(0) ||
         y0+SAME_POINT_TOLERANCE < min_coord(1) ||
         x1-SAME_POINT_TOLERANCE > max_coord(0) ||
         y1-SAME_POINT_TOLERANCE > max_coord(1))
     {
       cerr << "The edge from the file with lines is out of the mesh domain\n";
       cerr << "This edge is: (" << x0 << "," << y0 << ") -> (" 
            << x1 << "," << y1 << ") at line " << line_index << "\n";
       return 2;
     }
     
     const int i1 = get_point_index(points, Vertex(x0, y0));
     const int i2 = get_point_index(points, Vertex(x1, y1));

     Edge new_edge(i1, i2, NEW_LINE_ID);
     find_intersection_points(edges, points, new_edge);

#if defined(DEBUG)
     print_edges(edges);
#endif

     ++n_inserted_lines;
   }

   inlines.close();

   cout << "inserted lines: " << n_inserted_lines << "\n";

#if defined(DEBUG)
   print_edges(edges);
#endif

   for (size_t ei = 0; ei < edges.size(); ++ei)
     edges[ei].sort_intermediate_points(points);

#if defined(DEBUG)
   print_edges(edges);
#endif
 
   ofstream out_geo(result_file);
   if (!out_geo)
   {
     cerr << "\nCan not open file for resuling geometry: " << result_file << "\n" << endl;
     return 2;
   }
   make_geo_file(out_geo, min_coord, max_coord, edges, points);
   out_geo.close();

   cout << "file " << result_file << " created successfully\n\n";

   return 0;
}

