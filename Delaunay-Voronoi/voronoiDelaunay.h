#include<iostream>
#include<vector>
#include<algorithm>

#define double long double

using namespace std;

const double eps = 1e-5;

struct Point {
	double x, y;
	Point() : x(0), y(0) {}
	Point(const double& x, const double& y) : x(x), y(y) {}

	bool operator==(const Point& other) const {
		return abs(other.x - x) <= eps && abs(other.y - y) <= eps;
	}

	bool operator!=(const Point& other) const {
		return !(operator==(other));
	}

	friend ostream& operator<<(ostream& os, const Point& p) {
		os << "(" << p.x << ",  " << p.y << ")";
		return os;
	}
};

struct Edge {
	Point p1, p2;

	Edge(const Point& p1, const Point& p2) : p1(p1), p2(p2) {}

	bool operator==(const Edge& other) const {
		return (p1 == other.p1 && p2 == other.p2 || p2 == other.p1 && p1 == other.p2);
	}

	bool operator!=(const Edge& other) const {
		return !(operator==(other));
	}

	friend std::ostream& operator<<(std::ostream& os, const Edge& e) {
		os << e.p1 << ", " << e.p2;
		return os;
	}
};

struct Cell {
	Point site;
	vector<Edge> Edges;

	Cell(Point site) : site(site) {};
	Cell(Point site, vector<Edge> Edges) : site(site), Edges(Edges) { };
};

// Ax + By = C
struct Line {
	double A, B, C;
};


//returns the perpendicular bisector of points P, Q
Line findPB(Point P, Point Q) {
	double xm = (P.x + Q.x) / 2, ym = (P.y + Q.y) / 2;

	double A = Q.y - P.y;
	double B = P.x - Q.x;
	double C = A * P.x + B * P.y;

	//pb equation: A1 * x + B1* y = C1
	double A1 = -B;
	double B1 = A;
	double C1 = A * ym - B * xm;

	return { A1,B1,C1 };
}


//returns intersection point of the edge e with the line (e not parallel with line)
Point getIntersection(Edge e, Line line) {
	Point P = e.p1, Q = e.p2;

	// A1 * x + B1 * y = C - line containing edge e
	double A1 = Q.y - P.y;
	double B1 = P.x - Q.x;
	double C1 = A1 * Q.x + B1 * Q.y; 
	double A2 = line.A, B2 = line.B, C2 = line.C;
	double det = A1 * B2 - A2 * B1;

	// compute intersection point of the 2 lines
	double x = (B2 * C1 - B1 * C2) / det;
	double y = (A1 * C2 - A2 * C1) / det;
	return { x,y };
}

// returns 1 if point p lies on the left side of line and 0 otherwise
int getSide(Line line, Point p) {
	return line.A * p.x + line.B * p.y > line.C;
}	

struct VoronoiDelaunay {
	vector<Edge> Edges;
	vector<Edge> Delaunay;
};

VoronoiDelaunay getVDEdgesFromCells(const vector<Point> &Sites, const vector<Cell> &Cells);


//Input: vector of points, width and height of the window
//Output: Voronoi and Delaunay triangulation of the points
VoronoiDelaunay computeTriangulation(vector<Point> points, double width, double height) {

	// add randomness
	for (auto &point : points) {
		float r1 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 0.01)) + 0.0005;
		float r2 = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 0.01)) + 0.0005;

		point.x += r1;
		point.y += r2;
	}
	vector<Cell> Cells;

	//creating boundary cells
	double k = 1000, k2 = 30; 
	vector<Edge> boundary1{ Edge(Point(0,0), Point(0, -k * height)),
					    Edge(Point(0, -k * height), Point(-k * width, -k * height)),
					    Edge(Point(-k * width, -k * height), Point(-k * width, 0)),
					    Edge(Point(-k * width, 0), Point(0, 0)) };

	vector<Edge> boundary2{ Edge(Point(0, 0), Point(-k * width, 0)),
					    Edge(Point(-k * width, 0), Point(-k * width, k * height)),
					    Edge(Point(-k * width, k * height), Point(0, k * height)),
					    Edge(Point(0, k * height), Point(0, 0)) };

	vector<Edge> boundary3{ Edge(Point(0, 0), Point(0, k * height)),
					    Edge(Point(0, k * height), Point(k * height, k * height)),
					    Edge(Point(k * height, k * height), Point(k * height, 0)),
					    Edge(Point(k * height, 0), Point(0, 0)) };

	vector<Edge> boundary4{ Edge(Point(0, 0), Point(k * width, 0)),
					    Edge(Point(k * width, 0), Point(k * width, -k * height)),
					    Edge(Point(k * width, -k * height), Point(0, -k * height)),
					    Edge(Point(0, -k * height), Point(0, 0)) };


	Point p1(-k2 *width, -k2 * height);
	Point p2(-k2 * width, k2 * height);
	Point p4(k2 * width, -k2 * height);
	Point p3(k2 * width, k2 * height);

	Cells.push_back(Cell(p1, boundary1));
	Cells.push_back(Cell(p2, boundary2));
	Cells.push_back(Cell(p3, boundary3));
	Cells.push_back(Cell(p4, boundary4));

	for (const auto& site : points) {
		Cell currCell(site);		// initialize cell for current point
		for (auto& cell : Cells) {    // iterate over all cells
			Line pb = findPB(site, cell.site); 
			vector<Point> critPoints;
			int sideSite = getSide(pb, site);
			int inters = 0;
			
			for (auto edge = cell.Edges.begin(); edge != cell.Edges.end(); ) {
				// for each edge of the cell test the spatial relation with pb
				int sideP1 = getSide(pb, edge->p1);
				int sideP2 = getSide(pb, edge->p2);

				if (sideSite == sideP1 && sideP1 == sideP2) { // edge is on the side of current site
					edge = cell.Edges.erase(edge);	      // delete edge
					continue;
				} else if (sideP1 != sideP2) {			 // edge intersects pb
					Point intersPoint = getIntersection(*edge, pb);
					critPoints.push_back(intersPoint);		 // store the intersection point 

					if (sideSite == sideP1) edge->p1 = intersPoint;	//clip the edge
					else if (sideSite == sideP2) edge->p2 = intersPoint;
				}
				++edge;
			}

			
			if (critPoints.size()==2) {
				// create a new edge that connects the points of intersections and to the cells
				Edge newEdge(critPoints[0], critPoints[1]);
				cell.Edges.push_back(newEdge);
				currCell.Edges.push_back(newEdge);
			}
		}

		Cells.push_back(currCell);
	}
	
	//delete boundary cells
	Cells.erase(Cells.begin(), Cells.begin() + 4);
	return getVDEdgesFromCells(points, Cells);
}

// get the corresponding Voronoi diagram and delaunay edges from a vector of voronoi Cells
VoronoiDelaunay getVDEdgesFromCells(const vector<Point> &points, const vector<Cell> &Cells) {
	vector<pair<Edge, int> > allEdges;
	VoronoiDelaunay VD;

	//store all edges from Vonoi cells in the vector allEdges
	for (auto i = 0; i < Cells.size(); i++) {
		for (auto edge : Cells[i].Edges) allEdges.push_back({ edge, i });
	}

	vector<bool> isDuplicate(allEdges.size(), 0);

	for (auto it = allEdges.begin(); it != allEdges.end(); ++it) {
		for (auto it2 = it + 1; it2 != allEdges.end(); ++it2) { 
			int idx1 = it->second, idx2 = it2->second;

			//for each pair of edges from different cells check if they are equal
			if (idx1 != idx2 && it->first == it2->first) {
				VD.Delaunay.push_back(Edge(points[idx1], points[idx2]));  // add edge to triangulation edges
				isDuplicate[it2 - allEdges.begin()] = 1;			   // identify duplicate edges
			}
		}
	}

	for (int i = 0; i < allEdges.size(); ++i) {
		if (!isDuplicate[i]) VD.Edges.push_back(allEdges[i].first);
	}

	return VD;
}