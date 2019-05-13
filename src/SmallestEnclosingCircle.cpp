
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include "SmallestEnclosingCircle.h"

using std::size_t;
using std::vector;


/*---- Members of struct Point ----*/

Point Point::subtract(const Point &p) const {
	return Point{x - p.x, y - p.y};
}


double Point::distance(const Point &p) const {
	return std::hypot(x - p.x, y - p.y);
}


double Point::cross(const Point &p) const {
	return x * p.y - y * p.x;
}


/*---- Members of struct Circle ----*/

const Circle Circle::INVALID{Point{0, 0}, -1};

const double Circle::MULTIPLICATIVE_EPSILON = 1 + 1e-14;


bool Circle::contains(const Point &p) const {
	return c.distance(p) <= r * MULTIPLICATIVE_EPSILON;
}


bool Circle::contains(const vector<Point> &ps) const {
	for (const Point &p : ps) {
		if (!contains(p))
			return false;
	}
	return true;
}


/*---- Smallest enclosing circle algorithm ----*/

static Circle makeSmallestEnclosingCircleOnePoint(const vector<Point> &points, size_t end, const Point &p);
static Circle makeSmallestEnclosingCircleTwoPoints(const vector<Point> &points, size_t end, const Point &p, const Point &q);

static std::default_random_engine randGen((std::random_device())());


// Initially: No boundary points known
Circle makeSmallestEnclosingCircle(const vector<Point> &points) {
	
	// Progressively add points to circle or recompute circle
	Circle c(Circle::INVALID);
	for (size_t i = 0; i < points.size(); i++) {
		const Point &p = points.at(i);
		if (c.r < 0 || !c.contains(p))
			c = makeSmallestEnclosingCircleOnePoint(points, i + 1, p);
	}
	return c;
}


// One boundary point known
static Circle makeSmallestEnclosingCircleOnePoint(const vector<Point> &points, size_t end, const Point &p) {
	Circle c{p, 0};
	for (size_t i = 0; i < end; i++) {
		const Point &q = points.at(i);
		if (!c.contains(q)) {
			if (c.r == 0)
				c = makeDiameter(p, q);
			else
				c = makeSmallestEnclosingCircleTwoPoints(points, i + 1, p, q);
		}
	}
	return c;
}


// Two boundary points known
static Circle makeSmallestEnclosingCircleTwoPoints(const vector<Point> &points, size_t end, const Point &p, const Point &q) {
	Circle circ = makeDiameter(p, q);
	Circle left = Circle::INVALID;
	Circle right = Circle::INVALID;
	
	// For each point not in the two-point circle
	Point pq = q.subtract(p);
	for (size_t i = 0; i < end; i++) {
		const Point &r = points.at(i);
		if (circ.contains(r))
			continue;
		
		// Form a circumcircle and classify it on left or right side
		double cross = pq.cross(r.subtract(p));
		Circle c = makeCircumcircle(p, q, r);
		if (c.r < 0)
			continue;
		else if (cross > 0 && (left.r < 0 || pq.cross(c.c.subtract(p)) > pq.cross(left.c.subtract(p))))
			left = c;
		else if (cross < 0 && (right.r < 0 || pq.cross(c.c.subtract(p)) < pq.cross(right.c.subtract(p))))
			right = c;
	}
	
	// Select which circle to return
	if (left.r < 0 && right.r < 0)
		return circ;
	else if (left.r < 0)
		return right;
	else if (right.r < 0)
		return left;
	else
		return left.r <= right.r ? left : right;
}


Circle makeDiameter(const Point &a, const Point &b) {
	Point c{(a.x + b.x) / 2, (a.y + b.y) / 2};
	return Circle{c, std::max(c.distance(a), c.distance(b))};
}


Circle makeCircumcircle(const Point &a, const Point &b, const Point &c) {
	// Mathematical algorithm from Wikipedia: Circumscribed circle
	double ox = (std::min(std::min(a.x, b.x), c.x) + std::max(std::min(a.x, b.x), c.x)) / 2;
	double oy = (std::min(std::min(a.y, b.y), c.y) + std::max(std::min(a.y, b.y), c.y)) / 2;
	double ax = a.x - ox, ay = a.y - oy;
	double bx = b.x - ox, by = b.y - oy;
	double cx = c.x - ox, cy = c.y - oy;
	double d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2;
	if (d == 0)
		return Circle::INVALID;
	double x = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (ay - by)) / d;
	double y = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (bx - ax)) / d;
	Point p{ox + x, oy + y};
	double r = std::max(std::max(p.distance(a), p.distance(b)), p.distance(c));
	return Circle{p, r};
}



//------------------------

/*
 input: an object set
 output: the smallest circle that cover all object in the object set
 */
disk_t* MCC(obj_set_t* obj_set_v)
{
    
    obj_node_t* obj_node_v;
    disk_t* cir;
    
//    if(obj_set_v->obj_n == 1)
//    {
//        cir = alloc_disk(2);
//        cir->center->coord[0] = obj_set_v->head->next->obj_v->MBR[0].min;
//        cir->center->coord[1] = obj_set_v->head->next->obj_v->MBR[1].min;
//        cir->radius = 0.0;
//        return cir;
//    }
    
    vector<Point> points;

    
    obj_node_v = obj_set_v->head->next;
    while( obj_node_v != NULL)
    {
        points.push_back(Point{obj_node_v->obj_v->MBR[0].min, obj_node_v->obj_v->MBR[1].min});
        
        obj_node_v = obj_node_v->next;
    }
    
    Circle a = makeSmallestEnclosingCircle(points);
    
    cir = alloc_disk(2);
    cir->center->coord[0] = a.c.x;
    cir->center->coord[1] = a.c.y;
    cir->radius = a.r;
    
    vector<Point>().swap(points);
    

    return cir;
}




