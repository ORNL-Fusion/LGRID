//
//  vertex.cpp
//  lgrid
//
//  Created by Cianciosa, Mark R. on 10/12/15.
//  Copyright (c) 2015 Cianciosa, Mark R. All rights reserved.
//

#include "vertex.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

//------------------------------------------------------------------------------
///  @brief @ref vertex constructor.
///
///  Constructs a single vertex at a point.
///
///  @param[in] point_ Current point to the vertex.
//------------------------------------------------------------------------------
vertex::vertex(const vector3d point_) : point(point_) {
    this->next = this;
    this->last = this;
}

//------------------------------------------------------------------------------
///  @brief @ref vertex destructor.
///
///  Destructs the current vertext and calls the destructor on the next vertex
///  if it exists.
//------------------------------------------------------------------------------
vertex::~vertex() {
    this->last->next = nullptr;
    if (this->next != nullptr) {
        delete this->next;
    }
}

//------------------------------------------------------------------------------
///  @brief Insert a new @ref vertex.
///
///  Inserts a new vertex between the current instance the next instance.
//------------------------------------------------------------------------------
void vertex::insert(vertex *v) {
    this->last->next = v;
    v->last = this->last;
    this->last = v;
    v->next = this;
}

//------------------------------------------------------------------------------
///  @brief Normal to a line.
///
///  Gets the normal to a line defined by two verices.
///
///  @param[in] v1 First vertex.
///  @param[in] v2 Second vertex.
///  @return The normal at a face.
//------------------------------------------------------------------------------
const vector3d vertex::get_normal(const vertex &v1, const vertex &v2) {
    const vector3d line_vector = v2.point - v1.point;
    return vector3d(line_vector.y, -line_vector.x, 0.0).normalize();
}

//------------------------------------------------------------------------------
///  @brief Normal to a line.
///
///  Normal to a line defined by the intersection of two lines. This is the
///  average of the normals of the two lines.
///
///  @param[in] v1 First vertex.
///  @param[in] v2 Second vertex.
///  @param[in] v3 Third vertex.
///  @return The normal at a corner.
//------------------------------------------------------------------------------
const vector3d vertex::get_normal(const vertex &v1, const vertex &v2, const vertex &v3) {
    const vector3d norm_vector1 = get_normal(v1, v2);
    const vector3d norm_vector2 = get_normal(v2, v3);
    return (norm_vector2 + norm_vector1)/2.0;
}

// The next set of functions defines are used to find the distance from a point to a line.
//------------------------------------------------------------------------------
///  @brief a coefficient.
///
///  Coefficient to find used to find the distance from a point to a line.
///
///  @param[in] v1 First vertex.
///  @param[in] v2 Second vertex.
///  @return The a coefficient
//------------------------------------------------------------------------------
const double vertex::a(const vertex &v1, const vertex &v2) {
    return v2.point.y - v1.point.y;
}

//------------------------------------------------------------------------------
///  @brief b coefficient.
///
///  Coefficient to find used to find the distance from a point to a line.
///
///  @param[in] v1 First vertex.
///  @param[in] v2 Second vertex.
///  @return The b coefficient
//------------------------------------------------------------------------------
const double vertex::b(const vertex &v1, const vertex &v2) {
    return -(v2.point.x - v1.point.x);
}

//------------------------------------------------------------------------------
///  @brief c coefficient.
///
///  Coefficient to find used to find the distance from a point to a line.
///
///  @param[in] v1 First vertex.
///  @param[in] v2 Second vertex.
///  @return The c coefficient
//------------------------------------------------------------------------------
const double vertex::c(const vertex &v1, const vertex &v2) {
    return v2.point.x*v1.point.y - v1.point.x*v2.point.y;
}

// Find the minimum x and y position along the line to the point.
//------------------------------------------------------------------------------
///  @brief Minimum distance in the x direction.
///
///  The minimum x along the line to the point.
///
///  @param[in] v1 First vertex.
///  @param[in] v2 Second vertex.
///  @param[in] p  Point of interest.
///  @return The minimum x along the line to the point.
//------------------------------------------------------------------------------
const double vertex::x_min(const vertex &v1, const vertex &v2, const vector3d &p) {
    const double a_value = a(v1, v2);
    const double b_value = b(v1, v2);
    return (b_value*(b_value*p.x - a_value*p.y) - a_value*c(v1, v2))/(a_value*a_value + b_value*b_value);
}

//------------------------------------------------------------------------------
///  @brief Minimum distance in the y direction.
///
///  The minimum y along the line to the point.
///
///  @param[in] v1 First vertex.
///  @param[in] v2 Second vertex.
///  @param[in] p  Point of interest.
///  @return The minimum y along the line to the point.
//------------------------------------------------------------------------------
const double vertex::y_min(const vertex &v1, const vertex &v2, const vector3d &p) {
    const double a_value = a(v1, v2);
    const double b_value = b(v1, v2);
    return (a_value*(-b_value*p.x + a_value*p.y) - b_value*c(v1, v2))/(a_value*a_value + b_value*b_value);
}

//------------------------------------------------------------------------------
///  @brief Distance to a limiter face.
///
///  The length of the line from a point to the limiter face.
///
///  @param[in] v1 First vertex.
///  @param[in] v2 Second vertex.
///  @param[in] p  Point of interest.
///  @return The length to a from a point to the line.
//------------------------------------------------------------------------------
const double vertex::length(const vertex &v1, const vertex &v2, const vector3d &p) {
    const double a_value = a(v1, v2);
    const double b_value = b(v1, v2);
    return fabs(a_value*p.x + b_value*p.y + c(v1, v2))/sqrt(a_value*a_value + b_value*b_value);
}

//------------------------------------------------------------------------------
///  @brief Distance to a vertex.
///
///  The length of a line from a point to the vertex.
///
///  @param[in] v First vertex.
///  @param[in] p  Point of interest.
///  @return The length to a from a point to the line.
//------------------------------------------------------------------------------
const double vertex::length(const vertex &v, const vector3d &p) {
    return (v.point - p).length();
}

//------------------------------------------------------------------------------
///  @brief Finds if the line to the point is in range.
///
///  A point is in range when the line normal to the limiter face that
///  intersects the point is between the two verticies.
///
///  @param[in] v1 First vertex.
///  @param[in] v2 Second vertex.
///  @param[in] x  X position of the point.
///  @param[in] y  Y position of the point.
///  @return The length to a from a point to the line.
//------------------------------------------------------------------------------
const bool vertex::in_range(const vertex &v1, const vertex &v2, const double x, const double y) {
    return (std::min(v1.point.x, v2.point.x) <= x && std::max(v1.point.x, v2.point.x) >= x) && (std::min(v1.point.y, v2.point.y) <= y && std::max(v1.point.y, v2.point.y) >= y);
}

//------------------------------------------------------------------------------
///  @brief Finds which side of the limiter face the point is on.
///
///  A point is in range when the line normal to the limiter face that
///  intersects the point is between the two verticies.
///
///  @param[in] v First vertex.
///  @param[in] n Vector normal to the vertex.
///  @param[in] p Point of interest.
///  @return The length to a from a point to the line.
//------------------------------------------------------------------------------
const double vertex::direction(const vertex &v, const vector3d &n, const vector3d &p) {
    return std::signbit(vector3d::dot_product(v.point - p, n)) ? 1.0 : -1.0;
}

//------------------------------------------------------------------------------
///  @brief Find the distances to the limiter.
///
///  Find the minimum distance to reach in range limiter face and joint between
///  limiter faces.
///
///  @param[in]  start Inital vertex.
///  @param[in]  p     Point of interest.
///  @param[out] d     Minimum distances to the limiter face.
//------------------------------------------------------------------------------
void vertex::distance(vertex *start, const vector3d &p, std::vector<double> &d) {
    for (vertex *v = start; v->next != start; v = v->next) {
        const double x = x_min(*v, *(v->next), p);
        const double y = y_min(*v, *(v->next), p);
        if (in_range(*v, *(v->next), x, y)) {
            d.push_back(direction(*v, get_normal(*v, *(v->next)), p)*length(*v, *(v->next), p));
        } else {
            // Normal is outside the line segment from v to v->next. Check the two corners.
            d.push_back(direction(*v, get_normal(*(v->last), *v, *(v->next)), p)*length(*v, p));
            d.push_back(direction(*(v->next), get_normal(*v, *(v->next), *(v->next->next)), p)*length(*(v->next), p));
        }
    }
    
    // Get distance to last segment.
    const double x = x_min(*(start->last), *start, p);
    const double y = y_min(*(start->last), *start, p);
    if (in_range(*(start->last), *start, x, y)) {
        d.push_back(direction(*(start->last), get_normal(*(start->last), *start), p)*length(*(start->last), *start, p));
    } else {
        // Normal is outside the line segment from v to v->next. Check the two corners.
        d.push_back(direction(*(start->last), get_normal(*(start->last->last), *(start->last), *start), p)*length(*(start->last), p));
        d.push_back(direction(*start, get_normal(*(start->last), *start, *(start->next)), p)*length(*start, p));
    }
}

//------------------------------------------------------------------------------
///  @brief Print out all vertices.
///
///  Prints out information on each vertex.
///
///  @param[in] start Inital vertex.
//------------------------------------------------------------------------------
void vertex::print(vertex *start) {
    std::cout << start->point.x << " " << start->point.y << std::endl;
    for (vertex *v = start->next; v != start; v = v->next) {
        std::cout << v->point.x << " " << v->point.y << std::endl;
    }
}
